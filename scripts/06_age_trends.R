#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 6 - how gene expression changes with age, per sex and per
# developmental stage.
#
#   Rscript scripts/06_age_trends.R
#
# Two complementary views of the same question:
#   A. Spearman rho against age, then hierarchical clustering of the z-scored
#      age profiles into trajectory shapes (age_trend_clusters).
#   B. Pairwise Wilcoxon between ages with a permutation-calibrated |log2FC|
#      threshold, on RNA and on ADT (age_changing), for every stage.
#
# Ages below gates.min_cells_per_stratum are dropped per stratum rather than
# retiring the stratum; the ages each table used are in its `ages` column.
#
# Reads:  results/objects/gmp_neutrophils.rds
# Writes: results/tables/age_trend_<stage>_<sex>.csv
#         results/tables/age_changing_<stage>_<sex>_<assay>.csv
#         results/tables/age_changing_<stage>_shared_<assay>.csv
#         results/figures/pseudobulk_pca_by_stage.pdf, rho_distributions.pdf
# ---------------------------------------------------------------------------

# Run from the project root ("Rscript scripts/06_age_trends.R") or from inside
# scripts/; setup.R locates the project root from either.
source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()

age_levels <- cfg$analysis$age_levels
sexes <- cfg$analysis$sex_levels
gmp_neu <- read_object(cfg, "gmp_neutrophils.rds")
trend_assay <- matched_assay_for(gmp_neu, cfg, step = 6)
require_metadata(gmp_neu, c("stage", "age", "sex"), context = "step 6")

# --- pseudobulk overview ---------------------------------------------------
save_figure(plot_pseudobulk_pca_grid(gmp_neu, cfg), cfg,
            "pseudobulk_pca_by_stage.pdf", width = 12, height = 8)

#' One stage, one sex, restricted to the ages that clear the cell gate.
#'
#' Ages below gates.min_cells_per_stratum are dropped from the stratum rather
#' than retiring it: 9m is thin in both libraries, and requiring every age to
#' clear the gate skipped nine of ten strata. `min_ages` is how many ages the
#' analysis needs -- three for a trend, two for a pairwise contrast.
#'
#' The ages kept are attached as attr(, "ages") so callers pass the same set
#' to the analysis.
stage_subset <- function(obj, stage, sex, min_ages = 3) {
  keep <- as.character(obj$stage) == stage & obj$sex == sex &
    obj$age %in% age_levels
  keep <- !is.na(keep) & keep
  ages <- usable_ages(obj$age[keep], age_levels, cfg$gates$min_cells_per_stratum)
  label <- paste(stage, sex, sep = "/")
  dropped <- setdiff(unique(as.character(obj$age[keep])), ages)
  if (length(ages) < min_ages) {
    log_step(sprintf("  skipping %s: %d age(s) clear the %d-cell gate (%s); need %d",
                     label, length(ages), cfg$gates$min_cells_per_stratum,
                     paste(ages, collapse = ", "), min_ages))
    return(NULL)
  }
  if (length(dropped))
    log_step("  ", label, ": dropping ", paste(dropped, collapse = ", "),
             " (below the cell gate); using ", paste(ages, collapse = ", "))

  cells <- colnames(obj)[keep & obj$age %in% ages]
  sub <- subset(drop_graphs(obj), cells = cells)
  attr(sub, "ages") <- ages
  sub
}

# --- A. Spearman age trends and shape clusters ----------------------------
trend_results <- list()
excess <- list()

for (stage in cfg$analysis$stage_levels) {
  for (sex in sexes) {
    sub <- stage_subset(gmp_neu, stage, sex, min_ages = 3)
    if (is.null(sub)) next
    ages <- attr(sub, "ages")
    log_step("=== age trends: ", stage, " / ", sex, " ===")

    # Excess over a permutation null rather than a bare |rho| cutoff: a fixed
    # threshold has no null and is n-dependent, so its gene count mostly
    # reports how many cells the stratum had.
    res <- age_trend_excess(sub, cfg, age_levels = ages,
                            label = paste(stage, sex, sep = "/"),
                            assay = trend_assay)
    if (is.null(res)) next
    trend_results[[paste(stage, sex, sep = "_")]] <- res

    write_table(strip_inferential_columns(res$changing_genes), cfg,
                sprintf("age_trend_%s_%s.csv", stage, sex))
    excess[[paste(stage, sex, sep = "_")]] <- data.frame(
      stage = stage, sex = sex,
      ages = paste(ages, collapse = ","),
      n_cells_min = min(res$n_cells),
      null_threshold = res$null_threshold,
      n_exceeding_null = res$n_exceeding,
      # Read these two together. A large gene count next to a detection trend
      # above the threshold is one observation, not two.
      detection_rho = res$detection_rho,
      detection_confounded = res$detection_confounded)
  }
}

if (length(excess)) {
  summary_tbl <- do.call(rbind, excess)
  write_table(summary_tbl, cfg, "age_trend_excess_over_null.csv")
  log_step("genes exceeding the permutation null, by stratum:")
  print(summary_tbl, row.names = FALSE)
}

# --- B. Permutation-calibrated pairwise DE --------------------------------
# Every stage, both modalities. A stratum needs two ages that clear the gate;
# the protein (ADT) half is the depth-robust readout, so it matters most in
# exactly the strata where RNA depth is lowest.
de_assays <- c(RNA = trend_assay, ADT = "ADT")
for (stage in cfg$analysis$stage_levels) {
  for (modality in names(de_assays)) {
    assay <- de_assays[[modality]]
    per_sex <- list()
    for (sex in sexes) {
      sub <- stage_subset(gmp_neu, stage, sex, min_ages = 2)
      if (is.null(sub)) next
      res <- tryCatch(
        age_changing(sub, sex, cfg, age_levels = attr(sub, "ages"), assay = assay),
        error = function(e) {
          log_step("  ", stage, "/", sex, "/", modality, " failed: ", conditionMessage(e))
          NULL
        })
      if (is.null(res)) next
      per_sex[[sex]] <- res
      write_table(res, cfg, sprintf("age_changing_%s_%s_%s.csv", stage, sex, modality))
      print(table(res$shape, useNA = "ifany"))
    }

    # Compared on shape only where both sexes kept the same ages; otherwise the
    # shape labels have different lengths and "same shape" means nothing.
    if (length(per_sex) == 2) {
      shared <- compare_sexes(per_sex$male, per_sex$female)
      shared$same_ages <- per_sex$male$ages[1] == per_sex$female$ages[1]
      write_table(shared, cfg, sprintf("age_changing_%s_shared_%s.csv", stage, modality))
      log_step(sprintf("%s / %s -- male: %d | female: %d | shared: %d", stage, modality,
                       nrow(per_sex$male), nrow(per_sex$female), nrow(shared)))
    }
  }
}

# --- Are the age-trend effect sizes systematically larger in one sex? -----
# age_trend_excess() returns the per-gene table as `trend`; this read
# `trend_df`, which never exists, so the comparison below never ran.
mature_trends <- lapply(stats::setNames(sexes, sexes), function(sex)
  trend_results[[paste("mature", sex, sep = "_")]]$trend)

if (!any(vapply(mature_trends, is.null, logical(1)))) {
  clean <- lapply(mature_trends, stats::na.omit)
  save_figure(plot_rho_distributions(clean$male, clean$female, cfg$age_trend$rho_cutoff),
              cfg, "rho_distributions.pdf", width = 7, height = 5)

  for (sex in sexes)
    write_table(clean[[sex]], cfg, sprintf("mature_%s_age_rho.csv", sex))

  print(summary(abs(clean$male$rho)))
  print(summary(abs(clean$female$rho)))
  print(stats::wilcox.test(abs(clean$male$rho), abs(clean$female$rho),
                           alternative = "greater"))
}

log_step("step 6 complete")
