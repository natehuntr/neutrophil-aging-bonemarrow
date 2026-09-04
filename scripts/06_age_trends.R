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
#      threshold, on RNA and on ADT (age_changing).
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
require_metadata(gmp_neu, c("fine_neu_labels", "age", "sex"), context = "step 6")

# --- pseudobulk overview ---------------------------------------------------
save_figure(plot_pseudobulk_pca_grid(gmp_neu, cfg), cfg,
            "pseudobulk_pca_by_stage.pdf", width = 12, height = 8)

#' One stage, one sex: re-normalised so variable features and scaling reflect
#' that stratum rather than the whole object.
stage_subset <- function(obj, stage, sex) {
  keep <- as.character(obj$fine_neu_labels) == stage &
    obj$sex == sex &
    obj$age %in% age_levels
  keep <- !is.na(keep) & keep
  cells <- colnames(obj)[keep]
  if (length(cells) < 30) return(NULL)

  sub <- subset(obj, cells = cells)
  Seurat::DefaultAssay(sub) <- "RNA"
  sub <- Seurat::NormalizeData(sub, verbose = FALSE)
  sub <- Seurat::FindVariableFeatures(sub, verbose = FALSE)
  Seurat::ScaleData(sub, features = Seurat::VariableFeatures(sub), verbose = FALSE)
}

# --- A. Spearman age trends and shape clusters ----------------------------
trend_results <- list()

for (stage in cfg$analysis$stage_levels) {
  for (sex in sexes) {
    sub <- stage_subset(gmp_neu, stage, sex)
    if (is.null(sub)) {
      log_step("skipping ", stage, " / ", sex, ": too few cells")
      next
    }
    log_step("=== age trends: ", stage, " / ", sex, " ===")

    res <- age_trend_clusters(sub, cfg, plot = FALSE)
    trend_results[[paste(stage, sex, sep = "_")]] <- res

    write_table(age_trend_table(res), cfg,
                sprintf("age_trend_%s_%s.csv", stage, sex))
  }
}

# --- B. Permutation-calibrated pairwise DE --------------------------------
# Run on the two stages with enough cells in both sexes for the permutation
# null to mean anything.
for (stage in c("mature", "immature")) {
  changing <- list()

  for (assay in c("RNA", "ADT")) {
    per_sex <- list()
    for (sex in sexes) {
      sub <- stage_subset(gmp_neu, stage, sex)
      if (is.null(sub)) next
      per_sex[[sex]] <- age_changing(sub, sex, cfg, age_levels, assay = assay)
      write_table(per_sex[[sex]], cfg,
                  sprintf("age_changing_%s_%s_%s.csv", stage, sex, assay))
      print(table(per_sex[[sex]]$shape, useNA = "ifany"))
    }

    # The shared set is only interpretable if both arms passed the sanity
    # check inside age_changing(); a shared list built against a broken arm
    # is meaningless.
    if (length(per_sex) == 2) {
      shared <- compare_sexes(per_sex$male, per_sex$female)
      write_table(shared, cfg, sprintf("age_changing_%s_shared_%s.csv", stage, assay))
      log_step(sprintf("%s / %s -- male: %d | female: %d | shared: %d", stage, assay,
                       nrow(per_sex$male), nrow(per_sex$female), nrow(shared)))
    }
    changing[[assay]] <- per_sex
  }
}

# --- Are the age-trend effect sizes systematically larger in one sex? -----
mature_trends <- lapply(stats::setNames(sexes, sexes), function(sex)
  trend_results[[paste("mature", sex, sep = "_")]]$trend_df)

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
