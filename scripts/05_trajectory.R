#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 5 - pseudotime trajectories along granulopoiesis.
#
#   Rscript scripts/05_trajectory.R
#
# One monocle3 trajectory per sex x age, plus a combined trajectory used to
# fit per-sex tradeSeq GAMs with age as the condition.
#
# Reads:  results/objects/gmp_neutrophils.rds
# Writes: results/tables/<sex>_<age>_trajectory_moran.csv
#         results/tables/<sex>_<age>_pseudotime_spearman.csv
#         results/objects/combined_cds.rds, <sex>_GAM.rds
#         results/figures/<sex>_trajectory_genes_venn.pdf
#
# The trajectory root is chosen programmatically from the GMP cluster (see
# R/trajectory.R), so this script runs unattended. In an interactive session
# pass root_group = NULL to learn_trajectory() to pick the root by hand.
# ---------------------------------------------------------------------------

# Run from the project root ("Rscript scripts/05_trajectory.R") or from inside
# scripts/; setup.R locates the project root from either.
source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()
require_packages("monocle3", "SeuratWrappers", "tradeSeq", "SingleCellExperiment")

age_levels <- cfg$analysis$age_levels
gmp_neu <- read_object(cfg, "gmp_neutrophils.rds")

# --- 1. One trajectory per sex x age --------------------------------------
moran_lists <- list()

for (sex in cfg$analysis$sex_levels) {
  cells_by_age <- lapply(stats::setNames(age_levels, age_levels), function(age)
    colnames(gmp_neu)[gmp_neu$sex == sex & gmp_neu$age == age])

  too_small <- names(cells_by_age)[lengths(cells_by_age) < 50]
  if (length(too_small))
    warning(sex, ": fewer than 50 cells at ", paste(too_small, collapse = ", "),
            " -- those trajectories are not worth much")

  results <- trajectory_by_group(gmp_neu, cfg, cells_by_age, prefix = sex,
                                 root_group = "GMPs")

  # Significantly graph-associated genes, per age, for the overlap plot.
  sig <- lapply(results, function(r) rownames(r$moran)[r$moran$q_value < 0.05])
  moran_lists[[sex]] <- sig

  save_figure(plot_gene_venn(sig, title = paste0(sex, ": trajectory genes by age")),
              cfg, paste0(sex, "_trajectory_genes_venn.pdf"), width = 6, height = 6)

  for (age in names(sig)) {
    unique_genes <- genes_unique_to(sig, age)
    write_table(data.frame(gene = unique_genes), cfg,
                sprintf("%s_%s_trajectory_genes_unique.csv", sex, age))
  }
}

# --- 2. Combined trajectory across both sexes ------------------------------
log_step("=== combined trajectory ===")
combined_cds <- to_cds(gmp_neu)
combined_cds <- learn_trajectory(combined_cds, use_partition = FALSE, ncenter = 300,
                                 root_group = "GMPs")
save_object(combined_cds, cfg, "combined_cds.rds")

# Genes worth fitting: significant on any of the per-age trajectories.
all_sig_genes <- unique(unlist(moran_lists))
all_sig_genes <- intersect(all_sig_genes, rownames(SingleCellExperiment::counts(combined_cds)))
log_step(length(all_sig_genes), " genes carried into the GAM fits")
if (length(all_sig_genes) < 10)
  stop("only ", length(all_sig_genes), " genes were graph-associated in any ",
       "per-age trajectory; there is nothing to fit. Check the trajectories in ",
       "results/tables/*_trajectory_moran.csv before going further.")

# --- 3. Per-sex tradeSeq GAMs, conditioned on age -------------------------
# Each sex is fitted separately on the shared pseudotime axis; `conditions`
# is age, so the smoothers can be compared across timepoints within a sex.
for (sex in cfg$analysis$sex_levels) {
  log_step("fitting GAM: ", sex)
  cds_sex <- combined_cds[, SummarizedExperiment::colData(combined_cds)$sex == sex]

  pt <- monocle3::pseudotime(cds_sex)
  finite <- is.finite(pt)
  cds_sex <- cds_sex[, finite]
  pt <- pt[finite]

  ages_present <- table(SummarizedExperiment::colData(cds_sex)$age)
  if (any(ages_present < 20))
    warning(sex, ": ", paste(names(ages_present), ages_present, sep = "=",
                             collapse = ", "),
            " cells per age -- the condition smoothers will be unstable")

  gam <- tradeSeq::fitGAM(
    counts = SingleCellExperiment::counts(cds_sex)[all_sig_genes, , drop = FALSE],
    pseudotime = pt,
    cellWeights = rep(1, ncol(cds_sex)),
    nknots = 5,
    conditions = factor(SummarizedExperiment::colData(cds_sex)$age, levels = age_levels)
  )
  save_object(gam, cfg, paste0(sex, "_GAM.rds"))

  condition_res <- tradeSeq::conditionTest(gam, l2fc = 0)
  condition_res$padj <- stats::p.adjust(condition_res$pvalue, "fdr")
  condition_res$gene <- rownames(condition_res)
  write_table(condition_res[order(condition_res$padj), ], cfg,
              paste0(sex, "_gam_condition_test.csv"))
}

log_step("step 5 complete")
