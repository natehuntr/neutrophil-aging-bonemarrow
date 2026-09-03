#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 2 - reductions, reference labelling and the stem-cell/neutrophil subset.
#
#   Rscript scripts/02_annotate.R [sample_key ...]
#
# Reads:  results/objects/<key>_filtered.rds
# Writes: results/objects/<key>_annotated.rds   (all cells, labelled)
#         results/objects/<key>_neustem.rds     (stem cells + neutrophils)
# ---------------------------------------------------------------------------

# Run from the project root ("Rscript scripts/02_annotate.R") or from inside
# scripts/; setup.R locates the project root from either.
source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()
require_packages("SingleR", "celldex", "CytoTRACE2", "clustree",
                 "SingleCellExperiment", "SummarizedExperiment")

args <- commandArgs(trailingOnly = TRUE)
sample_keys <- if (length(args)) args else names(cfg$samples)

for (key in sample_keys) {
  sample_cfg <- sample_config(cfg, key)
  dims <- sample_cfg$dims
  log_step("=== annotating ", key, " ===")

  obj <- read_object(cfg, paste0(key, "_filtered.rds"))

  obj <- run_sct_reduction(obj, dims$sct_pca)
  obj <- run_rna_reduction(obj, dims$rna_pca)
  obj <- annotate_singler(obj, cfg)
  obj <- add_readable_labels(obj)

  save_figure(plot_elbows(obj), cfg, paste0(key, "_elbow.pdf"))
  save_object(obj, cfg, paste0(key, "_annotated.rds"))

  # --- stem cells + neutrophils ------------------------------------------
  keep <- grepl(cfg$annotation$neustem_label_pattern, obj$singleR_fine_label)
  neustem <- subset(obj, cells = colnames(obj)[keep])
  log_step(sprintf("stem/neutrophil subset: %d of %d cells", ncol(neustem), ncol(obj)))

  neustem <- run_rna_reduction(neustem, dims$neustem_rna_pca)
  neustem <- run_adt_reduction(neustem, dims$neustem_adt_pca)
  neustem <- run_wnn(neustem, dims$wnn_rna, dims$wnn_adt)

  # Resolution sweep over the WNN graph; clustree plots the stability of the
  # partitions so the chosen resolution is a decision, not a default.
  neustem <- sweep_resolutions(neustem, cfg$clustering$resolution_grid,
                               algorithm = cfg$clustering$algorithm)
  save_figure(clustree::clustree(neustem@meta.data, prefix = "clust_"),
              cfg, paste0(key, "_clustree.pdf"))
  neustem <- assign_clusters(neustem, cfg$clustering$chosen_resolution)

  neustem <- run_cytotrace2(neustem, cfg)

  save_object(neustem, cfg, paste0(key, "_neustem.rds"))
  rm(obj, neustem)
  gc(verbose = FALSE)
}

log_step("step 2 complete")
