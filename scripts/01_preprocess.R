#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 1 - per-sample preprocessing.
#
#   Rscript scripts/01_preprocess.R            # every sample in the config
#   Rscript scripts/01_preprocess.R male       # just one
#
# Reads the Cell Ranger filtered matrix, builds the RNA / ADT / HTO assays,
# demultiplexes the hashtags into ages, computes QC metrics, calls doublets
# and applies the QC thresholds.
#
# Writes: results/objects/<key>_filtered.rds
#         results/figures/<key>_qc.pdf, <key>_doublets.pdf
# ---------------------------------------------------------------------------

# Run from the project root ("Rscript scripts/01_preprocess.R") or from
# inside scripts/; setup.R finds the project root from there.
source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()
# SummarizedExperiment and SingleCellExperiment must be ATTACHED, not just
# installed: scDblFinder relies on S4 generics resolved through the search
# path. See the note in R/packages.R.
require_packages("scDblFinder", "SingleCellExperiment", "SummarizedExperiment")

args <- commandArgs(trailingOnly = TRUE)
sample_keys <- if (length(args)) args else names(cfg$samples)

for (key in sample_keys) {
  sample_cfg <- sample_config(cfg, key)
  log_step("=== preprocessing ", key, " (", sample_cfg$sample_id, ") ===")

  mats <- read_cite_sample(cfg, sample_cfg$sample_id)

  obj <- create_rna_object(mats)
  obj <- add_protein_assays(obj, mats, cfg, hashtags = names(sample_cfg$hashtags))
  obj <- demultiplex_hashtags(obj, sample_cfg)
  obj <- add_qc_metrics(obj, cfg)
  obj <- add_doublet_calls(obj)

  save_figure(plot_qc_metrics(obj), cfg, paste0(key, "_qc.pdf"))
  save_figure(plot_doublet_scores(obj), cfg, paste0(key, "_doublets.pdf"))

  obj <- filter_cells(obj, cfg, hashtags = names(sample_cfg$hashtags))
  save_object(obj, cfg, paste0(key, "_filtered.rds"))

  rm(mats, obj)
  gc(verbose = FALSE)
}

log_step("step 1 complete")
