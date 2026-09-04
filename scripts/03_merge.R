#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 3 - merge the two sexes and build the granulopoiesis object.
#
#   Rscript scripts/03_merge.R
#
# Two objects come out of this step:
#   bm_merged.rds        every cell from both sexes, jointly reduced
#   gmp_neutrophils.rds  GMPs + neutrophils only, module-scored, clustered
#                        into developmental stages. This is the object every
#                        downstream analysis works from.
#
# CytoTRACE2 potency is computed in step 2 on the stem/neutrophil subset only,
# so it is transferred by barcode here; cells outside that subset carry NA.
# ---------------------------------------------------------------------------

# Run from the project root ("Rscript scripts/03_merge.R") or from inside
# scripts/; setup.R locates the project root from either.
source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()

sample_keys <- names(cfg$samples)

annotated <- lapply(stats::setNames(sample_keys, sample_keys),
                    function(k) read_object(cfg, paste0(k, "_annotated.rds")))
neustem   <- lapply(stats::setNames(sample_keys, sample_keys),
                    function(k) read_object(cfg, paste0(k, "_neustem.rds")))

# --- carry the potency calls back onto the full objects --------------------
#' Copy columns from a subset back onto the full object, by cell name.
#'
#' Cells absent from the subset get NA. A requested column that the subset does
#' not have is reported rather than silently skipped: downstream steps select
#' on CytoTRACE2_Potency, and if it never arrives they fail with an empty
#' subset three steps from here.
transfer_metadata <- function(target, source, columns, label = "") {
  absent <- setdiff(columns, colnames(source@meta.data))
  if (length(absent))
    warning(label, ": column(s) not present on the subset and so not ",
            "transferred: ", paste(absent, collapse = ", "))

  idx <- match(colnames(source), colnames(target))
  matched <- !is.na(idx)
  if (!any(matched))
    stop(label, ": no cell names in common between the subset and the full ",
         "object, so nothing can be transferred. Were they built from ",
         "different runs of step 2?")

  for (col in intersect(columns, colnames(source@meta.data))) {
    values <- rep(NA, ncol(target))
    values[idx[matched]] <- source@meta.data[[col]][matched]
    target[[col]] <- values
  }
  target
}

for (key in sample_keys)
  annotated[[key]] <- transfer_metadata(annotated[[key]], neustem[[key]],
                                        CYTOTRACE_COLUMNS, label = key)

# Steps 4 and 7 select differentiated cells on this column, so confirm it
# arrived here rather than discovering it two steps downstream.
for (key in sample_keys) {
  coverage <- metadata_coverage(annotated[[key]], CYTOTRACE_COLUMNS)
  log_step(key, ": CytoTRACE2 coverage after transfer")
  print(coverage, row.names = FALSE)

  potency <- coverage$non_na[coverage$column == "CytoTRACE2_Potency"]
  if (is.na(potency) || potency == 0)
    stop(key, ": no cell has a CytoTRACE2_Potency value after the transfer.\n",
         "Steps 4 and 7 select on it, so they would fail with an empty subset.\n",
         "Re-run step 2 for this sample and check the CytoTRACE2 output.")
}

# --- all cells, both sexes -------------------------------------------------
log_step("=== merging all cells ===")
bm_merged <- merge(annotated[[1]], y = annotated[-1], merge.data = TRUE)
bm_merged <- join_all_layers(bm_merged)
bm_merged <- run_sct_reduction(bm_merged, cfg$samples[[1]]$dims$sct_pca)
bm_merged <- run_rna_reduction(bm_merged, cfg$samples[[1]]$dims$rna_pca)
bm_merged <- run_adt_reduction(bm_merged, cfg$samples[[1]]$dims$neustem_adt_pca)
save_object(bm_merged, cfg, "bm_merged.rds")

# --- GMPs + neutrophils, both sexes ---------------------------------------
log_step("=== merging GMPs and neutrophils ===")
granulo_labels <- c("Granulocyte-Monocyte Progenitor (GMP)", "Neutrophils")

subsets <- lapply(names(neustem), function(key) {
  obj <- neustem[[key]]
  select_cells(obj, list(
    "fine_label_readable is a GMP or neutrophil" =
      as.character(obj$fine_label_readable) %in% granulo_labels
  ), context = paste(key, "GMPs and neutrophils"))
})

gmp_neu <- merge(subsets[[1]], y = subsets[-1], merge.data = TRUE)
gmp_neu <- join_all_layers(gmp_neu)
gmp_neu <- run_rna_reduction(gmp_neu, cfg$samples[[1]]$dims$rna_pca)
gmp_neu <- add_module_scores(gmp_neu)

# Cluster on the RNA graph and name the clusters after the maturation stages.
gmp_neu <- Seurat::FindNeighbors(gmp_neu, dims = seq_len(cfg$samples[[1]]$dims$rna_pca),
                                 verbose = FALSE)
gmp_neu <- Seurat::FindClusters(gmp_neu, resolution = cfg$clustering$chosen_resolution,
                                verbose = FALSE)
gmp_neu <- add_stage_labels(gmp_neu, cfg)

log_step("cells per stage:")
print(table(gmp_neu$fine_neu_labels, gmp_neu$sex, useNA = "ifany"))

# The stage assignment is a manual mapping from cluster ids: these plots are
# how it gets checked, so they are always written.
module_plots <- lapply(paste0(names(NEUTROPHIL_MODULES), "_score"), function(m)
  Seurat::FeaturePlot(gmp_neu, reduction = "umap", features = m))
save_figure(patchwork::wrap_plots(module_plots), cfg, "gmp_neutrophil_modules.pdf")
save_figure(Seurat::DimPlot(gmp_neu, reduction = "umap", group.by = "seurat_clusters") +
              Seurat::DimPlot(gmp_neu, reduction = "umap", group.by = "fine_neu_labels"),
            cfg, "gmp_neutrophil_clusters.pdf", width = 12, height = 5)

save_object(gmp_neu, cfg, "gmp_neutrophils.rds")

log_step("step 3 complete")
