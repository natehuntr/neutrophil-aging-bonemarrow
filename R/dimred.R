# ---------------------------------------------------------------------------
# Normalisation, dimensionality reduction and clustering.
#
# The same four recipes (RNA, SCT, ADT, WNN) are applied to every object in
# the pipeline, so they live here once rather than being retyped per subset.
#
# Small strata are the reason for the guards below: RunPCA defaults to 50
# components, which errors when a subset has fewer cells (or, for ADT, fewer
# features) than that, and RunUMAP errors when asked for more dimensions than
# the PCA produced.
# ---------------------------------------------------------------------------

#' Components a PCA can actually return for this object.
safe_npcs <- function(n_cells, n_features, requested = 50) {
  max(2, min(requested, n_cells - 1, n_features - 1))
}

#' Dimensions available from a reduction, capped at what was asked for.
safe_dims <- function(obj, reduction, requested) {
  available <- ncol(Seurat::Embeddings(obj, reduction = reduction))
  if (requested > available)
    log_step(sprintf("  %s has %d dimensions; using those instead of %d",
                     reduction, available, requested))
  seq_len(min(requested, available))
}

#' Join the split layers every v5 assay gets after a merge.
#'
#' JoinLayers() only touches one assay, and ScaleData on an assay still split
#' into counts.1/counts.2 does not do what it looks like it does.
join_all_layers <- function(obj) {
  original <- Seurat::DefaultAssay(obj)
  for (assay in Seurat::Assays(obj)) {
    if (!inherits(obj[[assay]], "Assay5")) next
    Seurat::DefaultAssay(obj) <- assay
    obj <- Seurat::JoinLayers(obj)
  }
  Seurat::DefaultAssay(obj) <- original
  obj
}

#' Log-normalise, scale, PCA and UMAP on the RNA assay.
#'
#' Reductions are named "pca"/"umap" to match Seurat's defaults.
run_rna_reduction <- function(obj, dims, verbose = FALSE) {
  Seurat::DefaultAssay(obj) <- "RNA"
  obj <- Seurat::NormalizeData(obj, verbose = verbose)
  obj <- Seurat::FindVariableFeatures(obj, verbose = verbose)
  obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(obj), verbose = verbose)
  obj <- Seurat::RunPCA(obj, assay = "RNA", verbose = verbose,
                        npcs = safe_npcs(ncol(obj), length(Seurat::VariableFeatures(obj))))
  Seurat::RunUMAP(obj, assay = "RNA", reduction = "pca",
                  dims = safe_dims(obj, "pca", dims),
                  reduction.name = "umap", verbose = verbose)
}

#' SCTransform v2 normalisation with a matching PCA and UMAP.
#'
#' Used for visualisation and clustering; differential expression stays on the
#' log-normalised RNA assay.
run_sct_reduction <- function(obj, dims, verbose = FALSE) {
  Seurat::DefaultAssay(obj) <- "RNA"
  obj <- Seurat::SCTransform(obj, vst.flavor = "v2", conserve.memory = TRUE,
                             method = "glmGamPoi", verbose = verbose)
  obj <- Seurat::RunPCA(obj, assay = "SCT", reduction.name = "sct_pca", verbose = verbose,
                        npcs = safe_npcs(ncol(obj), nrow(obj[["SCT"]])))
  Seurat::RunUMAP(obj, assay = "SCT", reduction = "sct_pca",
                  dims = safe_dims(obj, "sct_pca", dims),
                  reduction.name = "sct_umap", verbose = verbose)
}

#' Scale the (already isotype-centred) ADT data and reduce it.
#'
#' All antibodies are used as features: the panel is small enough that
#' variable-feature selection would only add noise.
run_adt_reduction <- function(obj, dims, tsne = TRUE, perplexity = 50, verbose = FALSE) {
  Seurat::DefaultAssay(obj) <- "ADT"
  obj <- Seurat::ScaleData(obj, assay = "ADT", verbose = verbose)
  # The panel is typically ~30 antibodies, so the component count is bounded
  # by features rather than by cells.
  obj <- Seurat::RunPCA(obj, assay = "ADT", reduction.name = "adt_pca",
                        features = rownames(obj[["ADT"]]), verbose = verbose,
                        npcs = safe_npcs(ncol(obj), nrow(obj[["ADT"]])))
  adt_dims <- safe_dims(obj, "adt_pca", dims)
  obj <- Seurat::RunUMAP(obj, reduction = "adt_pca", dims = adt_dims,
                         reduction.name = "adt_umap", verbose = verbose)
  if (tsne) {
    # t-SNE requires perplexity < (n_cells - 1) / 3.
    perplexity <- min(perplexity, floor((ncol(obj) - 2) / 3))
    obj <- Seurat::RunTSNE(obj, reduction = "adt_pca", dims = adt_dims,
                           reduction.name = "adt_tsne", perplexity = perplexity)
  }
  Seurat::DefaultAssay(obj) <- "RNA"
  obj
}

#' Weighted-nearest-neighbour graph over the RNA and ADT PCAs.
run_wnn <- function(obj, rna_dims, adt_dims, verbose = FALSE) {
  Seurat::FindMultiModalNeighbors(
    obj,
    reduction.list = list("pca", "adt_pca"),
    dims.list = list(safe_dims(obj, "pca", rna_dims),
                     safe_dims(obj, "adt_pca", adt_dims)),
    verbose = verbose
  )
}

#' Cluster across a grid of resolutions so clustree can be used to pick one.
#'
#' Results land in metadata columns clust_<resolution>; the chosen resolution
#' is copied to `clusters` by `assign_clusters()`.
sweep_resolutions <- function(obj, resolutions, graph.name = "wsnn", algorithm = 4) {
  # algorithm 4 is Leiden, which needs the Python leidenalg package via
  # reticulate. Set clustering.algorithm to 1 (Louvain) in the config if that
  # is not installed.
  for (r in resolutions) {
    obj <- Seurat::FindClusters(obj, resolution = r, algorithm = algorithm,
                                graph.name = graph.name,
                                cluster.name = paste0("clust_", r),
                                verbose = FALSE)
  }
  obj
}

assign_clusters <- function(obj, resolution) {
  column <- paste0("clust_", resolution)
  if (!column %in% colnames(obj@meta.data))
    stop("no clustering at resolution ", resolution,
         "; run sweep_resolutions() with that value in the grid")
  obj$clusters <- obj[[column]][, 1]
  obj
}
