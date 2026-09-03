# ---------------------------------------------------------------------------
# Normalisation, dimensionality reduction and clustering.
#
# The same four recipes (RNA, SCT, ADT, WNN) are applied to every object in
# the pipeline, so they live here once rather than being retyped per subset.
# ---------------------------------------------------------------------------

#' Log-normalise, scale, PCA and UMAP on the RNA assay.
#'
#' Reductions are named "pca"/"umap" to match Seurat's defaults.
run_rna_reduction <- function(obj, dims, verbose = FALSE) {
  Seurat::DefaultAssay(obj) <- "RNA"
  obj <- Seurat::NormalizeData(obj, verbose = verbose)
  obj <- Seurat::FindVariableFeatures(obj, verbose = verbose)
  obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(obj), verbose = verbose)
  obj <- Seurat::RunPCA(obj, assay = "RNA", verbose = verbose)
  Seurat::RunUMAP(obj, assay = "RNA", reduction = "pca",
                  dims = seq_len(dims), reduction.name = "umap", verbose = verbose)
}

#' SCTransform v2 normalisation with a matching PCA and UMAP.
#'
#' Used for visualisation and clustering; differential expression stays on the
#' log-normalised RNA assay.
run_sct_reduction <- function(obj, dims, verbose = FALSE) {
  Seurat::DefaultAssay(obj) <- "RNA"
  obj <- Seurat::SCTransform(obj, vst.flavor = "v2", conserve.memory = TRUE,
                             method = "glmGamPoi", verbose = verbose)
  obj <- Seurat::RunPCA(obj, assay = "SCT", reduction.name = "sct_pca", verbose = verbose)
  Seurat::RunUMAP(obj, assay = "SCT", reduction = "sct_pca",
                  dims = seq_len(dims), reduction.name = "sct_umap", verbose = verbose)
}

#' Scale the (already isotype-centred) ADT data and reduce it.
#'
#' All antibodies are used as features: the panel is small enough that
#' variable-feature selection would only add noise.
run_adt_reduction <- function(obj, dims, tsne = TRUE, perplexity = 50, verbose = FALSE) {
  Seurat::DefaultAssay(obj) <- "ADT"
  obj <- Seurat::ScaleData(obj, assay = "ADT", verbose = verbose)
  obj <- Seurat::RunPCA(obj, assay = "ADT", reduction.name = "adt_pca",
                        features = rownames(obj[["ADT"]]), verbose = verbose)
  obj <- Seurat::RunUMAP(obj, reduction = "adt_pca", dims = seq_len(dims),
                         reduction.name = "adt_umap", verbose = verbose)
  if (tsne)
    obj <- Seurat::RunTSNE(obj, reduction = "adt_pca", dims = seq_len(dims),
                           reduction.name = "adt_tsne", perplexity = perplexity)
  Seurat::DefaultAssay(obj) <- "RNA"
  obj
}

#' Weighted-nearest-neighbour graph over the RNA and ADT PCAs.
run_wnn <- function(obj, rna_dims, adt_dims, verbose = FALSE) {
  Seurat::FindMultiModalNeighbors(
    obj,
    reduction.list = list("pca", "adt_pca"),
    dims.list = list(seq_len(rna_dims), seq_len(adt_dims)),
    verbose = verbose
  )
}

#' Cluster across a grid of resolutions so clustree can be used to pick one.
#'
#' Results land in metadata columns clust_<resolution>; the chosen resolution
#' is copied to `clusters` by `assign_clusters()`.
sweep_resolutions <- function(obj, resolutions, graph.name = "wsnn", algorithm = 4) {
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
