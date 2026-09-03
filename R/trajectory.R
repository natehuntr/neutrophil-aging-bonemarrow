# ---------------------------------------------------------------------------
# Monocle3 trajectories over the GMP -> mature neutrophil axis.
#
# INTERACTIVITY: monocle3::order_cells() and choose_graph_segments() open a
# Shiny window when called with no root. That is fine in RStudio but stalls a
# batch run, so `root_group` lets the root be chosen programmatically from a
# metadata column (the standard monocle3 recipe: pick the principal graph node
# closest to the cells of the earliest developmental stage).
# ---------------------------------------------------------------------------

#' Convert a Seurat object to a monocle3 cell_data_set and preprocess it.
to_cds <- function(obj, assay = "RNA", num_dim = 50) {
  obj <- Seurat::JoinLayers(obj)
  cds <- SeuratWrappers::as.cell_data_set(obj, assay = assay)
  cds <- monocle3::preprocess_cds(cds, num_dim = num_dim)
  monocle3::cluster_cells(cds)
}

#' The principal-graph node that most of `root_group`'s cells sit closest to.
pick_root_nodes <- function(cds, group_col, root_group) {
  groups <- SummarizedExperiment::colData(cds)[[group_col]]
  if (is.null(groups))
    stop("column '", group_col, "' is not present in colData(cds)")
  cell_ids <- which(as.character(groups) == root_group)
  if (!length(cell_ids))
    stop("no cells with ", group_col, " == '", root_group, "' to root the trajectory on")

  vertex_by_cell <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  vertex_by_cell <- as.matrix(vertex_by_cell)[colnames(cds), ]
  counts <- table(vertex_by_cell[cell_ids])
  igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(counts)))]
}

#' Learn a principal graph and order cells along it.
#'
#' @param root_group value of `root_group_col` marking the start of the
#'   trajectory. If NULL, order_cells() is called interactively.
learn_trajectory <- function(cds,
                             use_partition = FALSE,
                             ncenter = 300,
                             minimal_branch_len = NULL,
                             root_group = NULL,
                             root_group_col = "fine_neu_labels") {
  control <- list(prune_graph = TRUE)
  if (!is.null(ncenter)) control$ncenter <- ncenter
  if (!is.null(minimal_branch_len)) control$minimal_branch_len <- minimal_branch_len

  cds <- monocle3::learn_graph(cds,
                               use_partition = use_partition,
                               close_loop = FALSE,
                               learn_graph_control = control,
                               verbose = FALSE)

  if (is.null(root_group)) {
    if (!interactive())
      stop("root_group is NULL and the session is not interactive: ",
           "pass root_group (e.g. \"GMPs\") so the root can be chosen without a GUI")
    monocle3::order_cells(cds, reduction_method = "UMAP")
  } else {
    monocle3::order_cells(cds, reduction_method = "UMAP",
                          root_pr_nodes = pick_root_nodes(cds, root_group_col, root_group))
  }
}

#' Genes varying across the principal graph (Moran's I), most spatially
#' structured first.
trajectory_genes <- function(cds, cores = 4) {
  res <- monocle3::graph_test(cds, neighbor_graph = "principal_graph",
                              reduction_method = "UMAP", cores = cores)
  res[rev(order(res$morans_test_statistic)), ]
}

#' Spearman correlation of every gene against pseudotime.
#'
#' Cells with non-finite pseudotime (disconnected from the root) are dropped.
#' Uses the vectorised spearman_rows() from R/de.R.
pseudotime_spearman <- function(cds) {
  pt <- monocle3::pseudotime(cds)
  expr <- SingleCellExperiment::counts(cds)[, names(pt), drop = FALSE]

  valid <- is.finite(pt)
  expr <- expr[, valid, drop = FALSE]
  pt <- pt[valid]

  df <- spearman_rows(expr, pt)
  names(df)[names(df) == "pval"] <- "p_value"
  df$padj <- stats::p.adjust(df$p_value, method = "BH")
  df[order(-abs(df$rho)), ]
}

#' Run the whole per-group trajectory analysis and write both result tables.
#'
#' @param group_cells named list of cell barcodes, one entry per group to
#'   analyse (e.g. one per age).
trajectory_by_group <- function(obj, cfg, group_cells, prefix,
                                root_group = "GMPs",
                                root_group_col = "fine_neu_labels",
                                use_partition = FALSE,
                                ncenter = 300,
                                minimal_branch_len = 10,
                                cores = 4) {
  results <- list()
  for (nm in names(group_cells)) {
    log_step("trajectory: ", prefix, " ", nm)
    cds <- to_cds(subset(obj, cells = group_cells[[nm]]))
    cds <- learn_trajectory(cds,
                            use_partition = use_partition,
                            ncenter = ncenter,
                            minimal_branch_len = minimal_branch_len,
                            root_group = root_group,
                            root_group_col = root_group_col)

    moran <- trajectory_genes(cds, cores = cores)
    spearman <- pseudotime_spearman(cds)

    write_table(moran, cfg, sprintf("%s_%s_trajectory_moran.csv", prefix, nm),
                row.names = TRUE)
    write_table(spearman, cfg, sprintf("%s_%s_pseudotime_spearman.csv", prefix, nm))

    results[[nm]] <- list(cds = cds, moran = moran, spearman = spearman)
  }
  results
}

#' Genes significant in one group and in none of the others.
genes_unique_to <- function(gene_lists, target) {
  setdiff(gene_lists[[target]], unlist(gene_lists[names(gene_lists) != target]))
}
