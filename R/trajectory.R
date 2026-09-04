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
  obj <- join_layers(obj)
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
                                cores = cfg$compute$cores %||% 4) {
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

# ---------------------------------------------------------------------------
# Where cells sit along the trajectory.
#
# The gene-level trajectory models ask whether a gene's profile along
# pseudotime changes. They cannot see the other way development can change:
# cells piling up before a maturation step, or traversing a compressed range.
# That is a shift in the *distribution* of cells along pseudotime, and it is
# what the functions below measure.
#
# The same pseudoreplication caveat applies as everywhere else: a KS test over
# cells treats them as independent draws. The distances and median shifts are
# the result; the p-values rank comparisons and nothing more.
# ---------------------------------------------------------------------------

#' 1-Wasserstein (earth mover's) distance between two samples.
#'
#' The mean absolute gap between the two quantile functions. Unlike the KS
#' statistic it is in the units of the variable, so it can be compared across
#' comparisons and is not driven by a single crossing point.
wasserstein1 <- function(x, y, n_grid = 1000) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  probs <- (seq_len(n_grid) - 0.5) / n_grid
  mean(abs(stats::quantile(x, probs, names = FALSE, type = 7) -
           stats::quantile(y, probs, names = FALSE, type = 7)))
}

#' Per-group summary of a distribution.
distribution_summary <- function(values, groups, group_levels) {
  do.call(rbind, lapply(group_levels, function(g) {
    v <- values[groups == g]
    v <- v[is.finite(v)]
    if (!length(v)) return(NULL)
    data.frame(group = g, n = length(v),
               mean = mean(v), median = stats::median(v),
               q25 = unname(stats::quantile(v, 0.25)),
               q75 = unname(stats::quantile(v, 0.75)),
               iqr = unname(diff(stats::quantile(v, c(0.25, 0.75)))),
               stringsAsFactors = FALSE)
  }))
}

#' Pairwise comparison of a distribution across groups.
#'
#' Every group is compared against `reference` (the first level by default),
#' which keeps the comparisons interpretable as "shift away from baseline"
#' rather than an all-pairs matrix that has to be read as a whole.
#'
#' `frac_beyond_reference_median` is the share of the group's cells sitting
#' past the reference group's median. It is 0.5 for the reference by
#' construction, so a value of 0.65 at 18m reads directly as "cells have moved
#' further along the trajectory".
compare_distributions <- function(values, groups, group_levels,
                                  reference = group_levels[1]) {
  groups <- as.character(groups)
  ref_values <- values[groups == reference]
  ref_values <- ref_values[is.finite(ref_values)]
  if (!length(ref_values))
    stop("reference group '", reference, "' has no finite values")
  ref_median <- stats::median(ref_values)

  res <- do.call(rbind, lapply(setdiff(group_levels, reference), function(g) {
    v <- values[groups == g]
    v <- v[is.finite(v)]
    if (length(v) < 10) {
      log_step("  skipping ", g, ": ", length(v), " cells")
      return(NULL)
    }
    ks <- suppressWarnings(stats::ks.test(v, ref_values))
    data.frame(
      group = g, reference = reference, n = length(v), n_reference = length(ref_values),
      median_shift = stats::median(v) - ref_median,
      wasserstein = wasserstein1(v, ref_values),
      frac_beyond_reference_median = mean(v > ref_median),
      ks_statistic = unname(ks$statistic),
      ks_p = ks$p.value,
      stringsAsFactors = FALSE
    )
  }))

  if (!is.null(res)) res$ks_padj <- stats::p.adjust(res$ks_p, method = "BH")
  res
}

#' Compare a per-cell quantity across ages, separately within each sex.
#'
#' Within-sex is the design-clean comparison: all ages share a library, so an
#' age shift cannot be a batch shift. Comparing the same quantity between
#' sexes would compare two libraries and is deliberately not done here.
compare_across_ages_by_sex <- function(values, meta, cfg,
                                       age_levels = cfg$analysis$age_levels,
                                       age_col = "age", sex_col = "sex") {
  summaries <- list()
  shifts <- list()

  for (this_sex in cfg$analysis$sex_levels) {
    keep <- meta[[sex_col]] == this_sex & meta[[age_col]] %in% age_levels
    keep[is.na(keep)] <- FALSE
    if (sum(keep) < 30) {
      log_step("skipping ", this_sex, ": ", sum(keep), " cells")
      next
    }
    log_step("=== ", this_sex, " ===")

    v <- values[keep]
    g <- as.character(meta[[age_col]][keep])

    s <- distribution_summary(v, g, age_levels)
    if (!is.null(s)) { s$sex <- this_sex; summaries[[this_sex]] <- s }

    d <- compare_distributions(v, g, age_levels)
    if (!is.null(d)) { d$sex <- this_sex; shifts[[this_sex]] <- d }
  }

  list(summary = do.call(rbind, summaries),
       shifts = do.call(rbind, shifts))
}
