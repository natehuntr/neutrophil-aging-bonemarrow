# ---------------------------------------------------------------------------
# Reading Cell Ranger output.
# ---------------------------------------------------------------------------

#' Locate a sample's raw and filtered matrix directories.
sample_matrix_paths <- function(cfg, sample_id) {
  root <- cfg$paths$cellranger_root
  paths <- list(
    raw      = file.path(root, sample_id, cfg$paths$raw_matrix_dir),
    filtered = file.path(root, sample_id, cfg$paths$filtered_matrix_dir)
  )
  missing <- paths[!vapply(paths, dir.exists, logical(1))]
  if (length(missing))
    stop("matrix directory not found:\n  ", paste(unlist(missing), collapse = "\n  "),
         "\nCheck paths.cellranger_root in config/config.yml.")
  paths
}

#' Read one hashed CITE-seq sample.
#'
#' Returns the four matrices the pipeline needs, with the filtered gene matrix
#' padded out to the full raw gene list. Padding rather than intersecting keeps
#' the two matrices on identical row sets, which SoupChannel requires, without
#' silently dropping genes that Cell Ranger filtered out of the cell matrix.
read_cite_sample <- function(cfg, sample_id) {
  paths <- sample_matrix_paths(cfg, sample_id)
  log_step("reading ", sample_id)

  raw      <- Seurat::Read10X(paths$raw)
  filtered <- Seurat::Read10X(paths$filtered)

  raw_gex      <- raw[["Gene Expression"]]
  filtered_gex <- align_gene_rows(raw_gex, filtered[["Gene Expression"]])

  list(
    sample_id    = sample_id,
    raw_gex      = raw_gex,
    filtered_gex = filtered_gex,
    raw_adt      = raw[["Antibody Capture"]],
    filtered_adt = filtered[["Antibody Capture"]]
  )
}

#' Zero-pad `target` so its rows match `reference`, then reorder to match.
align_gene_rows <- function(reference, target) {
  missing_genes <- setdiff(rownames(reference), rownames(target))
  if (length(missing_genes) > 0) {
    pad <- Matrix::Matrix(0,
                          nrow = length(missing_genes),
                          ncol = ncol(target),
                          sparse = TRUE,
                          dimnames = list(missing_genes, colnames(target)))
    target <- rbind(target, pad)
  }
  target <- target[rownames(reference), , drop = FALSE]
  stopifnot(identical(rownames(reference), rownames(target)))
  target
}
