# ---------------------------------------------------------------------------
# Reading Cell Ranger output.
#
# Only the filtered (cell-containing) matrix is read. The raw matrix was
# needed for ambient-RNA estimation; that step has been removed, so nothing
# in the pipeline looks at empty droplets any more.
# ---------------------------------------------------------------------------

#' Locate a sample's filtered matrix directory.
sample_matrix_path <- function(cfg, sample_id) {
  path <- file.path(cfg$paths$cellranger_root, sample_id, cfg$paths$filtered_matrix_dir)
  if (!dir.exists(path))
    stop("matrix directory not found:\n  ", path,
         "\nCheck paths.cellranger_root in config/config.yml.")
  path
}

#' Read one hashed CITE-seq sample.
#'
#' Returns the gene-expression and antibody matrices on the sample's cell
#' barcodes.
read_cite_sample <- function(cfg, sample_id) {
  path <- sample_matrix_path(cfg, sample_id)
  log_step("reading ", sample_id)

  filtered <- Seurat::Read10X(path)
  expected <- c("Gene Expression", "Antibody Capture")
  missing <- setdiff(expected, names(filtered))
  if (length(missing))
    stop(sample_id, ": Read10X returned no '", paste(missing, collapse = "', '"),
         "' matrix. Is this a multimodal (CITE-seq) run?")

  list(
    sample_id = sample_id,
    gex       = filtered[["Gene Expression"]],
    adt       = filtered[["Antibody Capture"]]
  )
}
