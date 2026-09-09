# ---------------------------------------------------------------------------
# Reading Cell Ranger output.
#
# The raw (all-droplet) matrix is read again when ambient correction is on:
# SoupX estimates the ambient profile from the empty droplets, so the filtered
# matrix alone is not enough.
# ---------------------------------------------------------------------------

#' Locate a sample's matrix directories.
sample_matrix_paths <- function(cfg, sample_id, need_raw = FALSE) {
  root <- cfg$paths$cellranger_root
  paths <- list(filtered = file.path(root, sample_id, cfg$paths$filtered_matrix_dir))
  if (need_raw) paths$raw <- file.path(root, sample_id, cfg$paths$raw_matrix_dir)

  absent <- paths[!vapply(paths, dir.exists, logical(1))]
  if (length(absent))
    stop("matrix directory not found:\n  ", paste(unlist(absent), collapse = "\n  "),
         "\nCheck paths.cellranger_root in config/config.yml.",
         if (need_raw) paste0("\nThe raw matrix is needed because ambient.method is '",
                              cfg$ambient$method, "'; set it to 'none' to skip it.") else "")
  paths
}

#' Read one hashed CITE-seq sample.
read_cite_sample <- function(cfg, sample_id, need_raw = FALSE) {
  paths <- sample_matrix_paths(cfg, sample_id, need_raw)
  log_step("reading ", sample_id)

  filtered <- Seurat::Read10X(paths$filtered)
  expected <- c("Gene Expression", "Antibody Capture")
  absent <- setdiff(expected, names(filtered))
  if (length(absent))
    stop(sample_id, ": Read10X returned no '", paste(absent, collapse = "', '"),
         "' matrix. Is this a multimodal (CITE-seq) run?")

  out <- list(sample_id = sample_id,
              gex = filtered[["Gene Expression"]],
              adt = filtered[["Antibody Capture"]])

  if (need_raw) {
    log_step("  reading the raw matrix for ambient estimation")
    raw <- Seurat::Read10X(paths$raw)
    raw_gex <- raw[["Gene Expression"]]
    # SoupChannel needs both matrices on identical row sets. Padding rather
    # than intersecting keeps every gene Cell Ranger saw in the droplets.
    out$raw_gex <- raw_gex
    out$gex <- align_gene_rows(raw_gex, out$gex)
  }
  out
}

#' Zero-pad `target` so its rows match `reference`, then reorder to match.
align_gene_rows <- function(reference, target) {
  missing_genes <- setdiff(rownames(reference), rownames(target))
  if (length(missing_genes) > 0) {
    pad <- Matrix::Matrix(0, nrow = length(missing_genes), ncol = ncol(target),
                          sparse = TRUE,
                          dimnames = list(missing_genes, colnames(target)))
    target <- rbind(target, pad)
  }
  target <- target[rownames(reference), , drop = FALSE]
  stopifnot(identical(rownames(reference), rownames(target)))
  target
}
