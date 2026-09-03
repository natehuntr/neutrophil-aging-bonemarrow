# ---------------------------------------------------------------------------
# Package loading.
#
# Only the packages that are used with `library()` semantics (Seurat's generic
# dispatch, the pipe-friendly tidyverse verbs) are attached; everything else is
# called with `::` so it is obvious where each function comes from.
# ---------------------------------------------------------------------------

REQUIRED_PACKAGES <- c(
  "Seurat", "Matrix", "dplyr", "tidyr", "tibble", "purrr", "ggplot2",
  "patchwork", "matrixStats", "pheatmap", "yaml", "scales"
)

#' Attach the core packages, failing with a single actionable message if any
#' are missing rather than one error per library() call.
load_packages <- function(extra = character()) {
  needed <- c(REQUIRED_PACKAGES, extra)
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("missing packages: ", paste(missing, collapse = ", "),
         "\nRun: Rscript scripts/00_install_dependencies.R")

  suppressPackageStartupMessages(
    for (p in needed) library(p, character.only = TRUE))
  invisible(needed)
}

#' Attach the step-specific packages, failing early if any are absent.
#'
#' These are ATTACHED, not merely checked, and the difference matters. S4
#' method dispatch in the Bioconductor stack resolves generics through the
#' search path, so calling everything as `Pkg::fun()` is not enough: a method
#' that runs inside another package can still fail to find a generic it needs.
#'
#' The concrete case that motivated this: scDblFinder internally runs
#' `counts(sce) <- ...`, whose SingleCellExperiment method body calls
#' `assay<-` from SummarizedExperiment. With SummarizedExperiment loaded but
#' not attached, that lookup fails with
#'   Error in assay(object, "counts") <- value :
#'     could not find function "assay<-"
#' which points at a package the calling code never mentions. Attaching
#' SummarizedExperiment fixes it. The same applies to monocle3 and
#' SingleCellExperiment generics elsewhere in the pipeline.
require_packages <- function(...) {
  pkgs <- c(...)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("this step needs: ", paste(missing, collapse = ", "),
         "\nRun: Rscript scripts/00_install_dependencies.R")

  suppressPackageStartupMessages(
    for (p in pkgs) library(p, character.only = TRUE))
  invisible(pkgs)
}
