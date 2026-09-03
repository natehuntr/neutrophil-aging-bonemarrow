# ---------------------------------------------------------------------------
# Package loading.
#
# Only the packages that are used with `library()` semantics (Seurat's generic
# dispatch, the pipe-friendly tidyverse verbs) are attached; everything else is
# called with `::` so it is obvious where each function comes from.
# ---------------------------------------------------------------------------

REQUIRED_PACKAGES <- c(
  "Seurat", "Matrix", "dplyr", "tidyr", "tibble", "purrr", "ggplot2",
  "patchwork", "matrixStats", "pheatmap", "yaml"
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

#' Stop early if an optional, step-specific package is absent.
require_packages <- function(...) {
  pkgs <- c(...)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("this step needs: ", paste(missing, collapse = ", "),
         "\nRun: Rscript scripts/00_install_dependencies.R")
  invisible(TRUE)
}
