#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Install everything the pipeline needs. Run once per machine:
#
#   Rscript scripts/00_install_dependencies.R
#
# Nothing here is installed automatically by the analysis scripts: they fail
# with a clear message instead, so a long run never stops halfway to compile a
# package.
# ---------------------------------------------------------------------------

cran <- c(
  "yaml", "here", "dplyr", "tidyr", "tibble", "purrr", "readr", "glue",
  "ggplot2", "patchwork", "pheatmap", "matrixStats", "Matrix",
  "Seurat", "SeuratObject", "R.utils", "ggVennDiagram", "clustree",
  "msigdbr", "remotes", "BiocManager"
)

bioc <- c(
  "SoupX", "scDblFinder", "glmGamPoi", "SingleR", "celldex", "dittoSeq",
  "SingleCellExperiment", "SummarizedExperiment", "fgsea", "tradeSeq",
  "clusterProfiler", "org.Mm.eg.db", "EnsDb.Mmusculus.v79", "batchelor"
)

github <- c(
  monocle3       = "cole-trapnell-lab/monocle3",
  SeuratWrappers = "satijalab/seurat-wrappers",
  CytoTRACE2     = "digitalcytometry/cytotrace2/cytotrace2_r"
)

install_missing <- function(pkgs, installer) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) installer(missing)
  invisible(missing)
}

install_missing(cran, function(p) install.packages(p, repos = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
install_missing(bioc, function(p) BiocManager::install(p, ask = FALSE, update = FALSE))

for (pkg in names(github)) {
  if (!requireNamespace(pkg, quietly = TRUE))
    remotes::install_github(github[[pkg]])
}

cat("\nInstalled versions:\n")
print(utils::sessionInfo())
