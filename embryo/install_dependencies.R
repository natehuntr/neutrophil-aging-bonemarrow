#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Install every package embryonic_neutrophils.qmd needs. Run once per machine,
# from this folder:
#
#   Rscript install_dependencies.R
#
# Then record the exact versions so the analysis can be rebuilt later:
#
#   Rscript -e 'renv::snapshot(type = "all")'
#
# and commit the renv.lock it writes. On a new machine, renv::restore()
# reinstalls those exact versions instead of running this script.
# ---------------------------------------------------------------------------

options(Ncpus = max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "2"))))

# Pin CRAN to a dated snapshot so the same versions are installed on every
# machine. Override with CRAN_SNAPSHOT=... if your R is newer than this date.
snapshot <- Sys.getenv("CRAN_SNAPSHOT",
                       unset = "https://packagemanager.posit.co/cran/2025-06-01")
options(repos = c(CRAN = snapshot))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

install_missing <- function(pkgs, installer) {
  todo <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(todo)) installer(todo)
  invisible(todo)
}

cran <- c("Seurat", "SeuratObject", "Matrix", "data.table", "R.utils", "dplyr",
          "tidyr", "ggplot2", "patchwork", "pheatmap", "clustree", "dbscan",
          "leidenbase", "msigdbr", "shiny", "plotly", "igraph")
bioc <- c("SingleCellExperiment", "SummarizedExperiment", "scuttle", "scran",
          "scDblFinder", "celda", "SingleR", "celldex", "DESeq2", "fgsea",
          "dittoSeq", "ComplexHeatmap", "BiocNeighbors",
          # monocle3 dependencies from Bioconductor
          "BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma",
          "lme4", "S4Vectors", "batchelor", "HDF5Array", "terra", "ggrastr")

install_missing(cran, install.packages)
install_missing(bioc, function(p) BiocManager::install(p, update = FALSE, ask = FALSE))

# GitHub packages. These install the current default branch; renv.lock (see
# the top of this file) records the exact commit, which is what makes the
# install repeatable.
github <- c(
  monocle3       = "cole-trapnell-lab/monocle3",
  SeuratWrappers = "satijalab/seurat-wrappers",
  CellChat       = "jinworks/CellChat",
  CytoTRACE2     = "digitalcytometry/cytotrace2"
)
for (p in names(github)) {
  if (requireNamespace(p, quietly = TRUE)) next
  ref <- github[[p]]
  if (p == "CytoTRACE2") {
    remotes::install_github(ref, subdir = "cytotrace2_r", upgrade = "never")
  } else {
    remotes::install_github(ref, upgrade = "never")
  }
}

still_missing <- c(cran, bioc, names(github))
still_missing <- still_missing[!vapply(still_missing, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) {
  stop("Could not install: ", paste(still_missing, collapse = ", "))
}
cat("All packages installed. Now run: Rscript -e 'renv::snapshot(type = \"all\")'\n")
