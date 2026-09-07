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

# Compile in parallel where the allocation allows it. detectCores() would
# report the whole node rather than this job's share.
NCPUS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
options(Ncpus = max(1L, NCPUS))

# Install into the first writable library, creating it if the site config
# named one that does not exist yet.
lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(lib)) {
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
}
cat("installing into:", .libPaths()[1], "\n")
cat("parallel jobs  :", getOption("Ncpus"), "\n\n")

cran <- c(
  "yaml", "here", "dplyr", "tidyr", "tibble", "purrr", "readr", "glue",
  "ggplot2", "patchwork", "pheatmap", "matrixStats", "Matrix", "scales",
  "Seurat", "SeuratObject", "R.utils", "ggVennDiagram", "clustree",
  "msigdbr", "remotes", "BiocManager"
)

bioc <- c(
  "scDblFinder", "glmGamPoi", "SingleR", "celldex", "dittoSeq",
  "SingleCellExperiment", "SummarizedExperiment", "fgsea", "tradeSeq",
  "clusterProfiler", "org.Mm.eg.db", "EnsDb.Mmusculus.v79", "batchelor"
)

# GitHub packages. Installed from the codeload archive endpoint rather than
# through install_github(), because that uses the GitHub *API*, whose
# unauthenticated limit is 60 requests per hour PER IP -- shared by everyone on
# a cluster's outbound address. It runs out, and the install fails with
# "HTTP error 403. API rate limit exceeded" having downloaded nothing.
# codeload is a plain file download with no such limit.
#
# Setting GITHUB_PAT (in slurm/env.local.sh) raises the API limit and is worth
# doing, but this path does not need it.
github <- list(
  list(name = "monocle3",       repo = "cole-trapnell-lab/monocle3"),
  list(name = "SeuratWrappers", repo = "satijalab/seurat-wrappers"),
  list(name = "CytoTRACE2",     repo = "digitalcytometry/cytotrace2",
       subdir = "cytotrace2_r")
)

#' Install one package from GitHub, trying each route that avoids the API.
#'
#' Three routes, because they fail for different reasons and clusters differ in
#' what they allow out:
#'   1. install_github  -- only when GITHUB_PAT is set, since unauthenticated it
#'      is the API path that hit the rate limit in the first place;
#'   2. codeload tarball -- a plain file download, no API involved;
#'   3. install_git      -- git clone over HTTPS, for sites where the archive
#'      endpoint is blocked but git is not.
install_from_github <- function(name, repo, subdir = NULL,
                                branches = c("master", "main")) {
  if (requireNamespace(name, quietly = TRUE)) {
    message("  ", name, ": already installed")
    return(TRUE)
  }

  installed <- function() requireNamespace(name, quietly = TRUE)
  attempt <- function(label, expr) {
    message("  ", name, ": trying ", label)
    ok <- tryCatch({ force(expr); installed() },
                   error = function(e) { message("    ", conditionMessage(e)); FALSE },
                   warning = function(w) { message("    ", conditionMessage(w)); installed() })
    isTRUE(ok)
  }

  # 1. Authenticated API, if a token is available.
  if (nzchar(Sys.getenv("GITHUB_PAT"))) {
    if (attempt("install_github (GITHUB_PAT is set)",
                remotes::install_github(paste0(repo, if (!is.null(subdir)) paste0("/", subdir) else ""),
                                        upgrade = "never")))
      return(TRUE)
  }

  # 2. Source archive: a plain download, no API.
  for (branch in branches) {
    url <- sprintf("https://codeload.github.com/%s/tar.gz/refs/heads/%s", repo, branch)
    dest <- file.path(tempdir(), paste0(name, "-", branch, ".tar.gz"))

    got <- tryCatch({
      utils::download.file(url, dest, quiet = TRUE, mode = "wb")
      file.exists(dest) && file.size(dest) > 1000
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!isTRUE(got)) next

    exdir <- file.path(tempdir(), paste0(name, "-src"))
    unlink(exdir, recursive = TRUE)
    dir.create(exdir, recursive = TRUE)
    utils::untar(dest, exdir = exdir)
    root <- list.dirs(exdir, recursive = FALSE)[1]
    path <- if (is.null(subdir)) root else file.path(root, subdir)
    if (!dir.exists(path)) next

    if (attempt(paste0("source archive ", repo, "@", branch),
                remotes::install_local(path, upgrade = "never", dependencies = TRUE)))
      return(TRUE)
  }

  # 3. git clone, for sites where codeload is blocked but git is not.
  if (attempt("git clone",
              remotes::install_git(paste0("https://github.com/", repo, ".git"),
                                   subdir = subdir, upgrade = "never")))
    return(TRUE)

  message("  ", name, ": FAILED -- all routes exhausted")
  FALSE
}

for (spec in github)
  do.call(install_from_github, spec)

# ---------------------------------------------------------------------------
# Report. install.packages() warns rather than errors on a failure, so without
# this a run that installed nothing still looks like it worked.
# ---------------------------------------------------------------------------
required <- c(
  # what R/packages.R attaches
  "Seurat", "Matrix", "dplyr", "tidyr", "tibble", "purrr", "ggplot2",
  "patchwork", "matrixStats", "pheatmap", "yaml", "scales",
  # what individual steps need
  "scDblFinder", "SingleCellExperiment", "SummarizedExperiment", "SingleR",
  "celldex", "CytoTRACE2", "clustree", "glmGamPoi", "monocle3",
  "SeuratWrappers", "tradeSeq", "fgsea", "msigdbr", "clusterProfiler",
  "org.Mm.eg.db", "ggVennDiagram"
)
status <- vapply(required, requireNamespace, logical(1), quietly = TRUE)

cat("\n", strrep("=", 60), "\n", sep = "")
cat(sum(status), "of", length(status), "required packages available\n")
cat(strrep("=", 60), "\n")
if (any(!status)) {
  cat("MISSING:\n")
  cat(paste0("  ", names(status)[!status], collapse = "\n"), "\n")
  cat("\nSearch the install log for each name to see why.\n")
} else {
  cat("Everything the pipeline needs is installed.\n")
}

# BPCells and speedglm are optional dependencies of Seurat and monocle3
# respectively; both want a Matrix newer than this R can take. Neither is used
# by this pipeline, so their failure is not a problem.
cat("\nsessionInfo():\n")
print(utils::sessionInfo())

quit(status = as.integer(any(!status)))
