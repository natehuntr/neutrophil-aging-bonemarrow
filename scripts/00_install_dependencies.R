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
cat("parallel jobs  :", getOption("Ncpus"), "\n")

# ---------------------------------------------------------------------------
# Repositories.
#
# Current CRAN only carries the newest version of each package, and those
# increasingly require a newer R than a cluster module provides. R then reports
# "package 'Matrix' is not available for this version of R" and every package
# depending on it fails too -- one incompatible recommended package takes out
# twenty others.
#
# The snapshot is listed FIRST, and this matters more than it looks.
# BiocManager pins Bioconductor to the release for this R -- 3.18 for R 4.3.2,
# from late 2023. Taking CRAN dependencies from today instead pairs 2023
# Bioconductor sources with 2026 CRAN headers, and they no longer agree about
# the C++ standard: Bioconductor packages of that era pin CXX_STD = CXX11,
# while current BH and RcppArmadillo need C++14 or later. The result is
# hundreds of lines of "'is_final' has not been declared in 'std'" and
# "C++14 compiler required" from fgsea, glmGamPoi and treeio.
#
# A snapshot from the same era as the Bioconductor release keeps both halves
# consistent. Current CRAN stays as a fallback for anything the snapshot lacks.
SNAPSHOT <- Sys.getenv("CRAN_SNAPSHOT",
                       unset = "https://packagemanager.posit.co/cran/2024-04-15")
options(repos = c(SNAPSHOT = SNAPSHOT, CRAN = "https://cloud.r-project.org"))
cat("repositories   :\n")
cat(paste0("  ", names(getOption("repos")), ": ", getOption("repos")), sep = "\n")
cat("\n")

cran <- c(
  "yaml", "here", "dplyr", "tidyr", "tibble", "purrr", "readr", "glue",
  "ggplot2", "patchwork", "pheatmap", "matrixStats", "Matrix", "scales",
  "Seurat", "SeuratObject", "R.utils", "ggVennDiagram", "clustree",
  "msigdbr", "remotes", "BiocManager",
  # monocle3's dependencies, installed here rather than left to it, so the
  # API-free route below has everything it needs already present. sf and
  # spdep need GDAL/GEOS/PROJ; units needs UDUNITS. See slurm/env.local.sh.
  "sf", "units", "spdep", "terra", "leidenbase", "RhpcBLASctl", "lme4", "pscl"
)

# Packages CRAN has archived, which no repository serves any more. monocle3
# still declares speedglm as a dependency.
archived <- list(
  list(name = "speedglm",
       url = "https://cran.r-project.org/src/contrib/Archive/speedglm/speedglm_0.3-5.tar.gz")
)

bioc <- c(
  "scDblFinder", "glmGamPoi", "SingleR", "celldex", "dittoSeq",
  "SingleCellExperiment", "SummarizedExperiment", "fgsea", "tradeSeq",
  "clusterProfiler", "org.Mm.eg.db", "EnsDb.Mmusculus.v79", "batchelor"
)

#' Install only what is absent, and report what was attempted.
install_missing <- function(pkgs, installer, label) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  cat("\n", label, ": ", length(missing), " of ", length(pkgs),
      " to install\n", sep = "")
  if (!length(missing)) return(invisible(character()))
  cat("  ", paste(missing, collapse = ", "), "\n", sep = "")
  installer(missing)
  invisible(missing)
}

install_missing(cran, function(p) install.packages(p), "CRAN")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager picks the Bioconductor release matching this R. Its repositories
# are added to the CRAN ones already set, rather than replacing them, so the
# snapshot fallback still applies to Bioconductor's own CRAN dependencies.
if (requireNamespace("BiocManager", quietly = TRUE)) {
  cat("\nBioconductor version: ",
      as.character(BiocManager::version()), "\n", sep = "")
  install_missing(bioc,
                  function(p) BiocManager::install(p, ask = FALSE, update = FALSE),
                  "Bioconductor")
} else {
  warning("BiocManager could not be installed; every Bioconductor package ",
          "will be missing")
}

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

  # A token that is wrong is worse than no token at all: remotes sends it on
  # every request, including the ones that would otherwise be anonymous, and
  # GitHub answers 401 Bad credentials rather than serving the file. Check it
  # once and drop it if it is not accepted.
  if (nzchar(Sys.getenv("GITHUB_PAT"))) {
    valid <- tryCatch({
      con <- url("https://api.github.com/rate_limit", open = "r",
                 headers = c(Authorization = paste("Bearer", Sys.getenv("GITHUB_PAT"))))
      on.exit(close(con), add = TRUE)
      length(readLines(con, n = 1, warn = FALSE)) > 0
    }, error = function(e) FALSE)
    if (!isTRUE(valid)) {
      message("  GITHUB_PAT is set but GitHub rejects it -- ignoring it. ",
              "Check for a truncated or expired token in slurm/env.local.sh.")
      Sys.unsetenv("GITHUB_PAT")
    }
  }

  # 1. Authenticated API, if a usable token is available.
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

    # remotes reads the package's Remotes: field and resolves those entries
    # through the GitHub API even when the source is already on disk, so the
    # rate limit bites a local install too. R CMD INSTALL ignores Remotes
    # entirely; the dependencies it needs are in the CRAN list above.
    if (attempt(paste0("R CMD INSTALL of ", repo, "@", branch, " (ignores Remotes:)"),
                utils::install.packages(path, repos = NULL, type = "source")))
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

for (spec in archived) {
  if (requireNamespace(spec$name, quietly = TRUE)) next
  message("  ", spec$name, ": installing from the CRAN archive")
  tryCatch(install.packages(spec$url, repos = NULL, type = "source"),
           error = function(e) message("    ", conditionMessage(e)))
}

for (spec in github)
  do.call(install_from_github, spec)

# ---------------------------------------------------------------------------
# Report. install.packages() warns rather than errors on a failure, so without
# this a run that installed nothing still looks like it worked.
# ---------------------------------------------------------------------------
# Which steps need what. A package missing here blocks only the steps listed
# against it, and saying so is more useful than a flat list -- SeuratWrappers,
# for instance, is used only when building a trajectory from a Seurat object,
# so a finished step 5 does not need it again.
needed_by <- list(
  "Seurat"               = "all steps",
  "Matrix"               = "all steps",
  "dplyr"                = "all steps",
  "tidyr"                = "all steps",
  "tibble"               = "all steps",
  "purrr"                = "all steps",
  "ggplot2"              = "all steps",
  "patchwork"            = "all steps",
  "matrixStats"          = "all steps",
  "pheatmap"             = "6",
  "yaml"                 = "all steps",
  "scales"               = "9",
  "scDblFinder"          = "1",
  "SingleCellExperiment" = "1, 2, 5, 7, 9",
  "SummarizedExperiment" = "1, 2, 5, 7, 9",
  "SingleR"              = "2",
  "celldex"              = "2",
  "CytoTRACE2"           = "2",
  "clustree"             = "2",
  "clusterProfiler"      = "4",
  "org.Mm.eg.db"         = "4",
  "monocle3"             = "5, 7, 9",
  "SeuratWrappers"       = "5",
  "tradeSeq"             = "5",
  "glmGamPoi"            = "2, 3, 7",
  "fgsea"                = "8",
  "msigdbr"              = "8",
  "ggVennDiagram"        = "5"
)
required <- names(needed_by)
status <- vapply(required, requireNamespace, logical(1), quietly = TRUE)

cat("\n", strrep("=", 68), "\n", sep = "")
cat(sum(status), "of", length(status), "required packages available\n")
cat(strrep("=", 68), "\n")
if (any(!status)) {
  missing_tbl <- data.frame(package = required[!status],
                            blocks_steps = unlist(needed_by[!status]),
                            row.names = NULL)
  print(missing_tbl, row.names = FALSE)
  cat("\nSearch the install log for each name to see why.\n")
  blocked <- sort(unique(unlist(strsplit(missing_tbl$blocks_steps, ", "))))
  cat("Steps affected:", paste(blocked, collapse = ", "), "\n")
} else {
  cat("Everything the pipeline needs is installed.\n")
}

# BPCells and speedglm are optional dependencies of Seurat and monocle3
# respectively; both want a Matrix newer than this R can take. Neither is used
# by this pipeline, so their failure is not a problem.
cat("\nsessionInfo():\n")
print(utils::sessionInfo())

quit(status = as.integer(any(!status)))
