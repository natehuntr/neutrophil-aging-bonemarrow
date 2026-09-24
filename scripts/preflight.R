#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Preflight: can this pipeline run end to end, from nothing?
#
#   Rscript scripts/preflight.R
#
# A full run is roughly a day of compute, and most of the ways it fails are
# visible in seconds: a missing package, an input directory that is not where
# the config says, a config key a step reads but nobody wrote. This checks
# those before the first job is submitted rather than at hour nine.
#
# It does NOT check the analysis. It checks that every step could start.
# Exit status is 0 when everything a full run needs is present.
# ---------------------------------------------------------------------------

source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")

failures <- character()
warnings_seen <- character()

check <- function(label, ok, detail = "", fatal = TRUE) {
  status <- if (isTRUE(ok)) "ok  " else if (fatal) "FAIL" else "warn"
  cat(sprintf("  [%s] %-46s %s\n", status, label, detail))
  if (!isTRUE(ok)) {
    if (fatal) failures <<- c(failures, label)
    else warnings_seen <<- c(warnings_seen, label)
  }
  invisible(ok)
}

cat("\n--- configuration ---\n")
cfg <- tryCatch(load_config(), error = function(e) {
  check("config/config.yml parses", FALSE, conditionMessage(e)); NULL
})
if (is.null(cfg)) quit(status = 1)
check("config/config.yml parses", TRUE)

# Keys read by a step that has no default for them. A typo here surfaces as an
# obscure error deep in a script hours later.
required_keys <- list(
  "analysis$age_levels" = cfg$analysis$age_levels,
  "analysis$sex_levels" = cfg$analysis$sex_levels,
  "analysis$stage_levels" = cfg$analysis$stage_levels,
  "analysis$progenitor_levels" = cfg$analysis$progenitor_levels,
  "stage_assignment$method" = cfg$stage_assignment$method,
  "stage_assignment$rna_signatures" = cfg$stage_assignment$rna_signatures,
  "endpoint$reference" = cfg$endpoint$reference,
  "endpoint$endpoint" = cfg$endpoint$endpoint,
  "gates$min_cells_per_stratum" = cfg$gates$min_cells_per_stratum,
  "depth$match_on" = cfg$depth$match_on,
  "samples" = cfg$samples
)
for (key in names(required_keys))
  check(paste0("config key ", key), length(required_keys[[key]]) > 0)

signature_names <- names(cfg$stage_assignment$rna_signatures)
check("every stage_level has an RNA signature",
      all(unlist(cfg$analysis$stage_levels) %in% signature_names),
      paste("missing:", paste(setdiff(unlist(cfg$analysis$stage_levels),
                                      signature_names), collapse = ", ")))

check("endpoint ages are in analysis$age_levels",
      all(c(cfg$endpoint$reference, cfg$endpoint$endpoint) %in%
            unlist(cfg$analysis$age_levels)))

cat("\n--- modules ---\n")
modules <- tryCatch({ load_modules() },
                    error = function(e) { check("R/ modules source", FALSE,
                                                conditionMessage(e)); NULL })
if (!is.null(modules)) check("R/ modules source", TRUE,
                             paste(length(modules), "files"))

cat("\n--- packages ---\n")
by_step <- list(
  "1"  = c("Seurat", "SoupX", "scDblFinder", "scuttle", "SingleCellExperiment"),
  "2"  = c("SingleR", "celldex", "CytoTRACE2", "clustree"),
  "3"  = c("UCell", "DropletUtils"),
  "4"  = c("fgsea", "msigdbr", "UCell", "DropletUtils"),
  "5"  = c("monocle3", "tradeSeq", "SummarizedExperiment"),
  "6"  = c("matrixStats"),
  "7"  = c("glmGamPoi", "splines"),
  "8"  = c("fgsea", "msigdbr"),
  "9"  = c("monocle3"),
  "10" = c("fgsea", "msigdbr"),
  "11" = c("fgsea", "msigdbr", "DropletUtils")
)
for (step in names(by_step)) {
  pkgs <- by_step[[step]]
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  check(paste0("step ", step, " packages"), !length(missing),
        if (length(missing)) paste("missing:", paste(missing, collapse = ", ")) else "")
}

cat("\n--- inputs ---\n")
root <- cfg$paths$cellranger_root
check("cellranger_root exists", dir.exists(root), root)
if (dir.exists(root))
  for (key in names(cfg$samples)) {
    sample <- cfg$samples[[key]]
    for (which in c("raw_matrix_dir", "filtered_matrix_dir")) {
      path <- file.path(root, sample$sample_id, cfg$paths[[which]])
      check(paste0(key, " ", sub("_matrix_dir", "", which), " matrix"),
            dir.exists(path), path)
    }
  }

cat("\n--- step chain ---\n")
# Each step's script must exist and be registered in both runners, or
# submit_all will silently skip it.
steps <- c(1:11)
scripts <- list.files(project_path("scripts"), pattern = "^[0-9]+_.*[.]R$")
for (step in steps) {
  file <- grep(sprintf("^%02d_", step), scripts, value = TRUE)
  check(paste0("scripts/", sprintf("%02d", step), "_*.R present"),
        length(file) == 1, paste(file, collapse = ", "))
}
runner <- paste(readLines(project_path("scripts", "run_all.R")), collapse = "\n")
submitter <- paste(readLines(project_path("slurm", "submit_all.sh")), collapse = "\n")
for (step in steps) {
  check(paste0("step ", step, " registered in run_all.R"),
        grepl(sprintf('"%d" *=', step), runner), fatal = TRUE)
  check(paste0("step ", step, " has SLURM resources"),
        grepl(sprintf("\n *%d\\) echo \"--time", step), submitter), fatal = TRUE)
}

cat("\n", strrep("=", 66), "\n", sep = "")
if (length(failures)) {
  cat(length(failures), "check(s) failed:\n")
  cat(paste0("  - ", failures, collapse = "\n"), "\n")
  cat("A full run would not complete. Fix these first.\n")
  quit(status = 1)
}
if (length(warnings_seen))
  cat(length(warnings_seen), "warning(s); a run can proceed.\n")
cat("Ready: every step can start.\n")
cat(strrep("=", 66), "\n")
