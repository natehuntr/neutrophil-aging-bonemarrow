#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Run the whole pipeline in order.
#
#   Rscript scripts/run_all.R           # every step
#   Rscript scripts/run_all.R 3 4 5     # only those steps
#
# Each step is run in a fresh R process so a failure in one does not leave
# state behind, and memory is released between steps.
# ---------------------------------------------------------------------------

steps <- c(
  "1" = "01_preprocess.R",
  "2" = "02_annotate.R",
  "3" = "03_merge.R",
  "4" = "04_sex_differences.R",
  "5" = "05_trajectory.R",
  "6" = "06_age_trends.R",
  "7" = "07_glm_models.R",
  "8" = "08_gsea.R"
)

args <- commandArgs(trailingOnly = TRUE)
selected <- if (length(args)) steps[args] else steps
if (anyNA(selected))
  stop("unknown step(s): ", paste(args[is.na(steps[args])], collapse = ", "),
       "\nAvailable: ", paste(names(steps), collapse = ", "))

script_dir <- if (dir.exists("scripts")) "scripts" else "."
rscript <- file.path(R.home("bin"), "Rscript")

for (step in selected) {
  message("\n==================== ", step, " ====================")
  started <- Sys.time()
  status <- system2(rscript, file.path(script_dir, step))
  if (status != 0) stop(step, " failed with status ", status)
  message(step, " finished in ",
          round(difftime(Sys.time(), started, units = "mins"), 1), " min")
}

message("\nPipeline complete.")
