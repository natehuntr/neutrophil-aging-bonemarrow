#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Rebuild sex_candidates.csv from step 4's saved effects, without re-running it.
#
#   Rscript scripts/rescore_candidates.R
#
# Step 4 costs about 3.5 hours, and nearly all of it is the within-library
# null: 200 matched splits per stage, each rescoring every module. The ranking
# that consumes those numbers is arithmetic on a table.
#
# So when build_candidate_table() changes -- as it did when rows below the
# null floor stopped being ranked by effect size -- the expensive part does
# not need repeating. sex_module_effects.csv already holds every column the
# ranking reads, including the null percentiles the splits produced.
#
# Re-run step 4 itself if the CELLS change: new QC, new staging, a different
# matched assay. This script only re-applies the ranking.
#
# Reads:  results/tables/sex_module_effects.csv
# Writes: results/tables/sex_candidates.csv
# ---------------------------------------------------------------------------

source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()

path <- table_path(cfg, "sex_module_effects.csv")
if (!file.exists(path))
  stop("missing ", path, "\nRun step 4 once before rescoring its output.")

effects <- utils::read.csv(path, stringsAsFactors = FALSE)
log_step("read ", nrow(effects), " module effects from ", basename(path))

needed <- c("module", "effect_size", "ci_lower", "ci_upper", "direction",
            "null_percentile", "survives_downsampling", "stage", "n_min")
absent <- setdiff(needed, names(effects))
if (length(absent))
  stop("sex_module_effects.csv is missing: ", paste(absent, collapse = ", "),
       "\nIt predates this script; re-run step 4 to regenerate it.")

table_in <- data.frame(
  candidate = sub("_UCell$", "", effects$module),
  effect_size = effects$effect_size,
  ci_lower = effects$ci_lower,
  ci_upper = effects$ci_upper,
  direction = effects$direction,
  null_percentile = effects$null_percentile,
  survives_downsampling = effects$survives_downsampling,
  stage = effects$stage,
  n_cells_min = effects$n_min,
  stringsAsFactors = FALSE
)

candidate_table <- build_candidate_table(table_in, cfg)
write_table(candidate_table, cfg, "sex_candidates.csv")

believable <- candidate_table[!candidate_table$below_null_floor, ]
n_floor <- sum(candidate_table$below_null_floor, na.rm = TRUE)

if (nrow(believable)) {
  log_step("candidates clearing the within-library null floor:")
  print(utils::head(believable[, c("candidate", "stage", "effect_size",
                                   "evidence_tier", "null_percentile",
                                   "exceeds_null", "survives_downsampling",
                                   "orthogonal_assay")], 15), row.names = FALSE)
} else {
  log_step("NO candidate clears the within-library null floor.")
  log_step("  Every module separates the sexes less well than a random split ",
           "of the female library separates its own cells.")
}
log_step(n_floor, " of ", nrow(candidate_table),
         " rows sit below the null floor and are ranked last.")
log_step(attr(candidate_table, "caveat"))
