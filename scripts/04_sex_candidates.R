#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 4 - candidates for orthogonal validation from the sex contrast.
#
#   Rscript scripts/04_sex_candidates.R
#
# This step does NOT test for sex differences, and is built so it cannot be
# read as though it does. Sex is confounded with library, run, probe barcode
# and staining batch; one hashtag per age leaves no biological replication; the
# mice are inbred, so no demultiplexing recovers individuals. What comes out is
# a ranked, robustness-annotated candidate list, and every column exists to
# help decide whether a candidate is worth an orthogonal experiment.
#
# The order matters and is not negotiable:
#   1. gates            - refuse to proceed on a comparison that cannot work
#   2. depth matching   - the primary path, not a sensitivity check
#   3. module scores    - set-level, within ADT-defined stage
#   4. depth controls   - does the effect survive matching, and is it created
#                         by matching alone?
#   5. within-library null - the only quantitative floor available
#   6. candidate table  - effect sizes with CIs, no p-values
#
# Reads:  results/objects/gmp_neutrophils.rds
# Writes: results/tables/sex_candidates.csv        <- the deliverable
#         results/tables/sex_control_gates.csv
#         results/tables/sex_module_effects.csv
#         results/tables/within_library_null.csv
#         results/figures/confound_*.pdf
# ---------------------------------------------------------------------------

source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()
require_packages("glmGamPoi", "UCell", "DropletUtils", "SingleCellExperiment",
                 "SummarizedExperiment", "fgsea", "msigdbr")

obj <- read_object(cfg, "gmp_neutrophils.rds")
require_metadata(obj, c("sex", "age", "stage"), context = "step 4")
sexes <- cfg$analysis$sex_levels

# ===========================================================================
# 1. Gates
# ===========================================================================
counts <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
groups <- as.character(obj$sex)

panel <- check_sex_genes_in_panel(obj, cfg)
write_table(panel, cfg, "sex_gene_panel_check.csv")
log_step("sex-chromosome genes in the probe set: ",
         paste(panel$gene[panel$in_panel], collapse = ", "), " (of ",
         paste(panel$gene, collapse = ", "), ")")

# Gates that do not depend on depth matching are checked first, so a stratum
# problem is reported before minutes are spent thinning counts.
stage_sex <- table(paste(obj$sex, obj$stage), useNA = "no")
# require_all = FALSE: undersized stages are skipped by the per-stage loop
# below, exactly as steps 6 and 8 skip theirs. Only a complete absence of
# usable strata is worth stopping for.
pre_gates <- run_gates(list(gate_stratum_sizes(stage_sex, cfg, require_all = FALSE)),
                       cfg, "sex contrast, before matching")

# ===========================================================================
# 2. Depth matching, and the diagnostics that justify it
# ===========================================================================
log_step("depth ratio before matching: ",
         round(depth_ratio(counts, groups), 3), "x")
detection <- detection_rates(counts, groups)
write_confound_diagnostics(obj, cfg, detection = detection)

# Built in step 3 on depth.match_on (age_sex), which is finer than the sex
# grouping this step contrasts, so it is never less conservative here.
matched_assay <- matched_assay_for(obj, cfg, step = 4)
if (!identical(matched_assay, "RNAmatched"))
  stop("step 4 cannot run on unmatched counts: the depth-matched assay IS the ",
       "analysis. Put 4 back in depth.matched_steps.")
matched_counts <- Seurat::GetAssayData(obj, assay = matched_assay, layer = "counts")

# The depth gate belongs HERE, on the matched counts -- those are what every
# effect size below is computed from. Checking it before matching would fail
# on the raw libraries every time, which is the known starting condition
# rather than a finding.
gate_results <- rbind(
  pre_gates,
  run_gates(list(gate_depth_ratio(matched_counts, groups, cfg,
                                  label = "depth-matched counts")),
            cfg, "sex contrast, after matching"))
write_table(gate_results, cfg, "sex_control_gates.csv")

# ===========================================================================
# 3. Module scores, within ADT-defined stage
# ===========================================================================
# Set-level first: coherence across a module survives depth noise that no
# single gene survives.
gs <- gobp_pathways(rownames(obj), cfg)
gs <- gs[lengths(gs) >= cfg$gsea$min_set_size & lengths(gs) <= cfg$gsea$max_set_size]
# Only modules with a proposed orthogonal assay are worth scoring: a candidate
# that cannot be tested is not a candidate.
testable <- names(gs)[suggest_orthogonal_assay(names(gs)) != "none proposed -- not a candidate until one is"]
gs <- gs[testable]
log_step(length(gs), " gene sets have a proposed orthogonal assay")

Seurat::DefaultAssay(obj) <- "RNAmatched"
obj <- add_module_scores_ucell(obj, gs, assay = "RNAmatched")
# Name the columns from the sets that were just scored rather than grepping the
# suffix: other steps leave their own UCell scores on the object, and a grep
# would sweep those in as candidate modules.
score_cols <- intersect(paste0(names(gs), "_UCell"), colnames(obj@meta.data))
if (!length(score_cols)) stop("no module scores were written; nothing to rank")

# ===========================================================================
# 4-5. Per stage: effects, depth controls, and the null floor
# ===========================================================================
candidates <- list()
null_summaries <- list()

for (stage in cfg$analysis$stage_levels) {
  cells <- colnames(obj)[!is.na(obj$stage) & as.character(obj$stage) == stage]
  if (!length(cells)) next
  sub <- subset(obj, cells = cells)

  n_by_sex <- table(as.character(sub$sex))
  if (!stratum_is_usable(n_by_sex, stage, cfg)) next
  log_step("=== ", stage, " (", paste(names(n_by_sex), n_by_sex, sep = "=",
                                      collapse = ", "), ") ===")

  effects <- do.call(rbind, lapply(score_cols, function(col)
    module_effect(sub, col, "sex", sexes)))
  effects$stage <- stage

  # Does the same module separate two halves of ONE library? Anything it finds
  # there is the floor for what it finds between the sexes.
  null_effects <- within_library_null(
    sub, cfg,
    analysis = function(split_obj, assignment) {
      split_obj$.split <- assignment
      max(abs(vapply(score_cols, function(col)
        module_effect(split_obj, col, ".split", c("split_a", "split_b"))$effect_size,
        numeric(1))), na.rm = TRUE)
    },
    stage_col = "stage", target_n = min(n_by_sex))

  effects$null_floor <- null_effects$floor
  effects$null_percentile <- vapply(abs(effects$effect_size),
                                    null_percentile, numeric(1),
                                    null_result = null_effects)
  null_summaries[[stage]] <- compare_to_null(
    max(abs(effects$effect_size), na.rm = TRUE), null_effects, stage)

  # Did it survive matching? Compare the matched effect with the unmatched one.
  Seurat::DefaultAssay(sub) <- "RNA"
  sub <- add_module_scores_ucell(sub, gs, assay = "RNA", suffix = "_raw")
  raw_cols <- sub("_UCell$", "_raw", score_cols)
  unmatched <- do.call(rbind, lapply(intersect(raw_cols, colnames(sub@meta.data)),
    function(col) module_effect(sub, col, "sex", sexes)))

  effects$survives_downsampling <- mapply(
    survives_matching,
    observed = unmatched$effect_size[match(sub("_UCell$", "_raw", effects$module),
                                           unmatched$module)],
    matched = effects$effect_size)

  candidates[[stage]] <- effects
}

if (!length(candidates))
  stop("no stratum passed the minimum-cell gate, so there are no candidates.\n",
       "This is a statement about the design, not a pipeline failure: see the ",
       "stratum sizes in results/tables/sex_control_gates.csv.")

effects <- do.call(rbind, candidates)
write_table(strip_inferential_columns(effects), cfg, "sex_module_effects.csv")
write_table(do.call(rbind, null_summaries), cfg, "within_library_null.csv")

# ===========================================================================
# 6. The deliverable
# ===========================================================================
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

# Only rows that clear the within-library null belong in a "top candidates"
# list. A row below the null floor is one the data argues against -- a random
# split of one library separated the sexes better -- and printing it under
# that heading is the one thing this table exists not to do. The full ranking,
# floor rows included, is in the csv.
believable <- candidate_table[!candidate_table$below_null_floor, ]
n_floor <- sum(candidate_table$below_null_floor, na.rm = TRUE)

if (nrow(believable)) {
  log_step("top candidates (ranked by what should be believed, not by effect size):")
  print(utils::head(believable[, c("candidate", "stage", "effect_size",
                                   "evidence_tier", "null_percentile",
                                   "exceeds_null", "survives_downsampling",
                                   "orthogonal_assay")], 15), row.names = FALSE)
} else {
  log_step("NO candidate clears the within-library null floor.")
  log_step("  Every module scored here separates the sexes less well than a ",
           "random split of the female library separates its own cells. ",
           "That is a result: this contrast has nothing to hand an orthogonal ",
           "assay. See within_library_null.csv for the floors.")
}
log_step(n_floor, " of ", nrow(candidate_table),
         " rows sit below the null floor and are ranked last in the csv.")
log_step(attr(candidate_table, "caveat"))
log_step("step 4 complete")
