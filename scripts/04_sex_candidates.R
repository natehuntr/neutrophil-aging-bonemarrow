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

stage_sex <- table(paste(obj$sex, obj$stage), useNA = "no")
gates <- list(
  gate_depth_ratio(counts, groups, cfg),
  gate_stratum_sizes(stage_sex, cfg)
)
gate_results <- run_gates(gates, cfg, "sex contrast")
write_table(gate_results, cfg, "sex_control_gates.csv")

# ===========================================================================
# 2. Depth matching, and the diagnostics that justify it
# ===========================================================================
detection <- detection_rates(counts, groups)
write_confound_diagnostics(obj, cfg, detection = detection)

obj <- add_matched_assay(obj, cfg, group_col = "sex")
matched_counts <- Seurat::GetAssayData(obj, assay = "RNAmatched", layer = "counts")
log_step("depth ratio after matching: ",
         round(depth_ratio(matched_counts, groups), 3), "x")

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
score_cols <- grep("_UCell$", colnames(obj@meta.data), value = TRUE)

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

log_step("top candidates (ranked by what should be believed, not by effect size):")
print(utils::head(candidate_table[, c("candidate", "stage", "effect_size",
                                      "evidence_tier", "null_percentile",
                                      "survives_downsampling", "orthogonal_assay")], 15),
      row.names = FALSE)
log_step(attr(candidate_table, "caveat"))
log_step("step 4 complete")
