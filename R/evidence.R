# ---------------------------------------------------------------------------
# Evidence framing.
#
# This pipeline does not test sex differences. Sex is completely confounded
# with library, run, probe barcode and staining batch, there is one hashtag per
# age, and the mice are inbred, so no demultiplexing route recovers individual
# animals. What it produces is a ranked, robustness-annotated candidate list
# for orthogonal validation.
#
# The functions here exist to make that framing structural rather than a
# caveat: a p-value column that survives into a table gets read as evidence no
# matter what the footnote says, so the tables are built without one.
# ---------------------------------------------------------------------------

#' One-line statement of what a cross-sex number is and is not.
#'
#' Attached to every cross-sex table that gets written, so the constraint
#' travels with the file rather than living in a methods paragraph.
sex_contrast_caveat <- function(cfg) {
  sprintf(paste0("Sex is confounded with %s. %d mice pooled per hashtag (%d total): ",
                 "pooling improves precision, not accuracy -- the confounded quantity is ",
                 "measured very well and remains confounded -- and between-animal variance ",
                 "is unmeasurable, so precision cannot be quantified either. No inferential ",
                 "statistic here describes mice."),
          paste(cfg$design$sex_confounded_with, collapse = ", "),
          cfg$design$mice_per_hashtag, cfg$design$total_mice)
}

#' Annotate results by whether they run with or against the known bias.
#'
#' The confound has a direction: the deeper library captures more per cell, so
#' it is the one that will look "high" for anything depth-sensitive. A result
#' pointing the other way had to overcome that, and is worth more than a
#' nominally identical result pointing along it. This single column changes how
#' every downstream table reads.
#'
#' @param direction character vector naming the sex each row is higher in.
annotate_bias <- function(results, direction, cfg) {
  deeper <- cfg$design$deeper_library
  direction <- as.character(direction)

  results$direction <- direction
  results$against_bias <- !is.na(direction) & direction != deeper
  results$evidence_tier <- ifelse(is.na(direction), NA_character_,
                           ifelse(results$against_bias,
                                  "resists confound", "aligned with confound"))
  attr(results, "bias_direction") <- deeper
  results
}

#' Refuse to analyse a stratum too small to mean anything.
#'
#' Returns TRUE when the stratum is usable. Logs the reason when it is not, so
#' a missing row in the output is traceable to a decision rather than a bug.
stratum_is_usable <- function(counts, label, cfg, min_cells = NULL) {
  min_cells <- min_cells %||% cfg$gates$min_cells_per_stratum
  counts <- counts[!is.na(counts)]

  if (!length(counts)) {
    log_step("  skipping ", label, ": no cells")
    return(FALSE)
  }
  if (min(counts) < min_cells) {
    log_step(sprintf("  skipping %s: smallest group has %d cells (gate is %d)",
                     label, min(counts), min_cells))
    return(FALSE)
  }
  TRUE
}

#' Drop the columns that would be read as inference.
#'
#' Per-cell p-values answer a question about cells within one library, which is
#' not the question being asked of them. They are removed rather than kept with
#' a caveat.
strip_inferential_columns <- function(df, keep = character()) {
  inferential <- setdiff(
    grep("^(p_val|pval|pvalue|p\\.value|padj|adj_pval|p_val_adj|q_value|qvalue|FDR|fdr)$",
         names(df), value = TRUE, ignore.case = TRUE),
    keep)
  if (length(inferential)) {
    df <- df[, setdiff(names(df), inferential), drop = FALSE]
    # Set after subsetting: `[.data.frame` drops non-standard attributes.
    attr(df, "removed_columns") <- inferential
  }
  df
}

#' Bootstrap confidence interval for a per-cell statistic.
#'
#' Resamples cells, so the interval describes sampling of cells within the
#' library and nothing wider. Reported because an effect size without one
#' invites a reader to treat the point estimate as exact.
bootstrap_ci <- function(values, statistic = mean, n_boot = 1000, conf = 0.95) {
  values <- values[is.finite(values)]
  if (length(values) < 3) return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))

  estimates <- vapply(seq_len(n_boot), function(i)
    statistic(sample(values, length(values), replace = TRUE)), numeric(1))
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  c(estimate = statistic(values),
    stats::setNames(stats::quantile(estimates, probs, names = FALSE), c("lower", "upper")))
}

#' Which orthogonal assay could test a candidate.
#'
#' A candidate with no proposed orthogonal test is a number, not a finding.
#' Matching is by keyword against the pathway or module name; anything
#' unmatched is reported as such rather than guessed at.
ASSAY_KEYWORDS <- list(
  "mitochondrial content"  = c("oxidative phosphoryl", "oxphos", "respiratory chain",
                               "electron transport", "mitochond", "atp synth",
                               "tricarboxylic", "aerobic respiration"),
  "ROS"                    = c("reactive oxygen", "superoxide", "oxidative stress",
                               "respiratory burst", "nadph oxidase", "peroxide"),
  "NETs"                   = c("extracellular trap", "chromatin decondens",
                               "histone citrullin", "netosis"),
  "phagocytosis"           = c("phagocyt", "engulf", "opson", "fc receptor",
                               "complement receptor"),
  "glucose uptake"         = c("glycoly", "glucose", "hexose", "carbohydrate metab",
                               "pyruvate"),
  "CBC"                    = c("granulocyte differentiation", "myeloid differentiation",
                               "granulopoiesis", "cell cycle", "proliferat",
                               "myeloid cell development")
)

suggest_orthogonal_assay <- function(candidate, keywords = ASSAY_KEYWORDS) {
  lowered <- tolower(gsub("_", " ", candidate))
  vapply(lowered, function(text) {
    hits <- names(keywords)[vapply(keywords, function(terms)
      any(vapply(terms, grepl, logical(1), x = text, fixed = TRUE)), logical(1))]
    if (!length(hits)) "none proposed -- not a candidate until one is" else hits[1]
  }, character(1), USE.NAMES = FALSE)
}

#' Assemble the pipeline's actual deliverable.
#'
#' One row per candidate, shaped so that a reader cannot mistake it for a
#' results table. Every column answers a question a validation experiment would
#' ask; `orthogonal_assay` is the one that makes the row useful.
build_candidate_table <- function(candidates, cfg) {
  required <- c("candidate", "effect_size", "ci_lower", "ci_upper", "direction",
                "null_percentile", "survives_downsampling", "stage", "n_cells_min")
  missing <- setdiff(required, names(candidates))
  if (length(missing))
    stop("candidate table is missing: ", paste(missing, collapse = ", "))

  out <- annotate_bias(candidates, candidates$direction, cfg)
  out$orthogonal_assay <- suggest_orthogonal_assay(out$candidate)
  out$gate_min_cells <- out$n_cells_min >= cfg$gates$min_cells_per_stratum

  # Ranked by what should be believed, not by significance: something that
  # resists the confound, sits outside the null and survives depth matching
  # comes first, whatever its nominal effect size.
  out$rank_score <- (as.integer(out$against_bias) * 4) +
    (as.integer(out$null_percentile >= cfg$null_model$quantile * 100) * 2) +
    as.integer(out$survives_downsampling %in% TRUE)

  out <- out[order(-out$rank_score, -abs(out$effect_size)), ]
  out <- out[, c("candidate", "stage", "effect_size", "ci_lower", "ci_upper",
                 "direction", "against_bias", "evidence_tier", "null_percentile",
                 "survives_downsampling", "n_cells_min", "gate_min_cells",
                 "orthogonal_assay", "rank_score")]
  attr(out, "caveat") <- sex_contrast_caveat(cfg)
  out
}
