# ---------------------------------------------------------------------------
# The within-library empirical null.
#
# Sex is confounded with library, so sex labels cannot be permuted: shuffling
# them shuffles the library too, and the "null" would still contain the whole
# confound. The alternative is to ask what this pipeline produces from a
# comparison that is known to be empty -- a random split of cells inside ONE
# library, matched to the real contrast on stage composition and cell number.
#
# The distribution of hits across many such splits is the false-positive floor
# for this pipeline on this data. A real effect has to clear it.
#
# WHAT THIS IS NOT. It is a floor, not a significance test. Because each
# library pools several animals, a random cell split samples across mice, so it
# captures cell-level sampling plus some between-animal variation -- more than
# a one-mouse design would give, but still less than a true between-animal
# null, which cannot be built here because cells cannot be traced to animals.
# Describing it as anything more invites exactly the criticism it exists to
# pre-empt.
# ---------------------------------------------------------------------------

#' Split cells within one library into two pseudo-groups.
#'
#' Matched on stage composition and on the size of the smaller real group, so
#' the split poses the same question, at the same power, as the real contrast.
matched_split <- function(meta, stage_col = "stage", target_n = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stages <- as.character(meta[[stage_col]])
  assignment <- rep(NA_character_, nrow(meta))

  for (stage in unique(stages[!is.na(stages)])) {
    idx <- which(stages == stage)
    # Cap each stage at twice the target so the two halves can both reach it.
    if (!is.null(target_n) && length(idx) > 2 * target_n)
      idx <- sample(idx, 2 * target_n)
    shuffled <- sample(idx)
    half <- length(shuffled) %/% 2
    if (half < 1) next
    assignment[shuffled[seq_len(half)]] <- "split_a"
    assignment[shuffled[(half + 1):(2 * half)]] <- "split_b"
  }
  assignment
}

#' Run an analysis repeatedly on empty comparisons.
#'
#' @param analysis a function of (object, assignment) returning a single number
#'   -- the count of hits, or any other headline quantity being calibrated.
#' @return the null values, plus the floor at the configured quantile.
within_library_null <- function(obj, cfg, analysis,
                                stage_col = "stage",
                                sex_col = "sex",
                                n_splits = NULL,
                                target_n = NULL) {
  n_splits <- n_splits %||% cfg$null_model$n_splits

  # Build the null in the deeper library: it is the one that produces more
  # spurious detection, so its floor is the conservative one.
  library_sex <- if (identical(cfg$null_model$library, "deeper"))
    cfg$design$deeper_library else cfg$null_model$library
  cells <- colnames(obj)[as.character(obj[[sex_col]][, 1]) == library_sex]
  if (length(cells) < 4 * (target_n %||% 50))
    warning("the ", library_sex, " library has ", length(cells),
            " cells; the null will be noisy")
  sub <- subset(obj, cells = cells)

  log_step(sprintf("within-library null: %d splits inside the %s library (%d cells)",
                   n_splits, library_sex, ncol(sub)))

  values <- vapply(seq_len(n_splits), function(i) {
    assignment <- matched_split(sub@meta.data, stage_col, target_n, seed = cfg$seed + i)
    ok <- !is.na(assignment)
    if (sum(ok) < 2) return(NA_real_)
    tryCatch(as.numeric(analysis(subset(sub, cells = colnames(sub)[ok]), assignment[ok])),
             error = function(e) NA_real_)
  }, numeric(1))

  finite <- values[is.finite(values)]
  floor_value <- if (length(finite)) stats::quantile(finite, cfg$null_model$quantile,
                                                     names = FALSE) else NA_real_

  log_step(sprintf("  null: median %.1f, %g%% floor %.1f (%d/%d splits usable)",
                   stats::median(finite), 100 * cfg$null_model$quantile,
                   floor_value, length(finite), n_splits))

  structure(list(values = values, floor = floor_value,
                 quantile = cfg$null_model$quantile,
                 library = library_sex, n_splits = n_splits),
            class = "within_library_null")
}

#' Where an observed value sits in the null distribution.
#'
#' Reported as a percentile rather than a p-value, deliberately: it says how
#' unusual the value is for this pipeline on this data, which is a weaker and
#' more honest claim than a significance test.
null_percentile <- function(observed, null_result) {
  values <- null_result$values[is.finite(null_result$values)]
  if (!length(values) || is.na(observed)) return(NA_real_)
  100 * mean(values <= observed)
}

#' Compare an observed count against the floor, in the form a table wants.
compare_to_null <- function(observed, null_result, label = "") {
  pct <- null_percentile(observed, null_result)
  data.frame(
    comparison = label,
    n_observed = observed,
    null_median = stats::median(null_result$values, na.rm = TRUE),
    null_floor = null_result$floor,
    null_percentile = pct,
    exceeds_floor = !is.na(observed) && !is.na(null_result$floor) &&
      observed > null_result$floor,
    stringsAsFactors = FALSE
  )
}

#' A sentence for the report, so the floor is never quoted as a p-value.
describe_null <- function(null_result, cfg) {
  sprintf(paste0("False-positive floor from %d random splits within the %s library, ",
                 "matched on stage composition and cell number. This is the number of ",
                 "hits this pipeline yields from a comparison known to be empty; a real ",
                 "effect must exceed it. It is a floor, not a significance test: the ",
                 "split samples across the %d animals pooled per hashtag, so it is still ",
                 "anticonservative relative to a between-animal null, which this design ",
                 "cannot provide."),
          null_result$n_splits, null_result$library, cfg$design$mice_per_hashtag)
}
