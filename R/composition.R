# ---------------------------------------------------------------------------
# Stage composition: does the *mix* of developmental stages shift with age?
#
# This is the most direct reading of "how does neutrophil development change
# with age": a shift toward immature stages is the classic emergency-
# granulopoiesis signature of an aged marrow, and it is invisible to any
# analysis that compares expression within a stage.
#
# A CAVEAT THAT APPLIES TO EVERY TEST IN THIS FILE. There is one library per
# sex and one hashtag per age, so cells are pseudoreplicates: these tests treat
# each cell as an independent draw, which it is not. The p-values therefore
# answer "is this cell population's composition different?" and not "is the
# composition of mice of this age different?", and they are anticonservative
# for the second question by roughly the square root of the cells per group.
# Read the effect sizes (proportions, log-odds slopes) as the result and the
# p-values as a ranking device.
#
# The age comparison is at least free of batch confounding, since all ages
# share a library. The sex comparison is not: see docs/does-it-answer-the-question.md.
# ---------------------------------------------------------------------------

#' Cell counts and proportions per stage, within each age x sex group.
#'
#' Proportions are within-group (they sum to 1 over stages for each age x sex),
#' which is what makes them comparable when group sizes differ several-fold.
stage_composition <- function(obj, cfg,
                              stage_col = "fine_neu_labels",
                              age_col = "age",
                              sex_col = "sex",
                              age_levels = cfg$analysis$age_levels) {
  meta <- data.frame(
    stage = factor(as.character(obj[[stage_col]][, 1]), levels = cfg$analysis$stage_levels),
    age   = factor(as.character(obj[[age_col]][, 1]), levels = age_levels),
    sex   = factor(as.character(obj[[sex_col]][, 1]), levels = cfg$analysis$sex_levels)
  )
  meta <- meta[stats::complete.cases(meta), ]

  counts <- as.data.frame(table(meta$sex, meta$age, meta$stage),
                          responseName = "n", stringsAsFactors = FALSE)
  names(counts)[1:3] <- c("sex", "age", "stage")

  totals <- stats::aggregate(n ~ sex + age, counts, sum)
  names(totals)[3] <- "n_group"

  out <- merge(counts, totals, by = c("sex", "age"))
  out$proportion <- out$n / out$n_group
  ci <- wilson_interval(out$n, out$n_group)
  out$ci_low <- ci$lower
  out$ci_high <- ci$upper

  out$age <- factor(out$age, levels = age_levels)
  out$stage <- factor(out$stage, levels = cfg$analysis$stage_levels)
  out[order(out$sex, out$age, out$stage), ]
}

#' Wilson score interval for a binomial proportion.
#'
#' Preferred over the normal approximation because stage proportions are often
#' near 0 or 1, where the normal interval runs outside [0, 1].
wilson_interval <- function(successes, total, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- successes / total
  denom <- 1 + z^2 / total
  centre <- (p + z^2 / (2 * total)) / denom
  half <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) / denom
  list(lower = pmax(0, centre - half), upper = pmin(1, centre + half))
}

#' Centred log-ratio of a composition.
#'
#' Proportions are constrained to sum to 1, so one stage rising forces the
#' others down; CLR removes that constraint so stages can be compared without
#' the closure artefact.
clr_transform <- function(proportions, pseudocount = 0.5, n_group = NULL) {
  p <- proportions
  if (!is.null(n_group)) p <- (proportions * n_group + pseudocount) / (n_group + pseudocount)
  p[p <= 0] <- .Machine$double.eps
  log(p) - mean(log(p))
}

#' Does the stage distribution depend on age? One test per sex.
#'
#' Two things are reported per sex:
#'   - an omnibus chi-square of stage x age, i.e. "does the mix change at all";
#'   - per stage, the log-odds slope of that stage against age from a binomial
#'     GLM on the aggregated counts, i.e. "which way and how fast".
#'
#' Age enters the slope model as its rank in `age_levels`, so the slope is per
#' timepoint step, not per month.
composition_age_test <- function(composition, cfg, age_levels = cfg$analysis$age_levels) {
  results <- list()

  for (this_sex in levels(composition$sex)) {
    sub <- composition[composition$sex == this_sex, ]
    if (!nrow(sub)) next

    tab <- stats::xtabs(n ~ stage + age, sub)
    tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
    omnibus <- suppressWarnings(stats::chisq.test(tab))

    log_step(sprintf("%s: stage x age chi-square = %.1f, df = %d, p = %.3g",
                     this_sex, omnibus$statistic, omnibus$parameter, omnibus$p.value))
    if (any(omnibus$expected < 5))
      log_step("  (some expected counts < 5; treat the omnibus p as approximate)")

    per_stage <- do.call(rbind, lapply(levels(sub$stage), function(this_stage) {
      rows <- sub[sub$stage == this_stage, ]
      rows <- rows[order(match(as.character(rows$age), age_levels)), ]
      if (nrow(rows) < 3 || sum(rows$n) < 20) return(NULL)

      age_rank <- match(as.character(rows$age), age_levels)
      fit <- stats::glm(cbind(rows$n, rows$n_group - rows$n) ~ age_rank,
                        family = stats::binomial())
      coefs <- summary(fit)$coefficients

      data.frame(
        sex = this_sex,
        stage = this_stage,
        n_total = sum(rows$n),
        prop_first = rows$proportion[1],
        prop_last = rows$proportion[nrow(rows)],
        log_odds_per_step = unname(coefs["age_rank", "Estimate"]),
        se = unname(coefs["age_rank", "Std. Error"]),
        p_value = unname(coefs["age_rank", "Pr(>|z|)"]),
        stringsAsFactors = FALSE
      )
    }))

    if (!is.null(per_stage)) {
      per_stage$padj <- stats::p.adjust(per_stage$p_value, method = "BH")
      per_stage$fold_change <- per_stage$prop_last / per_stage$prop_first
    }

    results[[this_sex]] <- list(omnibus = omnibus, per_stage = per_stage)
  }

  list(
    per_sex = results,
    trends = do.call(rbind, lapply(results, `[[`, "per_stage"))
  )
}

#' Does the stage distribution differ between the sexes, within each age?
#'
#' Reported for completeness and flagged, not to be believed on its own: with
#' one library per sex, a composition difference between sexes is also a
#' difference between libraries, and dissociation efficiency differs between
#' runs in exactly the way that would move mature-cell fractions.
composition_sex_test <- function(composition, cfg, age_levels = cfg$analysis$age_levels) {
  do.call(rbind, lapply(age_levels, function(this_age) {
    sub <- composition[composition$age == this_age, ]
    tab <- stats::xtabs(n ~ stage + sex, sub)
    tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)

    test <- suppressWarnings(stats::chisq.test(tab))
    data.frame(age = this_age,
               chisq = unname(test$statistic),
               df = unname(test$parameter),
               p_value = test$p.value,
               caveat = "confounded with library: one 10x run per sex",
               stringsAsFactors = FALSE)
  }))
}

#' The composition as a wide, CLR-transformed matrix, one row per age x sex.
#'
#' This is the form to use when comparing the size of a shift across stages,
#' since CLR values are not forced to sum to a constant.
composition_clr_matrix <- function(composition) {
  groups <- unique(composition[, c("sex", "age")])
  out <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
    rows <- composition[composition$sex == groups$sex[i] &
                          composition$age == groups$age[i], ]
    rows <- rows[order(rows$stage), ]
    stats::setNames(clr_transform(rows$proportion, n_group = rows$n_group),
                    as.character(rows$stage))
  }))
  rownames(out) <- paste(groups$sex, groups$age, sep = "_")
  out
}
