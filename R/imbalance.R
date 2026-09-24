# ---------------------------------------------------------------------------
# Population imbalance: which populations are over- or under-represented, and
# in which age or sex.
#
# A proportion table says how big each population is. It does not say whether
# that size is surprising. A stage holding 40% of cells at every age is not
# imbalanced; one holding 40% at 18m and 12% at 3m is, and both look like
# "40%" in a composition table read one row at a time.
#
# The score here is log2(observed / expected), where expected comes from the
# margins of the table being examined -- what the population's share would be
# if it were distributed across ages (or sexes) in the same ratio as every
# other population. Zero means exactly as expected. +1 means twice its share.
#
# Two contrasts, and they are NOT equivalent:
#
#   across ages WITHIN a sex   all ages share a library, so this is clean.
#   across sexes WITHIN an age  sex is confounded with library, run, probe
#                               barcode and staining batch; every such row
#                               carries that caveat in the table itself.
# ---------------------------------------------------------------------------

#' Counts of each population per age and sex.
#'
#' The population column is generalised (`stage`, `fine_label_readable`, ...)
#' so the same machinery serves maturation stages and the progenitor
#' compartment upstream of them.
population_counts <- function(obj, cfg, population_col = "stage",
                              levels = NULL,
                              age_col = "age", sex_col = "sex",
                              age_levels = cfg$analysis$age_levels) {
  values <- as.character(obj[[population_col]][, 1])
  levels <- levels %||% sort(unique(values[!is.na(values)]))

  meta <- data.frame(
    population = factor(values, levels = levels),
    age = factor(as.character(obj[[age_col]][, 1]), levels = age_levels),
    sex = factor(as.character(obj[[sex_col]][, 1]), levels = cfg$analysis$sex_levels)
  )
  meta <- meta[stats::complete.cases(meta), ]

  out <- as.data.frame(table(meta$sex, meta$age, meta$population),
                       responseName = "n", stringsAsFactors = FALSE)
  names(out)[1:3] <- c("sex", "age", "population")
  out$age <- factor(out$age, levels = age_levels)
  out$population <- factor(out$population, levels = levels)
  out[order(out$sex, out$age, out$population), ]
}

#' log2(observed / expected) and standardised residuals for one table.
#'
#' `tab` is populations (rows) x groups (columns). Expected is the outer
#' product of the margins over the grand total: the count each cell would hold
#' if population and group were independent.
#'
#' The standardised Pearson residual divides by the full variance rather than
#' sqrt(expected), so residuals from cells with very different margins are
#' comparable -- a plain Pearson residual is systematically larger for small
#' expected counts and makes rare populations look more imbalanced than they
#' are.
table_imbalance <- function(tab, pseudocount = 0.5) {
  total <- sum(tab)
  row_m <- rowSums(tab)
  col_m <- colSums(tab)
  expected <- outer(row_m, col_m) / total

  # The standardised residual is undefined where its denominator vanishes:
  # a population with no cells anywhere in the block (expected 0), or one
  # holding every cell in it (row margin == total, so 1 - row_m/total == 0).
  # Both arise on real progenitor data -- rare lineages absent from a library,
  # and blocks left with a single population after the cell-count filter --
  # and both produce NaN, which reads downstream as a computed value. NA says
  # "not defined here" instead.
  denom <- expected * outer(1 - row_m / total, 1 - col_m / total)
  resid <- (tab - expected) / sqrt(denom)
  resid[!is.finite(resid)] <- NA_real_

  # The pseudocount keeps an absent population finite. With counts in the
  # hundreds it shifts log2 by well under 0.01; in a cell holding zero it is
  # the difference between a number and -Inf.
  log2_ratio <- log2((tab + pseudocount) / (expected + pseudocount))

  list(expected = expected, log2_ratio = log2_ratio, std_residual = resid)
}

#' Bootstrap interval for log2(observed / expected).
#'
#' Cells are resampled within each group from the observed proportions, and
#' the whole table -- margins included -- is rebuilt each iteration, because
#' the expected counts depend on the margins and would otherwise be treated as
#' fixed when they are estimated from the same cells.
#'
#' This is a cell-level interval. It does NOT capture between-animal variance,
#' which is unmeasurable in this design, so it understates the true
#' uncertainty. It says how much of the score survives cell sampling alone.
bootstrap_imbalance <- function(tab, n_boot = 1000, seed = 42, conf = 0.95) {
  set.seed(seed)
  col_totals <- colSums(tab)

  draws <- vapply(seq_len(n_boot), function(i) {
    resampled <- vapply(seq_along(col_totals), function(j) {
      if (col_totals[j] == 0) return(rep(0, nrow(tab)))
      probs <- tab[, j] / col_totals[j]
      as.numeric(stats::rmultinom(1, col_totals[j], probs))
    }, numeric(nrow(tab)))
    # vapply returns a bare vector when FUN.VALUE has length 1, i.e. when the
    # block holds a single population, and dimnames<- then fails on a
    # non-array. Rebuild the matrix rather than relying on the shape.
    resampled <- matrix(resampled, nrow = nrow(tab), ncol = length(col_totals),
                        dimnames = dimnames(tab))
    as.numeric(table_imbalance(resampled)$log2_ratio)
  }, numeric(length(tab)))

  alpha <- (1 - conf) / 2
  list(lower = matrix(apply(draws, 1, stats::quantile, alpha, na.rm = TRUE),
                      nrow = nrow(tab), dimnames = dimnames(tab)),
       upper = matrix(apply(draws, 1, stats::quantile, 1 - alpha, na.rm = TRUE),
                      nrow = nrow(tab), dimnames = dimnames(tab)))
}

#' Imbalance across one factor, held within the levels of another.
#'
#' @param across the column whose levels are compared ("age" or "sex").
#' @param within the column whose levels are each analysed separately.
imbalance_scores <- function(counts, across = "age", within = "sex",
                             n_boot = 1000, seed = 42,
                             caveat = "") {
  blocks <- split(counts, counts[[within]], drop = TRUE)

  out <- lapply(names(blocks), function(level) {
    block <- blocks[[level]]
    tab <- stats::xtabs(stats::reformulate(c("population", across), "n"), block)
    tab <- matrix(as.numeric(tab), nrow = nrow(tab),
                  dimnames = dimnames(tab))

    # A group holding no cells is not a group with an expected share of zero;
    # it is a group this contrast has nothing to say about. Emitting zeros for
    # it would read as "exactly as expected".
    empty <- colSums(tab) == 0
    if (any(empty)) {
      log_step("  ", level, ": no cells in ",
               paste(colnames(tab)[empty], collapse = ", "),
               " -- excluded from this contrast")
      tab <- tab[, !empty, drop = FALSE]
    }
    if (ncol(tab) < 2 || sum(tab) == 0) return(NULL)

    # Imbalance is relative: with one population left there is nothing for it
    # to be imbalanced against, and every score is identically zero. Saying so
    # beats reporting a table of zeros.
    if (nrow(tab) < 2) {
      log_step("  ", level, ": only one population (",
               rownames(tab)[1], ") -- imbalance is undefined with nothing ",
               "to compare against")
      return(NULL)
    }

    stats_tbl <- table_imbalance(tab)
    boot <- bootstrap_imbalance(tab, n_boot = n_boot, seed = seed)
    grid <- expand.grid(population = rownames(tab), group = colnames(tab),
                        stringsAsFactors = FALSE)

    data.frame(
      within_level = level,
      group = grid$group,
      population = grid$population,
      n = as.numeric(tab),
      n_in_group = rep(colSums(tab), each = nrow(tab)),
      observed_prop = as.numeric(tab) / rep(colSums(tab), each = nrow(tab)),
      expected_prop = as.numeric(stats_tbl$expected) /
        rep(colSums(tab), each = nrow(tab)),
      log2_obs_exp = as.numeric(stats_tbl$log2_ratio),
      ci_low = as.numeric(boot$lower),
      ci_high = as.numeric(boot$upper),
      std_residual = as.numeric(stats_tbl$std_residual),
      row.names = NULL, stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, Filter(Negate(is.null), out))
  if (is.null(out)) return(NULL)

  names(out)[names(out) == "within_level"] <- within
  names(out)[names(out) == "group"] <- across
  # A score whose interval spans zero is not an imbalance this data can see.
  out$excludes_zero <- (out$ci_low > 0) | (out$ci_high < 0)
  out$contrast <- paste0("across ", across, " within ", within)
  out$caveat <- caveat
  out[order(-abs(out$log2_obs_exp)), ]
}

#' Both contrasts, with the confounded one labelled as such.
#'
#' Reported together deliberately. The age contrast is the interpretable one
#' and the sex contrast is the one a reader will reach for first, so the
#' caveat travels in the same table rather than in a paragraph elsewhere.
imbalance_report <- function(counts, cfg, n_boot = NULL) {
  n_boot <- n_boot %||% cfg$imbalance$n_boot %||% 1000

  by_age <- imbalance_scores(counts, across = "age", within = "sex",
                             n_boot = n_boot, seed = cfg$seed,
                             caveat = "ages share a library: interpretable")
  by_sex <- imbalance_scores(counts, across = "sex", within = "age",
                             n_boot = n_boot, seed = cfg$seed,
                             caveat = sex_contrast_caveat(cfg))

  harmonise <- function(df, across) {
    if (is.null(df)) return(NULL)
    df$across_level <- as.character(df[[across]])
    df$within_level <- as.character(df[[setdiff(c("age", "sex"), across)]])
    df$across_factor <- across
    df[, c("contrast", "across_factor", "within_level", "across_level",
           "population", "n", "n_in_group", "observed_prop", "expected_prop",
           "log2_obs_exp", "ci_low", "ci_high", "std_residual",
           "excludes_zero", "caveat")]
  }

  do.call(rbind, Filter(Negate(is.null),
                        list(harmonise(by_age, "age"), harmonise(by_sex, "sex"))))
}

#' Print the populations furthest from expectation, age contrast first.
report_imbalance <- function(scores, top_n = 12) {
  if (is.null(scores) || !nrow(scores)) {
    log_step("no imbalance scores to report")
    return(invisible(NULL))
  }
  for (contrast in unique(scores$contrast)) {
    block <- scores[scores$contrast == contrast, ]
    log_step(contrast, " -- ", block$caveat[1])
    shown <- utils::head(block[block$excludes_zero, ], top_n)
    if (!nrow(shown)) {
      log_step("  no population's interval excludes zero")
      next
    }
    print(shown[, c("within_level", "across_level", "population", "n",
                    "observed_prop", "expected_prop", "log2_obs_exp",
                    "ci_low", "ci_high")], row.names = FALSE)
  }
  invisible(scores)
}
