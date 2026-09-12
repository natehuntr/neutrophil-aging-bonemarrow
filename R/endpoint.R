# ---------------------------------------------------------------------------
# Endpoint contrast, then interpolation.
#
# Two questions asked in order, because they need different amounts of data:
#
#   DISCOVERY  3m vs 18m within one stage and one sex. The extremes carry the
#              largest difference and the most cells, so this is where genes
#              are found. Calibrated against a permutation of the two labels,
#              so the threshold comes from the stratum rather than a constant.
#
#   SHAPE      Do 9m and 12m sit between the endpoints, for those same genes,
#              in those same cells? A gene that moves monotonically across all
#              four timepoints is a better aging hypothesis than one that
#              differs at the extremes and does something else in between.
#
# The shape check is NOT independent evidence: it reuses the cells the
# discovery contrast came from, and the intermediate strata are often small.
# It reorders a candidate list; it does not confirm anything.
# ---------------------------------------------------------------------------

#' Mean normalised expression per group, genes x groups.
group_means <- function(expr, groups, levels = NULL) {
  levels <- levels %||% sort(unique(groups[!is.na(groups)]))
  out <- vapply(levels, function(g) {
    idx <- which(!is.na(groups) & groups == g)
    if (!length(idx)) return(rep(NA_real_, nrow(expr)))
    Matrix::rowMeans(expr[, idx, drop = FALSE])
  }, numeric(nrow(expr)))
  rownames(out) <- rownames(expr)
  out
}

#' Per-gene difference in mean normalised expression between two groups.
endpoint_effect <- function(expr, groups, reference, endpoint) {
  means <- group_means(expr, groups, c(reference, endpoint))
  means[, 2] - means[, 1]
}

#' Null for the endpoint difference, from shuffling the two labels.
#'
#' The statistic is the MAXIMUM absolute difference across genes, so the 95th
#' percentile is a family-wise threshold for the stratum: it answers "how big
#' does the biggest difference get when the labels mean nothing here".
endpoint_null <- function(expr, groups, reference, endpoint,
                          n_perm = 100, seed = 42) {
  keep <- which(!is.na(groups) & groups %in% c(reference, endpoint))
  sub <- expr[, keep, drop = FALSE]
  labels <- groups[keep]

  vapply(seq_len(n_perm), function(i) {
    set.seed(seed + i)
    max(abs(endpoint_effect(sub, sample(labels), reference, endpoint)),
        na.rm = TRUE)
  }, numeric(1))
}

#' Where each intermediate age sits on the 3m -> 18m line, per gene.
#'
#' `fraction_of_endpoint` is (mean_age - mean_reference) / (mean_endpoint -
#' mean_reference). Between 0 and 1 means the timepoint interpolates; below 0
#' means it moves the other way first; above 1 means it overshoots and comes
#' back. The denominator is the discovery effect, which is non-zero by
#' construction for every gene reaching this function.
interpolation_profile <- function(expr, ages, genes, cfg,
                                  reference = "3m", endpoint = "18m",
                                  age_levels = cfg$analysis$age_levels) {
  expr <- expr[genes, , drop = FALSE]
  means <- group_means(expr, ages, age_levels)
  present <- colnames(means)[colSums(!is.na(means)) > 0]

  delta <- means[, endpoint] - means[, reference]
  middles <- setdiff(intersect(age_levels, present), c(reference, endpoint))

  out <- data.frame(gene = genes, stringsAsFactors = FALSE)
  for (age in age_levels)
    if (age %in% present) out[[paste0("mean_", age)]] <- unname(means[, age])

  fractions <- list()
  for (age in middles) {
    frac <- (means[, age] - means[, reference]) / delta
    out[[paste0("frac_", age)]] <- unname(frac)
    fractions[[age]] <- frac
  }

  if (!length(fractions)) {
    out$shape <- "no intermediate ages"
    return(out)
  }

  frac_mat <- do.call(cbind, fractions)
  # Ordered by age, an interpolating gene has fractions that rise from 0 to 1
  # without leaving [0, 1] and without going backwards.
  in_range <- apply(frac_mat, 1, function(f) all(is.finite(f) & f >= 0 & f <= 1))
  ordered <- if (ncol(frac_mat) < 2) rep(TRUE, nrow(frac_mat)) else
    apply(frac_mat, 1, function(f) all(is.finite(f)) && !is.unsorted(f))

  out$shape <- ifelse(in_range & ordered, "monotonic",
               ifelse(in_range, "between endpoints, out of order",
               ifelse(apply(frac_mat, 1, function(f) any(is.finite(f) & f < 0)),
                      "reverses before endpoint", "overshoots endpoint")))
  out
}

#' Discovery contrast plus shape check for one stage x sex stratum.
endpoint_contrast <- function(obj, cfg, label = "", assay = "RNA",
                              age_col = "age",
                              reference = cfg$endpoint$reference,
                              endpoint = cfg$endpoint$endpoint,
                              n_perm = cfg$endpoint$n_perm %||% 100) {
  ages <- as.character(obj[[age_col]][, 1])
  n_by_age <- table(ages)
  endpoints_n <- n_by_age[intersect(c(reference, endpoint), names(n_by_age))]

  if (length(endpoints_n) < 2 || !stratum_is_usable(endpoints_n, label, cfg))
    return(NULL)

  expr <- Seurat::GetAssayData(obj, assay = assay, layer = "data")
  # A gene absent from the stratum cannot differ across it, and leaving it in
  # only inflates the multiple-testing surface the permutation null covers.
  expressed <- Matrix::rowSums(expr > 0) >= cfg$endpoint$min_cells_detected
  expr <- expr[expressed, , drop = FALSE]
  log_step(sprintf("  %s: %d genes detected in >= %d cells", label, nrow(expr),
                   cfg$endpoint$min_cells_detected))

  effect <- endpoint_effect(expr, ages, reference, endpoint)
  null_max <- endpoint_null(expr, ages, reference, endpoint, n_perm, cfg$seed)
  threshold <- stats::quantile(null_max, cfg$endpoint$quantile %||% 0.95,
                               names = FALSE, na.rm = TRUE)

  hits <- names(which(abs(effect) > threshold))
  log_step(sprintf("  %s: permutation threshold %.4f; %d of %d genes exceed it",
                   label, threshold, length(hits), length(effect)))
  if (!length(hits)) return(list(threshold = threshold, null_max = null_max,
                                 genes = NULL, n_by_age = n_by_age))

  shape <- interpolation_profile(expr, ages, hits, cfg,
                                 reference = reference, endpoint = endpoint)
  shape$effect <- unname(effect[hits])
  shape$direction <- ifelse(shape$effect > 0,
                            paste0("up at ", endpoint), paste0("down at ", endpoint))
  shape$excess_over_null <- abs(shape$effect) - threshold

  # Monotonic first: same discovery evidence, better-behaved across the ages
  # the contrast never used.
  shape <- shape[order(shape$shape != "monotonic", -shape$excess_over_null), ]
  log_step("  ", label, ": shape of those genes across 9m/12m -- ",
           paste(sprintf("%s %d", names(table(shape$shape)), table(shape$shape)),
                 collapse = ", "))

  list(threshold = threshold, null_max = null_max, genes = shape,
       n_by_age = n_by_age)
}
