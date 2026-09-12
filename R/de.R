# ---------------------------------------------------------------------------
# Differential expression and age-trend analysis.
#
# Three independent approaches are kept, because they answer different
# questions and were all used in the original analysis:
#
#   1. age_trend_clusters()  - Spearman rho vs age per gene, then cluster the
#                              z-scored age profiles into trajectory shapes.
#   2. age_changing()        - pairwise Wilcoxon between ages with a
#                              permutation-calibrated |log2FC| threshold.
#   3. glm_lrt() + calibrate_threshold() - negative-binomial GLM likelihood
#                              ratio tests with an empirical-FDR cutoff.
# ---------------------------------------------------------------------------

# --- 1. Spearman age trends ------------------------------------------------

#' Spearman correlation of every row of `expr` against `y`, vectorised.
#'
#' apply(expr, 1, cor.test) coerces the whole sparse matrix to dense and calls
#' cor.test tens of thousands of times, which dominates the runtime of every
#' age-trend step. Spearman is Pearson on ranks, so ranking once and calling
#' cor() gives identical rho for a fraction of the cost.
#'
#' The p-value is the same asymptotic t approximation cor.test() uses for
#' method = "spearman", exact = FALSE, so results match that call rather than
#' merely resembling it. Rows with no variance get NA, as they would there.
spearman_rows <- function(expr, y) {
  n <- length(y)
  stopifnot(ncol(expr) == n, n > 3)

  y_rank <- rank(y)
  expr_rank <- t(apply(as.matrix(expr), 1, rank))

  rho <- as.numeric(suppressWarnings(stats::cor(t(expr_rank), y_rank)))
  names(rho) <- rownames(expr)

  # t = rho * sqrt((n - 2) / (1 - rho^2)), two-sided on n - 2 df.
  denom <- 1 - rho^2
  t_stat <- ifelse(denom <= 0, NA_real_, rho * sqrt((n - 2) / denom))
  pval <- 2 * stats::pt(abs(t_stat), df = n - 2, lower.tail = FALSE)

  data.frame(gene = names(rho), rho = unname(rho), pval = pval,
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Per-gene Spearman correlation of expression against age.
#'
#' Age is converted to a rank via `age_levels`, so the test is for a monotonic
#' trend across the ordered timepoints, not a linear fit to months.
age_trend_test <- function(obj, age_levels, age_col = "age", assay = "RNA") {
  age_char <- as.character(obj[[age_col]][, 1])
  extra <- setdiff(unique(age_char), age_levels)
  if (length(extra))
    stop("ages present that are not in age_levels: ", paste(extra, collapse = ", "))

  age_numeric <- as.numeric(factor(age_char, levels = age_levels))
  expr <- Seurat::GetAssayData(obj, assay = assay, layer = "data")

  df <- spearman_rows(expr, age_numeric)
  df$padj <- stats::p.adjust(df$pval, method = "BH")
  rownames(df) <- df$gene
  df
}

#' Signed -log10(p) ranking vector for preranked GSEA.
gsea_ranking <- function(trend_df, min_p = 1e-300) {
  p_floor <- ifelse(is.na(trend_df$pval), NA_real_, pmax(trend_df$pval, min_p))
  score <- sign(trend_df$rho) * (-log10(p_floor))
  score[is.na(score)] <- 0

  ord <- order(-score, -ifelse(is.na(trend_df$rho), 0, trend_df$rho), trend_df$gene)
  ranks <- stats::setNames(score[ord], trend_df$gene[ord])
  ranks[!duplicated(names(ranks))]
}

#' Mean expression per age, z-scored across ages, in age order.
#'
#' AverageExpression prefixes group names that start with a digit with "g"
#' ("3m" -> "g3m"), so columns are matched on the trailing age string rather
#' than on the prefix.
age_profile_matrix <- function(obj, genes, age_levels, age_col = "age", assay = "RNA") {
  avg <- Seurat::AverageExpression(obj, group.by = age_col,
                                   assays = assay, layer = "data")[[assay]]
  col_age <- sub("^g", "", colnames(avg))
  age_order <- match(age_levels, col_age)
  if (anyNA(age_order))
    stop("could not match AverageExpression columns to age_levels; columns were: ",
         paste(colnames(avg), collapse = ", "))

  avg <- avg[, age_order, drop = FALSE]
  colnames(avg) <- age_levels
  avg <- avg[intersect(genes, rownames(avg)), , drop = FALSE]

  z <- t(scale(t(avg)))
  z[stats::complete.cases(z), , drop = FALSE]   # drop zero-variance genes
}

#' Name a cluster from its own centroid: e.g. "high_middle_low".
shape_name <- function(centroid, tol) {
  paste(ifelse(centroid > tol, "high",
               ifelse(centroid < -tol, "low", "middle")), collapse = "_")
}

#' Age-trend testing plus trajectory-shape clustering.
#'
#' Returns the genome-wide trend table and GSEA ranking regardless; the
#' clustering is skipped when `do_cluster` is FALSE or too few genes pass.
age_trend_clusters <- function(obj, cfg,
                               age_col = "age",
                               age_levels = cfg$analysis$age_levels,
                               rho_cutoff = cfg$age_trend$rho_cutoff,
                               n_clusters = cfg$age_trend$n_clusters,
                               tol = cfg$age_trend$shape_tolerance,
                               assay = "RNA",
                               do_cluster = TRUE,
                               plot = TRUE) {

  trend_df <- age_trend_test(obj, age_levels, age_col, assay)
  ranks <- gsea_ranking(trend_df)
  log_step(sprintf("genome-wide ranking: %d genes", length(ranks)))

  changing <- trend_df[!is.na(trend_df$rho) & abs(trend_df$rho) > rho_cutoff, ]
  changing <- changing[order(changing$padj), ]
  log_step(sprintf("%d of %d genes pass |rho| > %.2f",
                   nrow(changing), nrow(trend_df), rho_cutoff))

  empty <- list(gsea_rank = ranks, trend_df = trend_df, changing_genes = changing,
                z_scored = NULL, heatmap = NULL, centroids = NULL,
                gene_cluster_df = NULL)
  if (!do_cluster) return(empty)
  if (nrow(changing) < n_clusters) {
    warning("fewer changing genes than n_clusters -- returning ranking only")
    return(empty)
  }

  z <- age_profile_matrix(obj, changing$gene, age_levels, age_col, assay)
  hm <- pheatmap::pheatmap(z, cluster_cols = FALSE, cutree_rows = n_clusters,
                           show_rownames = TRUE, silent = !plot)
  row_clusters <- stats::cutree(hm$tree_row, k = n_clusters)

  centroids <- sapply(sort(unique(row_clusters)), function(k)
    colMeans(z[row_clusters == k, age_levels, drop = FALSE]))
  colnames(centroids) <- sort(unique(row_clusters))
  shapes <- apply(centroids, 2, shape_name, tol = tol)

  if (any(duplicated(shapes)))
    warning("two clusters share a shape name -- consider a smaller tol or fewer clusters")

  gene_cluster_df <- data.frame(
    gene = names(row_clusters),
    cluster_number = unname(row_clusters),
    cluster_name = factor(shapes[as.character(row_clusters)], levels = unique(shapes)),
    row.names = NULL
  )

  c(empty[c("gsea_rank", "trend_df", "changing_genes")],
    list(z_scored = z, heatmap = hm, centroids = centroids,
         gene_cluster_df = gene_cluster_df))
}

#' Join the cluster assignment onto the trend statistics for reporting.
age_trend_table <- function(result) {
  if (is.null(result$gene_cluster_df)) return(result$changing_genes)
  merge(result$gene_cluster_df, result$changing_genes, by = "gene", sort = FALSE)
}

# --- 2. Permutation-calibrated pairwise DE ---------------------------------

#' Genes detected in enough cells of *every* age group.
#'
#' Requiring detection in all groups is what guarantees each gene has all
#' three pairwise comparisons, and is what stops a gene seen in 3 of 26 cells
#' from producing a 256-fold "change".
detectable_genes <- function(obj, cfg, age_levels, age_col = "age") {
  cnt <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
  ages <- as.character(obj[[age_col]][, 1])

  ok <- vapply(age_levels, function(a) {
    idx <- ages == a
    need <- max(cfg$permutation_de$min_cells,
                ceiling(cfg$permutation_de$min_frac * sum(idx)))
    Matrix::rowSums(cnt[, idx, drop = FALSE] > 0) >= need
  }, logical(nrow(cnt)))

  keep <- rownames(cnt)[rowSums(ok) == length(age_levels)]
  log_step(sprintf("  detection filter: %d of %d genes retained", length(keep), nrow(cnt)))
  keep
}

#' Antibodies whose median normalised signal clears the floor in every age.
adequate_adt <- function(obj, cfg, age_levels, age_col = "age") {
  dat <- Seurat::GetAssayData(obj, assay = "ADT", layer = "data")
  ages <- as.character(obj[[age_col]][, 1])

  ok <- vapply(age_levels, function(a) {
    idx <- ages == a
    apply(dat[, idx, drop = FALSE], 1, stats::median) > cfg$permutation_de$adt_min_median
  }, logical(nrow(dat)))

  keep <- rownames(dat)[rowSums(ok) == length(age_levels)]
  log_step(sprintf("  ADT filter: %d of %d antibodies retained", length(keep), nrow(dat)))
  keep
}

age_pairs <- function(age_levels) {
  utils::combn(age_levels, 2, simplify = FALSE)   # older listed second
}

#' All pairwise age comparisons, on real or shuffled age labels.
#'
#' `features` is fixed from the real data so permuted runs test exactly the
#' same gene set; without that the null is not calibrated like for like.
#' logfc.threshold and min.pct are 0 because filtering already happened.
run_pairwise_ages <- function(obj, features, age_levels, assay = "RNA",
                              permute = FALSE, seed = 1, age_col = "age") {
  ages <- as.character(obj[[age_col]][, 1])
  if (permute) { set.seed(seed); ages <- sample(ages) }
  obj$.age_use <- factor(ages, levels = age_levels)
  Seurat::Idents(obj) <- ".age_use"

  purrr::map_dfr(age_pairs(age_levels), function(p) {
    markers <- Seurat::FindMarkers(obj, ident.1 = p[2], ident.2 = p[1],
                                   assay = assay, features = features,
                                   logfc.threshold = 0, min.pct = 0,
                                   test.use = "wilcox", verbose = FALSE)
    data.frame(
      feature    = rownames(markers),
      comparison = paste0(p[1], "_to_", p[2]),
      log2FC     = markers$avg_log2FC,
      pct_diff   = if (all(c("pct.1", "pct.2") %in% names(markers)))
                     markers$pct.1 - markers$pct.2 else NA_real_,
      row.names  = NULL,
      stringsAsFactors = FALSE
    )
  })
}

#' Per-comparison |log2FC| threshold from shuffled age labels.
calibrate_pairwise <- function(obj, features, cfg, age_levels, assay = "RNA") {
  purrr::map_dfr(seq_len(cfg$permutation_de$n_perm), function(i)
    run_pairwise_ages(obj, features, age_levels, assay, permute = TRUE, seed = i)) |>
    dplyr::group_by(.data$comparison) |>
    dplyr::summarise(threshold = stats::quantile(abs(.data$log2FC),
                                                 cfg$permutation_de$quantile,
                                                 na.rm = TRUE),
                     .groups = "drop")
}

#' max() that returns NA instead of -Inf when everything is missing.
max_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) NA_real_ else max(x)
}

#' Classify the 3-point trajectory shape from the early and net changes.
#'
#' Positions are centred the same way a z-score would be, so the labels mean
#' the same thing as the ones from age_trend_clusters().
classify_shape <- function(early, net, thr_early, thr_net, age_levels, tol = NULL) {
  if (is.null(tol)) tol <- mean(c(thr_early, thr_net)) / 2

  positions <- cbind(0, early, net)
  colnames(positions) <- age_levels
  positions <- positions - rowMeans(positions)

  labels <- ifelse(positions > tol, "high", ifelse(positions < -tol, "low", "middle"))
  out <- apply(labels, 1, paste, collapse = "_")
  out[is.na(early) | is.na(net)] <- NA_character_
  out
}

#' Age-changing features with permutation-calibrated thresholds.
#'
#' Works on either assay: RNA uses the detection filter, ADT the median-signal
#' filter. The original notebook had two copies of this that overwrote each
#' other in the global environment; this is the single parameterised version.
age_changing <- function(obj, label, cfg, age_levels = cfg$analysis$age_levels,
                         assay = c("RNA", "ADT")) {
  assay <- match.arg(assay)
  log_step("=== ", label, " / ", assay, " ===")

  features <- if (assay == "RNA") detectable_genes(obj, cfg, age_levels)
              else adequate_adt(obj, cfg, age_levels)
  min_features <- if (assay == "RNA") 50 else 5
  if (length(features) < min_features)
    warning(label, ": only ", length(features), " features pass the filter -- ",
            "this arm is likely underpowered")

  observed <- run_pairwise_ages(obj, features, age_levels, assay)
  thresholds <- calibrate_pairwise(obj, features, cfg, age_levels, assay)
  log_step("  permutation thresholds (|log2FC|):")
  print(as.data.frame(thresholds), row.names = FALSE)

  cmp <- function(from, to) paste0(from, "_to_", to)
  n <- length(age_levels)
  thr <- function(name) thresholds$threshold[thresholds$comparison == name]
  thr_early <- thr(cmp(age_levels[1], age_levels[2]))
  thr_late  <- thr(cmp(age_levels[2], age_levels[n]))
  thr_net   <- thr(cmp(age_levels[1], age_levels[n]))

  wide <- observed |>
    dplyr::select("feature", "comparison", "log2FC") |>
    tidyr::pivot_wider(names_from = "comparison", values_from = "log2FC") |>
    as.data.frame()
  # Rename by position rather than with tidyselect, so the column names stay
  # driven by age_levels instead of being hard-coded.
  names(wide)[match(c(cmp(age_levels[1], age_levels[2]),
                      cmp(age_levels[2], age_levels[n]),
                      cmp(age_levels[1], age_levels[n])), names(wide))] <-
    c("early", "late", "net")

  pct <- observed |>
    dplyr::group_by(.data$feature) |>
    dplyr::summarise(max_abs_pct_diff = max_or_na(abs(.data$pct_diff)),
                     .groups = "drop")

  res <- wide |>
    dplyr::left_join(pct, by = "feature") |>
    dplyr::mutate(
      max_abs_log2FC = pmax(abs(.data$early), abs(.data$late), abs(.data$net)),
      passes = abs(.data$early) > thr_early | abs(.data$late) > thr_late |
               abs(.data$net) > thr_net,
      shape = classify_shape(.data$early, .data$net, thr_early, thr_net, age_levels),
      sex = label
    ) |>
    dplyr::filter(.data$passes) |>
    dplyr::arrange(dplyr::desc(.data$max_abs_log2FC))

  # Sparsity artefacts show up as absurd fold changes; say so rather than
  # letting them through silently.
  med <- stats::median(res$max_abs_log2FC)
  log_step(sprintf("  %d features | median |log2FC| = %.2f | %d NA shapes",
                   nrow(res), med, sum(is.na(res$shape))))
  if (assay == "RNA" && !is.na(med) && med > 2)
    warning(label, ": median |log2FC| = ", round(med, 2),
            " -- still sparsity-driven; raise permutation_de.min_cells")

  res
}

#' Join the two sexes' age-changing tables and flag matching shapes.
compare_sexes <- function(male_res, female_res) {
  dplyr::inner_join(male_res, female_res, by = "feature", suffix = c("_M", "_F")) |>
    dplyr::select("feature", "shape_M", "shape_F",
                  "early_M", "late_M", "net_M", "early_F", "late_F", "net_F",
                  "max_abs_log2FC_M", "max_abs_log2FC_F") |>
    dplyr::mutate(same_shape = .data$shape_M == .data$shape_F)
}

# --- 3. Negative-binomial GLM with empirical FDR ---------------------------

#' Likelihood ratio test of `full` against `reduced` with glmGamPoi.
glm_lrt <- function(counts, meta, full, reduced) {
  fit <- glmGamPoi::glm_gp(counts, design = full, col_data = meta,
                           size_factors = "normed_sum", on_disk = FALSE)
  list(fit = fit, res = glmGamPoi::test_de(fit, reduced_design = reduced))
}

#' p-value null built by refitting the model on shuffled metadata.
permutation_null <- function(counts, meta, full, reduced, shuffle, n_perm) {
  unlist(lapply(seq_len(n_perm), function(i) {
    glm_lrt(counts, shuffle(meta), full, reduced)$res$pval
  }))
}

#' Largest p-value cutoff whose empirical FDR is still at or below target.
calibrate_threshold <- function(obs_p, null_p, target = 0.05) {
  grid <- 10^seq(-10, log10(0.05), length.out = 300)
  fdr <- vapply(grid, function(cut)
    (mean(null_p < cut) * length(obs_p)) / max(1, sum(obs_p < cut)), numeric(1))
  ok <- which(fdr <= target)
  if (!length(ok)) return(list(threshold = 0, fdr = NA_real_))
  i <- max(ok)
  list(threshold = grid[i], fdr = fdr[i])
}

#' Shuffle age within the object (breaks the age effect, keeps depth structure).
shuffle_age <- function(meta) { meta$age <- sample(meta$age); meta }

#' Shuffle sex within each age (keeps the age main effect intact).
shuffle_sex_within_age <- function(meta) {
  shuffled <- stats::ave(as.character(meta$sex), meta$age, FUN = sample)
  meta$sex <- factor(shuffled, levels = levels(meta$sex))
  meta
}

#' Shuffle the (age, sex) pair together, leaving pseudotime fixed.
shuffle_age_sex <- function(meta) {
  i <- sample(nrow(meta))
  meta$age <- meta$age[i]
  meta$sex <- meta$sex[i]
  meta
}

#' Fitted group means on the log scale, from a treatment-coded ~ age fit.
#'
#' Column 1 of Beta is the intercept, i.e. the reference age; each later
#' column is that level's offset from it.
glm_group_means <- function(fit, age_levels) {
  B <- fit$Beta
  if (ncol(B) != length(age_levels))
    stop("fit has ", ncol(B), " coefficients but ", length(age_levels),
         " age levels were given; is the design ~ age with treatment coding?")
  means <- cbind(B[, 1], B[, 1] + B[, -1, drop = FALSE])
  colnames(means) <- age_levels
  means
}

#' Significant genes at the calibrated threshold, with z-scored fitted means.
glm_hits_table <- function(out, calibration, z = NULL) {
  tb <- out$res[out$res$pval < calibration$threshold,
                c("name", "pval", "adj_pval", "f_statistic")]
  tb <- tb[order(tb$pval), ]
  if (!is.null(z)) tb <- cbind(tb, z[match(tb$name, rownames(z)), , drop = FALSE])
  tb
}

# ---------------------------------------------------------------------------
# Depth-aware effect sizes.
#
# avg_log2FC is incoherent under depth asymmetry. Seurat's LogNormalize
# rescales each cell's total, but rescaling cannot restore a transcript that
# was never captured, so detection rates stay biased and a fold change computed
# from them inherits that bias -- which is how a gene can be called "up in
# male" while being detected in fewer male cells.
#
# The replacement puts depth in the model rather than in the normalisation: a
# negative-binomial GLM with library size as an explicit offset. The p-values
# it produces are still per-cell and still answer a question about cells within
# one library, so they are not carried into any reported table; the coefficient
# and its standard error are what this is for.
# ---------------------------------------------------------------------------

#' NB-GLM log fold change with an explicit size-factor offset.
#'
#' @return one row per gene: the coefficient on `contrast`, its standard error
#'   and a Wald confidence interval. No p-value column, deliberately.
glm_effect_size <- function(counts, meta, contrast = "sex", covariates = character(),
                            conf = 0.95) {
  require_packages("glmGamPoi")

  meta <- as.data.frame(meta)
  terms <- c(contrast, covariates)
  design <- stats::as.formula(paste("~", paste(terms, collapse = " + ")))

  keep <- stats::complete.cases(meta[, terms, drop = FALSE])
  if (sum(keep) < ncol(counts))
    log_step(sprintf("  dropping %d cells with missing %s", sum(!keep),
                     paste(terms, collapse = "/")))
  counts <- counts[, keep, drop = FALSE]
  meta <- droplevels(meta[keep, , drop = FALSE])

  fit <- glmGamPoi::glm_gp(counts, design = design, col_data = meta,
                           size_factors = "normed_sum", on_disk = FALSE)

  # The contrast column is the one that is not the intercept and not a covariate.
  coef_names <- colnames(fit$Beta)
  target <- grep(paste0("^", contrast), coef_names, value = TRUE)[1]
  if (is.na(target))
    stop("no coefficient matching '", contrast, "' in: ", paste(coef_names, collapse = ", "))

  estimate <- fit$Beta[, target]
  # glm_gp reports overdispersion, not per-coefficient SEs, so the Wald SE comes
  # from the deviance-based test it can run.
  res <- glmGamPoi::test_de(fit, contrast = target)
  se <- abs(estimate) / stats::qnorm(1 - res$pval / 2)
  se[!is.finite(se)] <- NA_real_
  z <- stats::qnorm(1 - (1 - conf) / 2)

  data.frame(
    gene = rownames(counts),
    log2FC = estimate / log(2),
    ci_lower = (estimate - z * se) / log(2),
    ci_upper = (estimate + z * se) / log(2),
    # Kept for diagnostics only, never for a reported table.
    .pval_diagnostic = res$pval,
    row.names = NULL, stringsAsFactors = FALSE
  )
}

#' Detection rate per group, and the direction it points.
#'
#' The cheapest check on whether a fold change is real: a gene called higher in
#' one group while being detected in fewer of that group's cells is describing
#' depth, not biology.
detection_rates <- function(counts, groups) {
  detected <- lapply(split(seq_along(groups), groups), function(idx)
    Matrix::rowMeans(counts[, idx, drop = FALSE] > 0))
  out <- as.data.frame(detected)
  names(out) <- paste0("pct_", names(detected))
  out$gene <- rownames(counts)
  out
}

#' Flag effect sizes that disagree with their own detection rates.
flag_detection_conflict <- function(effects, detection, higher_in) {
  cols <- grep("^pct_", names(detection), value = TRUE)
  if (length(cols) != 2) return(effects)

  merged <- merge(effects, detection, by = "gene", sort = FALSE)
  higher_col <- paste0("pct_", higher_in)
  other_col <- setdiff(cols, higher_col)

  merged$detected_more_in_claimed_group <-
    merged[[higher_col]] >= merged[[other_col]]
  conflicts <- sum(merged$log2FC > 0 & !merged$detected_more_in_claimed_group,
                   na.rm = TRUE)
  up <- sum(merged$log2FC > 0, na.rm = TRUE)
  if (up > 0)
    log_step(sprintf("  %d of %d genes called up in %s are detected in FEWER %s cells (%.0f%%)",
                     conflicts, up, higher_in, higher_in, 100 * conflicts / up))
  merged
}

# ---------------------------------------------------------------------------
# Age trends, calibrated against a permutation null.
#
# A bare |rho| cutoff has no null and is n-dependent, so the count of "changing"
# genes it returns is largely a statement about how many cells the stratum had.
# Permuting age within the stratum gives the rho distribution this stratum
# produces by chance, and what gets reported is the excess over it.
# ---------------------------------------------------------------------------

#' Null distribution of |rho| from permuted age labels.
age_trend_null <- function(obj, age_levels, n_perm = 100, age_col = "age",
                           assay = "RNA", seed = 42) {
  expr <- Seurat::GetAssayData(obj, assay = assay, layer = "data")
  ages <- as.character(obj[[age_col]][, 1])

  vapply(seq_len(n_perm), function(i) {
    set.seed(seed + i)
    shuffled <- sample(ages)
    numeric_age <- as.numeric(factor(shuffled, levels = age_levels))
    max(abs(spearman_rows(expr, numeric_age)$rho), na.rm = TRUE)
  }, numeric(1))
}

#' Genes whose age trend exceeds what permutation produces for this stratum.
#'
#' Reports excess over null rather than a raw count, so strata of different
#' sizes can be compared. A stratum below the gate returns NULL rather than a
#' number, since a count from 6 cells is not a smaller version of a count from
#' 600 -- it is not a measurement.
age_trend_excess <- function(obj, cfg, age_levels = cfg$analysis$age_levels,
                             age_col = "age", assay = "RNA", n_perm = 100,
                             label = "") {
  n_by_age <- table(as.character(obj[[age_col]][, 1]))
  if (!stratum_is_usable(n_by_age, label, cfg)) return(NULL)

  observed <- age_trend_test(obj, age_levels, age_col, assay)
  null_max <- age_trend_null(obj, age_levels, n_perm, age_col, assay, cfg$seed)
  threshold <- stats::quantile(null_max, 0.95, names = FALSE, na.rm = TRUE)

  n_observed <- sum(abs(observed$rho) > threshold, na.rm = TRUE)
  log_step(sprintf("  %s: |rho| threshold from permutation = %.3f; %d genes exceed it",
                   label, threshold, n_observed))

  # What the permutation null CANNOT protect against.
  #
  # Shuffling age labels destroys any association between age and library
  # complexity, so the null describes a stratum in which complexity does not
  # trend with age. If it does trend in the real data -- each age is a
  # separate hashtag, and hashtags differ in how many genes they detect --
  # then every gene whose detection follows complexity carries that trend into
  # the observed rho, clears a threshold built without it, and is counted as
  # signal. A large gene count and a strong detection trend are the same
  # observation reported twice.
  detected <- Matrix::colSums(
    Seurat::GetAssayData(obj, assay = assay, layer = "counts") > 0)
  numeric_age <- as.numeric(factor(as.character(obj[[age_col]][, 1]),
                                   levels = age_levels))
  detection_rho <- suppressWarnings(
    stats::cor(detected, numeric_age, method = "spearman", use = "complete.obs"))

  if (is.finite(detection_rho) && abs(detection_rho) > threshold)
    log_step(sprintf(
      paste0("  %s: WARNING genes detected trends with age at rho = %.3f, above ",
             "this stratum's own threshold (%.3f). The %d genes above may be ",
             "that one trend rather than %d independent findings, and ",
             "permuting age cannot calibrate against it."),
      label, detection_rho, threshold, n_observed, n_observed))
  else
    log_step(sprintf(
      paste0("  %s: genes detected trends with age at rho = %.3f, below the ",
             "%.3f threshold, so the gene count is not a complexity trend."),
      label, detection_rho, threshold))

  list(trend = observed,
       null_threshold = threshold,
       null_max = null_max,
       n_exceeding = n_observed,
       detection_rho = detection_rho,
       detection_confounded = is.finite(detection_rho) &&
         abs(detection_rho) > threshold,
       changing_genes = observed[which(abs(observed$rho) > threshold), ],
       n_cells = as.integer(n_by_age))
}

# ---------------------------------------------------------------------------
# Set-level scoring.
#
# Coherence across a module survives depth noise that no individual gene
# survives, which is why a pathway can be the strongest result in a dataset
# where one or two of its genes reach significance. Module scores lead; gene
# hits illustrate.
# ---------------------------------------------------------------------------

#' Rank-based module scores (UCell), which depend on gene ordering within a
#' cell rather than on absolute counts, and so are far less depth-sensitive
#' than a mean-expression score.
add_module_scores_ucell <- function(obj, signatures, assay = "RNA", suffix = "_UCell") {
  require_packages("UCell")
  UCell::AddModuleScore_UCell(obj, features = signatures, assay = assay, name = suffix)
}

#' Compare a module score between two groups, as an effect size with a CI.
module_effect <- function(obj, score_col, group_col, groups, n_boot = 1000) {
  values <- obj@meta.data[[score_col]]
  labels <- as.character(obj[[group_col]][, 1])

  a <- values[labels == groups[1]]
  b <- values[labels == groups[2]]
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 3 || length(b) < 3)
    return(data.frame(module = score_col, effect_size = NA_real_,
                      ci_lower = NA_real_, ci_upper = NA_real_,
                      direction = NA_character_, n_min = min(length(a), length(b))))

  # Standardised mean difference, so modules with different score ranges are
  # comparable to each other.
  pooled_sd <- sqrt(((length(a) - 1) * stats::var(a) + (length(b) - 1) * stats::var(b)) /
                      (length(a) + length(b) - 2))
  effect <- function(idx_a, idx_b) (mean(idx_a) - mean(idx_b)) / pooled_sd

  boots <- vapply(seq_len(n_boot), function(i)
    effect(sample(a, length(a), TRUE), sample(b, length(b), TRUE)), numeric(1))

  data.frame(module = score_col,
             effect_size = effect(a, b),
             ci_lower = stats::quantile(boots, 0.025, names = FALSE),
             ci_upper = stats::quantile(boots, 0.975, names = FALSE),
             direction = if (mean(a) > mean(b)) groups[1] else groups[2],
             n_min = min(length(a), length(b)),
             stringsAsFactors = FALSE)
}
