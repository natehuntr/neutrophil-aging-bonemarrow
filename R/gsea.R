# ---------------------------------------------------------------------------
# Sex x age interaction GSEA.
#
# Ranking metric: Fisher-z contrast of per-cell Spearman rho (expression vs
# numeric age), computed independently within each sex on one shared gene
# universe.
#
#   z_diff = (atanh(rho_M) - atanh(rho_F)) / sqrt(1/(n_M - 3) + 1/(n_F - 3))
#
# z_diff > 0 means ONLY that the age trajectory is more positive in males than
# in females. It does NOT mean the gene is higher in males, nor that males age
# "faster". Direction is resolved per pathway by describe_leading_edge(),
# which goes back to the two underlying rho values for every leading-edge gene.
#
# Caveat kept from the original analysis: cells are pseudoreplicates of pooled
# animals, so z_diff is not per-gene evidence. Its job is to produce a sensibly
# scaled ranking, which is all GSEA needs.
# ---------------------------------------------------------------------------

# Genes flagged (never auto-removed) in leading edges.
SEX_CHR_GENES <- c("Xist", "Tsix", "Ddx3y", "Uty", "Eif2s3y", "Kdm5d", "Jpx", "Ftx")

# Highly abundant granule transcripts, for the compositional sensitivity run.
GRANULE_GENES <- c("Ngp", "Lcn2", "Ltf", "Camp", "Chil3", "Serpinb1a",
                   "Wfdc21", "Mmp8", "Retnlg", "S100a8", "S100a9")

#' Spearman rho of every detected gene against numeric age, within one sex.
#'
#' Spearman is computed as Pearson on ranks, which is far faster than
#' cor.test() per gene and gives the same rho.
compute_rho_by_sex <- function(obj, sex_value, cfg,
                               sex_col = "sex", age_col = "age", assay = "RNA") {
  cells <- colnames(obj)[as.character(obj[[sex_col]][, 1]) == sex_value]
  if (!length(cells)) stop("no cells matched ", sex_col, " == '", sex_value, "'")
  sub <- subset(obj, cells = cells)

  expr <- Seurat::GetAssayData(sub, assay = assay, layer = "data")
  age <- suppressWarnings(
    as.numeric(gsub("[^0-9.]", "", as.character(sub[[age_col]][, 1]))))
  if (anyNA(age)) stop("could not parse ages for ", sex_value)
  stopifnot(length(age) == ncol(expr))

  pct <- Matrix::rowMeans(expr > 0)
  expr <- expr[names(pct)[pct >= cfg$gsea$min_pct_cells], , drop = FALSE]

  age_ranks <- rank(age)
  expr_ranks <- t(apply(as.matrix(expr), 1, rank))
  rho <- as.numeric(stats::cor(t(expr_ranks), age_ranks))

  tibble::tibble(gene = rownames(expr), rho = rho, n = ncol(expr))
}

#' Fisher-z contrast of the two per-sex rho tables, with a per-gene label
#' describing what each sex's trend actually does.
fisher_z_contrast <- function(rho_male, rho_female, rho_flat = 0.10, clamp = 0.999) {
  trend <- function(rho) dplyr::case_when(rho >  rho_flat ~ "up",
                                          rho < -rho_flat ~ "down",
                                          TRUE            ~ "flat")

  dplyr::inner_join(
      dplyr::select(rho_male,   gene = "gene", rho_M = "rho", n_M = "n"),
      dplyr::select(rho_female, gene = "gene", rho_F = "rho", n_F = "n"),
      by = "gene") |>
    dplyr::mutate(
      z_M = atanh(pmax(pmin(.data$rho_M, clamp), -clamp)),
      z_F = atanh(pmax(pmin(.data$rho_F, clamp), -clamp)),
      se  = sqrt(1 / (.data$n_M - 3) + 1 / (.data$n_F - 3)),
      z_diff = (.data$z_M - .data$z_F) / .data$se,
      trend_M = trend(.data$rho_M),
      trend_F = trend(.data$rho_F),
      pattern = dplyr::case_when(
        trend_M == "up"   & trend_F == "down" ~ "opposing: up in M, down in F",
        trend_M == "down" & trend_F == "up"   ~ "opposing: down in M, up in F",
        trend_M == "flat" & trend_F == "down" ~ "female-only: declines in F",
        trend_M == "flat" & trend_F == "up"   ~ "female-only: rises in F",
        trend_M == "down" & trend_F == "flat" ~ "male-only: declines in M",
        trend_M == "up"   & trend_F == "flat" ~ "male-only: rises in M",
        trend_M == "up"   & trend_F == "up"   ~ "shared: rises in both",
        trend_M == "down" & trend_F == "down" ~ "shared: declines in both",
        TRUE                                  ~ "no trend in either"),
      divergent = grepl("^opposing", .data$pattern)) |>
    dplyr::filter(is.finite(.data$z_diff)) |>
    dplyr::arrange(dplyr::desc(.data$z_diff))
}

#' GO:BP gene sets restricted to the tested universe, so set sizes are honest.
gobp_pathways <- function(universe, cfg) {
  gs <- msigdbr::msigdbr(species = "Mus musculus",
                         collection = cfg$gsea$collection,
                         subcollection = cfg$gsea$subcollection)
  pathways <- split(gs$gene_symbol, gs$gs_name)
  pathways <- lapply(pathways, function(g) intersect(unique(g), universe))
  pathways <- pathways[lengths(pathways) >= cfg$gsea$min_set_size]
  log_step("gene sets tested: ", length(pathways))
  pathways
}

#' What the leading edge of one pathway is actually made of.
describe_leading_edge <- function(leading_edge, lookup) {
  tab <- lookup[lookup$gene %in% leading_edge, ]
  if (!nrow(tab))
    return(tibble::tibble(dominant_pattern = NA_character_, dominant_n = NA_integer_,
                          dominant_frac = NA_real_, divergent_frac = NA_real_,
                          mean_rho_M = NA_real_, mean_rho_F = NA_real_))
  counts <- sort(table(tab$pattern), decreasing = TRUE)
  tibble::tibble(
    dominant_pattern = names(counts)[1],
    dominant_n       = as.integer(counts[1]),
    dominant_frac    = as.numeric(counts[1]) / nrow(tab),
    divergent_frac   = mean(tab$divergent),
    mean_rho_M       = mean(tab$rho_M),
    mean_rho_F       = mean(tab$rho_F)
  )
}

#' Preranked GSEA on the interaction statistic, annotated so the sign of NES
#' can be read correctly.
run_interaction_gsea <- function(rank_tbl, pathways, cfg) {
  ranks <- stats::setNames(rank_tbl$z_diff, rank_tbl$gene)

  res <- fgsea::fgseaMultilevel(pathways, ranks,
                                minSize = cfg$gsea$min_set_size,
                                maxSize = cfg$gsea$max_set_size, eps = 0)

  ordered <- res[order(res$pval), ]
  sig <- ordered[ordered$padj < 0.05, ]
  # GO:BP is heavily redundant; collapse nested and overlapping sets.
  main <- if (nrow(sig) > 0) fgsea::collapsePathways(sig, pathways, ranks)
          else list(mainPathways = character(0))

  lookup <- rank_tbl[, c("gene", "rho_M", "rho_F", "pattern", "divergent")]
  le_summary <- purrr::map_dfr(res$leadingEdge, describe_leading_edge, lookup = lookup)

  out <- dplyr::bind_cols(tibble::as_tibble(res), le_summary) |>
    dplyr::mutate(
      independent = .data$pathway %in% main$mainPathways,
      # Neutral restatement of the NES sign: no claim about which sex is worse.
      trajectory_more_positive_in = ifelse(.data$NES > 0, "male", "female"),
      signal_type = dplyr::case_when(
        .data$divergent_frac >= 0.5                    ~ "opposing trajectories",
        grepl("^male-only",   .data$dominant_pattern)  ~ "male-specific change",
        grepl("^female-only", .data$dominant_pattern)  ~ "female-specific change",
        grepl("^shared",      .data$dominant_pattern)  ~ "shared direction, differing magnitude",
        TRUE                                           ~ "mixed / weak"),
      leadingEdge_n = lengths(.data$leadingEdge),
      interpretation = sprintf(
        "trajectory more positive in %s; leading edge mostly '%s' (%d/%d genes); mean rho M=%.2f F=%.2f",
        .data$trajectory_more_positive_in, .data$dominant_pattern,
        .data$dominant_n, .data$leadingEdge_n, .data$mean_rho_M, .data$mean_rho_F),
      sex_chr_in_LE = vapply(.data$leadingEdge,
                             function(g) paste(intersect(g, SEX_CHR_GENES), collapse = ";"),
                             character(1)),
      granule_in_LE = vapply(.data$leadingEdge,
                             function(g) paste(intersect(g, GRANULE_GENES), collapse = ";"),
                             character(1)),
      leadingEdge = vapply(.data$leadingEdge, paste, character(1), collapse = ";")) |>
    dplyr::select("pathway", "NES", "pval", "padj", "independent", "signal_type",
                  "trajectory_more_positive_in", "dominant_pattern", "dominant_n",
                  "dominant_frac", "divergent_frac", "mean_rho_M", "mean_rho_F",
                  "leadingEdge_n", "sex_chr_in_LE", "granule_in_LE",
                  "interpretation", "leadingEdge") |>
    dplyr::arrange(.data$padj, dplyr::desc(abs(.data$NES)))

  list(result = out, ranks = ranks, main_pathways = main$mainPathways)
}

#' Sensitivity run with the abundant granule transcripts removed.
#'
#' Guards against the 18m compositional-dilution signature driving secretory
#' and translation sets. Pathways surviving both runs are the ones to report.
granule_sensitivity <- function(gsea_out, pathways, cfg) {
  ranks_sens <- gsea_out$ranks[!names(gsea_out$ranks) %in% GRANULE_GENES]
  res_sens <- fgsea::fgseaMultilevel(lapply(pathways, setdiff, GRANULE_GENES),
                                     ranks_sens,
                                     minSize = cfg$gsea$min_set_size,
                                     maxSize = cfg$gsea$max_set_size, eps = 0)

  gsea_out$result |>
    dplyr::select("pathway", NES_full = "NES", padj_full = "padj", "independent") |>
    dplyr::left_join(tibble::as_tibble(res_sens) |>
                       dplyr::select("pathway", NES_sens = "NES", padj_sens = "padj"),
                     by = "pathway") |>
    dplyr::mutate(robust = .data$padj_full < 0.05 & .data$padj_sens < 0.05 &
                    sign(.data$NES_full) == sign(.data$NES_sens)) |>
    dplyr::arrange(dplyr::desc(.data$robust), .data$padj_full)
}

#' Male-vs-female Wilcoxon within one age.
#'
#' Ranked on avg_log2FC, not on the p-value: cell numbers differ ~10-fold
#' between the arms at 12m, and a Wilcoxon p is a function of that imbalance
#' whereas log2FC is not.
sex_de_at_age <- function(obj, this_age, cfg,
                          sex_col = "sex", age_col = "age", assay = "RNA") {
  cells <- colnames(obj)[as.character(obj[[age_col]][, 1]) == this_age]
  sub <- subset(obj, cells = cells)
  Seurat::Idents(sub) <- sub[[sex_col]][, 1]
  log_step("age ", this_age, ": ", paste(table(Seurat::Idents(sub)), collapse = " / "), " cells")

  markers <- Seurat::FindMarkers(sub, ident.1 = "male", ident.2 = "female",
                                 assay = assay, layer = "data",
                                 logfc.threshold = 0, min.pct = cfg$gsea$min_pct_cells,
                                 test.use = "wilcox", verbose = FALSE)
  tibble::rownames_to_column(markers, "gene")
}

#' Per-age sex contrasts lined up against the interaction statistic.
#'
#' All four NES columns are computed on one shared universe; without that they
#' sit on different backgrounds and are not comparable, which is the whole
#' point of the comparison.
compare_perage_to_interaction <- function(obj, rank_tbl, pathways, cfg,
                                          age_levels = cfg$analysis$age_levels) {
  ranks <- stats::setNames(rank_tbl$z_diff, rank_tbl$gene)
  de_list <- lapply(stats::setNames(age_levels, age_levels),
                    function(a) sex_de_at_age(obj, a, cfg))

  universe <- Reduce(intersect, c(lapply(de_list, `[[`, "gene"), list(names(ranks))))
  log_step("shared universe: ", length(universe), " genes")

  rank_at_age <- lapply(de_list, function(d) {
    d <- d[d$gene %in% universe, ]
    d <- d[order(-d$avg_log2FC), ]
    stats::setNames(d$avg_log2FC, d$gene)
  })

  pathways_u <- lapply(pathways, intersect, universe)
  pathways_u <- pathways_u[lengths(pathways_u) >= cfg$gsea$min_set_size]

  gsea_by_age <- lapply(names(rank_at_age), function(nm) {
    fgsea::fgseaMultilevel(pathways_u, rank_at_age[[nm]],
                           minSize = cfg$gsea$min_set_size,
                           maxSize = cfg$gsea$max_set_size, eps = 0) |>
      tibble::as_tibble() |>
      dplyr::select("pathway", "NES", "pval", "padj") |>
      dplyr::rename_with(~ paste0(.x, "_", nm), c("NES", "pval", "padj"))
  })

  res_int <- fgsea::fgseaMultilevel(pathways_u, ranks[universe],
                                    minSize = cfg$gsea$min_set_size,
                                    maxSize = cfg$gsea$max_set_size, eps = 0) |>
    tibble::as_tibble() |>
    dplyr::select("pathway", NES_int = "NES", padj_int = "padj")

  combined <- Reduce(function(a, b) dplyr::full_join(a, b, by = "pathway"), gsea_by_age)
  combined <- dplyr::full_join(combined, res_int, by = "pathway")

  # These summaries are built in base R: the column names are driven by
  # age_levels, so tidy evaluation would only obscure what is happening.
  first <- age_levels[1]
  last  <- age_levels[length(age_levels)]
  nes_cols  <- paste0("NES_",  age_levels)
  padj_cols <- paste0("padj_", age_levels)

  nes_mat  <- as.matrix(combined[, nes_cols,  drop = FALSE])
  padj_mat <- as.matrix(combined[, padj_cols, drop = FALSE])

  # how much the sex difference moves across the lifespan
  combined$delta_last_first <- combined[[paste0("NES_", last)]] -
                               combined[[paste0("NES_", first)]]
  # is the direction of the sex difference stable across all ages?
  combined$consistent_sign <- apply(sign(nes_mat), 1,
                                    function(x) length(unique(x)) == 1)
  # does the per-age picture agree with the interaction statistic?
  combined$agrees_with_int <- sign(combined$delta_last_first) == sign(combined$NES_int)
  combined$sig_any_age <- apply(padj_mat, 1, function(x) min(x, na.rm = TRUE)) < 0.05
  combined$sig_int <- combined$padj_int < 0.05

  combined$class <- with(combined, dplyr::case_when(
    sig_int & sig_any_age & agrees_with_int ~ "sex-differential ageing (supported)",
    sig_int & !agrees_with_int              ~ "interaction only - check middle-age leverage",
    !sig_int & sig_any_age & consistent_sign &
      abs(delta_last_first) < 0.5           ~ "constitutive sex difference",
    sig_any_age                             ~ "age-specific sex difference",
    TRUE                                    ~ "ns"))

  combined <- combined[order(combined$padj_int), ]
  combined
}

#' Jaccard overlap between the significant sets from each age and the
#' interaction test.
significant_set_overlap <- function(combined, age_levels) {
  sets <- lapply(stats::setNames(age_levels, age_levels), function(a)
    combined$pathway[which(combined[[paste0("padj_", a)]] < 0.05)])
  sets$interaction <- combined$pathway[which(combined$padj_int < 0.05)]

  jaccard <- outer(seq_along(sets), seq_along(sets), Vectorize(function(i, j) {
    a <- sets[[i]]; b <- sets[[j]]
    if (!length(union(a, b))) return(NA_real_)
    length(intersect(a, b)) / length(union(a, b))
  }))
  dimnames(jaccard) <- list(names(sets), names(sets))
  list(sizes = lengths(sets), jaccard = round(jaccard, 2))
}
