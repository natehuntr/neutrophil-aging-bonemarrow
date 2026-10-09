# ---------------------------------------------------------------------------
# Summary figures, drawn from the results TABLES alone.
#
# Every function here takes data frames read from results/tables and returns a
# ggplot. Nothing reads a Seurat object, so step 13 runs in seconds, on a
# laptop, on a results bundle -- and a figure can be redrawn without re-running
# the analysis that produced its numbers.
#
# One figure per claim the results support or refuse, numbered in the order a
# reader needs them: what the data can be trusted for first (depth, staging),
# then composition, position, expression, protein, and the sex contrast last.
# ---------------------------------------------------------------------------

AGE_ORDER <- c("3m", "9m", "12m", "18m")
STAGE_ORDER <- c("GMPs", "proNeu", "preNeu", "immature", "mature")

#' Factor ages and stages in their biological order.
order_levels <- function(df) {
  if ("age" %in% names(df)) df$age <- factor(as.character(df$age),
                                              levels = intersect(AGE_ORDER, unique(as.character(df$age))))
  if ("stage" %in% names(df)) df$stage <- factor(as.character(df$stage),
                                                  levels = unique(c(intersect(STAGE_ORDER, df$stage), df$stage)))
  df
}

#' Short readable progenitor names.
short_progenitor <- function(x) {
  map <- c("Long-Term HSC (LT-HSC)" = "LT-HSC",
           "Long-Term HSC (CD34- Flt3-)" = "LT-HSC (CD34- Flt3-)",
           "Short-Term HSC (CD34+ Flt3-)" = "ST-HSC (CD34+ Flt3-)",
           "Short-Term HSC (Sca1+ Lin-, SLAM)" = "ST-HSC (SLAM)",
           "HSC subset (CD150- CD48-, SLAM MPP)" = "SLAM MPP",
           "Multipotent Progenitor (CD34+ Flt3-)" = "MPP (CD34+ Flt3-)",
           "Multipotent Lymphoid Progenitor" = "MLP",
           "Common Myeloid Progenitor (CMP)" = "CMP",
           "Granulocyte-Monocyte Progenitor (GMP)" = "GMP",
           "Megakaryocyte-Erythroid Progenitor (MEP)" = "MEP",
           "Common Lymphoid Progenitor/pro-B cell" = "CLP / pro-B",
           "Macrophage-DC Progenitor (MDP)" = "MDP",
           "Common Dendritic Cell Progenitor (CDP)" = "CDP")
  out <- unname(map[as.character(x)])
  ifelse(is.na(out), as.character(x), out)
}

caption_cells <- "Intervals are on cells: they show sampling of captured cells, not variation between mice."

# --- 1. Depth --------------------------------------------------------------

#' Median UMIs per cell by stage, age and sex, with the gate's tolerance.
#'
#' The figure that explains which comparisons are depth-controlled. Lines that
#' slope within a library are depth varying between hashtags; a stage sitting
#' far below the others is a stage global matching could not reach.
fig_depth_by_stage <- function(depth) {
  d <- depth[depth$stage != "ALL STAGES (composition-confounded)", ]
  d <- order_levels(d)
  ggplot2::ggplot(d, ggplot2::aes(.data$age, .data$median_umi, colour = .data$sex,
                                  group = .data$sex)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n_cells)) +
    ggplot2::scale_size_area(max_size = 4, name = "Cells") +
    ggplot2::facet_wrap(~ stage, nrow = 1) +
    ggplot2::scale_y_log10(labels = scales::comma,
                           breaks = c(250, 500, 1000, 2000, 4000, 8000)) +
    scale_sex() +
    ggplot2::labs(x = NULL, y = "Median UMIs per cell (log scale)",
                  title = "Sequencing depth by stage, age and library",
                  subtitle = paste("Counts the analyses read", if ("assay" %in% names(d))
                    paste0("(", unique(d$assay)[1], ")") else "",
                    "-- a slope within one library is depth changing between hashtags")) +
    theme_bm()
}

#' Depth ratio inside each matching stratum, before and after matching.
fig_matching_targets <- function(targets) {
  t <- targets
  t$stratum <- factor(t$stratum, levels = rev(unique(t$stratum)))
  long <- rbind(data.frame(t[, c("assay", "stratum")], when = "before", ratio = t$ratio_before),
                data.frame(t[, c("assay", "stratum")], when = "after", ratio = t$ratio_after))
  long$when <- factor(long$when, levels = c("before", "after"))
  ggplot2::ggplot(long, ggplot2::aes(.data$ratio, .data$stratum)) +
    ggplot2::geom_vline(xintercept = 1.3, colour = MUTED, linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(group = .data$stratum), colour = GRID, linewidth = 1.2) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$when), size = 2.4) +
    ggplot2::scale_colour_manual(values = c(before = "#86b6ef", after = "#0d366b"), name = NULL) +
    ggplot2::facet_wrap(~ assay, scales = "free_y") +
    ggplot2::labs(x = "Deepest / shallowest group median UMIs", y = NULL,
                  title = "Depth matching, stratum by stratum",
                  subtitle = "Dashed: the 1.3x gate. Each stratum is thinned toward its own shallowest group.") +
    theme_bm()
}

# --- 2. Staging -------------------------------------------------------------

#' RNA stage against protein (ADT) stage, as row fractions.
fig_stage_confusion <- function(confusion) {
  names(confusion)[1:3] <- c("rna", "adt", "n")
  confusion$adt[is.na(confusion$adt) | confusion$adt == ""] <- "unassigned"
  totals <- stats::aggregate(n ~ rna, confusion, sum)
  confusion$frac <- confusion$n / totals$n[match(confusion$rna, totals$rna)]
  confusion$rna <- factor(confusion$rna, levels = rev(intersect(STAGE_ORDER, confusion$rna)))
  confusion$adt <- factor(confusion$adt, levels = c(intersect(STAGE_ORDER, confusion$adt), "unassigned"))
  ggplot2::ggplot(confusion, ggplot2::aes(.data$adt, .data$rna, fill = .data$frac)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = .data$n,
                                    colour = .data$frac > 0.5), size = 3.2, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK)) +
    ggplot2::scale_fill_gradient(low = "#f4f7fb", high = "#104281", name = "Share of\nRNA stage",
                                 labels = scales::percent) +
    ggplot2::labs(x = "Stage from surface protein (ADT panel)", y = "Stage from RNA signatures",
                  title = "Do RNA and protein agree on maturation stage?",
                  subtitle = "Numbers are cells; shading is the share of each RNA stage. The diagonal is agreement.") +
    theme_bm() + ggplot2::theme(panel.grid = ggplot2::element_blank())
}

# --- 3. Composition ---------------------------------------------------------

#' Stage mix per age, one panel per library.
fig_stage_composition <- function(comp) {
  comp <- order_levels(comp)
  labels <- comp[comp$proportion >= 0.08, ]
  ggplot2::ggplot(comp, ggplot2::aes(.data$age, .data$proportion, fill = .data$stage)) +
    ggplot2::geom_col(width = 0.68, colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(data = labels, ggplot2::aes(label = scales::percent(.data$proportion, 1),
                                                   colour = .data$stage %in% c("immature", "mature", "preNeu")),
                       position = ggplot2::position_stack(vjust = 0.5), size = 3, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK)) +
    ggplot2::geom_text(data = unique(comp[, c("sex", "age", "n_group")]),
                       ggplot2::aes(.data$age, -0.04, label = paste0("n=", .data$n_group)),
                       inherit.aes = FALSE, size = 2.8, colour = MUTED) +
    ggplot2::facet_wrap(~ sex) +
    scale_stage() +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = NULL, y = "Share of staged granulocytic cells",
                  title = "Maturation stage mix by age",
                  subtitle = "Ages share a library, so differences within a panel are not batch. Panels are different libraries.") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    theme_bm()
}

#' Each stage's share against age, with Wilson intervals and the fitted slope.
fig_stage_trends <- function(comp, trends = NULL) {
  comp <- order_levels(comp)
  p <- ggplot2::ggplot(comp, ggplot2::aes(.data$age, .data$proportion, colour = .data$sex,
                                          group = .data$sex)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high), size = 0.3) +
    ggplot2::facet_wrap(~ stage, scales = "free_y", nrow = 1) +
    scale_sex() +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = NULL, y = "Share of staged cells",
                  title = "Stage shares across age, per library", caption = caption_cells) +
    theme_bm()
  if (!is.null(trends) && nrow(trends)) {
    lab <- order_levels(trends)
    lab$text <- sprintf("%s %+.2f", substr(lab$sex, 1, 1), lab$log_odds_per_step)
    lab <- stats::aggregate(text ~ stage, lab, paste, collapse = "   ")
    p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(x = -Inf, y = Inf, label = .data$text),
                                inherit.aes = FALSE, hjust = -0.05, vjust = 1.4, size = 2.7,
                                colour = MUTED) +
      ggplot2::labs(subtitle = "Top of each panel: log-odds change per age step, f = female, m = male")
  }
  p
}

#' Over- and under-representation of a population per age, within each sex.
#'
#' log2(observed / expected) from the sex's own age-pooled mix, with cell-level
#' bootstrap intervals. Used for maturation stages and for progenitors.
fig_imbalance <- function(imb, title, population_label = identity, order = NULL) {
  d <- imb[imb$contrast == "across age within sex", ]
  if (!nrow(d)) return(NULL)
  d$population <- population_label(d$population)
  d$age <- factor(d$across_level, levels = intersect(AGE_ORDER, d$across_level))
  d$sex <- d$within_level
  ord <- stats::aggregate(log2_obs_exp ~ population, d[d$age == utils::tail(levels(d$age), 1), ],
                          mean)
  d$population <- factor(d$population, levels = if (!is.null(order))
    rev(intersect(order, d$population)) else ord$population[order(ord$log2_obs_exp)])
  d$sig <- ifelse(d$excludes_zero, "interval excludes 0", "interval spans 0")
  ggplot2::ggplot(d, ggplot2::aes(.data$log2_obs_exp, .data$population, colour = .data$sex)) +
    ggplot2::geom_vline(xintercept = 0, colour = MUTED) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = .data$ci_low, xmax = .data$ci_high),
                            position = ggplot2::position_dodge(width = 0.6), linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$sig), size = 2.2,
                        position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_shape_manual(values = c("interval excludes 0" = 16, "interval spans 0" = 1),
                                name = NULL) +
    ggplot2::facet_wrap(~ age, nrow = 1) +
    scale_sex() +
    ggplot2::coord_cartesian(xlim = c(-3, 3)) +
    ggplot2::labs(x = "log2(observed / expected share)", y = NULL, title = title,
                  subtitle = "Expected = the population's share across all ages of the same library. Right of 0: more than expected.",
                  caption = caption_cells) +
    theme_bm()
}

#' Progenitor shares across age for the populations that move.
fig_progenitor_trends <- function(comp, min_cells = 20) {
  comp <- order_levels(comp)
  comp$population <- short_progenitor(comp$population)
  big <- stats::aggregate(n ~ population, comp, max)
  comp <- comp[comp$population %in% big$population[big$n >= min_cells], ]
  ggplot2::ggplot(comp, ggplot2::aes(.data$age, .data$proportion, colour = .data$sex,
                                     group = .data$sex)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high), size = 0.25) +
    ggplot2::facet_wrap(~ population, scales = "free_y") +
    scale_sex() +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    ggplot2::labs(x = NULL, y = "Share of the progenitor compartment",
                  title = "Progenitor compartment across age",
                  subtitle = paste0("Populations reaching ", min_cells, " cells in some group; SingleR/ImmGen labels"),
                  caption = caption_cells) +
    theme_bm()
}

# --- 4. Position along the trajectory --------------------------------------

#' Median shift from 3m, pooled over stages and within each stage.
#'
#' Read the two rows together: a pooled shift with no within-stage shift is
#' the stage mix moving; a within-stage shift is cells themselves moving.
fig_distribution_shifts <- function(pooled, within, what = "Pseudotime") {
  pooled$stage <- "all stages pooled"
  cols <- c("group", "sex", "stage", "median_shift", "frac_beyond_reference_median", "n")
  d <- rbind(pooled[, cols], if (!is.null(within)) within[, cols])
  d$stage <- factor(d$stage, levels = c("all stages pooled", intersect(STAGE_ORDER, d$stage)))
  d$age <- factor(d$group, levels = intersect(AGE_ORDER, d$group))
  ggplot2::ggplot(d, ggplot2::aes(.data$frac_beyond_reference_median, .data$age,
                                  colour = .data$sex)) +
    ggplot2::geom_vline(xintercept = 0.5, colour = MUTED) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n),
                        position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::scale_size_area(max_size = 4, name = "Cells") +
    ggplot2::facet_wrap(~ stage, nrow = 1) +
    scale_sex() +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    ggplot2::labs(x = paste0("Share of cells past the 3m median ", tolower(what)), y = NULL,
                  title = paste(what, "relative to 3m, pooled and within stage"),
                  subtitle = "50% = no shift. Pooled shifts with flat within-stage panels are the stage mix, not the cells.") +
    theme_bm()
}

# --- 5. Expression ----------------------------------------------------------

#' Endpoint (3m vs 18m) gene counts per stratum, with the depth check.
fig_endpoint_summary <- function(summary, progenitor = NULL) {
  s <- summary
  s$stratum <- paste(s$stage, s$sex, sep = " / ")
  s$kind <- "maturation stage"
  if (!is.null(progenitor) && nrow(progenitor)) {
    p <- progenitor
    p$stratum <- paste(short_progenitor(p$population), p$sex, sep = " / ")
    p$kind <- "progenitor"
    p$n_pathways <- NA
    s <- rbind(s[, c("stratum", "kind", "n_genes", "n_monotonic", "detection_skewed")],
               p[, c("stratum", "kind", "n_genes", "n_monotonic", "detection_skewed")])
  }
  s$stratum <- factor(s$stratum, levels = s$stratum[order(s$n_genes)])
  s$flag <- ifelse(s$detection_skewed %in% TRUE, "detection differs between endpoints",
                   "endpoints matched")
  long <- rbind(data.frame(stratum = s$stratum, kind = s$kind, flag = s$flag,
                           what = "all genes past the threshold", n = s$n_genes),
                data.frame(stratum = s$stratum, kind = s$kind, flag = s$flag,
                           what = "of which monotonic across 4 ages", n = s$n_monotonic))
  ggplot2::ggplot(long, ggplot2::aes(.data$n, .data$stratum, fill = .data$what)) +
    ggplot2::geom_col(position = ggplot2::position_identity(), width = 0.7) +
    ggplot2::geom_point(data = long[long$flag != "endpoints matched" & long$what != "of which monotonic across 4 ages", ],
                        ggplot2::aes(x = -max(long$n) * 0.04), shape = 4, size = 2.5,
                        colour = "#e34948", inherit.aes = TRUE, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c("all genes past the threshold" = "#86b6ef",
                                          "of which monotonic across 4 ages" = "#104281"),
                               name = NULL) +
    ggplot2::facet_grid(kind ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(x = "Genes changing 3m -> 18m (permutation FDR)", y = NULL,
                  title = "How much changes between 3m and 18m, per stratum",
                  subtitle = "Red x: endpoints differ in genes detected per cell, so the list may be that difference") +
    theme_bm()
}

#' Mean expression across all four ages for the genes each stratum found.
fig_endpoint_profiles <- function(genes, top_n = 6) {
  mean_cols <- grep("^mean_", names(genes), value = TRUE)
  genes$stratum <- paste(genes$stage, genes$sex, sep = " / ")
  top <- do.call(rbind, lapply(split(genes, genes$stratum), function(g)
    utils::head(g[order(-abs(g$effect)), ], top_n)))
  long <- do.call(rbind, lapply(mean_cols, function(col)
    data.frame(stratum = top$stratum, gene = top$gene, direction = top$direction,
               shape = top$shape, age = sub("^mean_", "", col), mean = top[[col]])))
  long <- long[is.finite(long$mean), ]
  long$age <- factor(long$age, levels = intersect(AGE_ORDER, long$age))
  ends <- long[long$age == utils::tail(levels(long$age), 1), ]
  dir_cols <- c(DIVERGING[["low"]], DIVERGING[["high"]])
  names(dir_cols) <- c(grep("^down", unique(long$direction), value = TRUE)[1],
                       grep("^up", unique(long$direction), value = TRUE)[1])
  ggplot2::ggplot(long, ggplot2::aes(.data$age, .data$mean, group = .data$gene,
                                     colour = .data$direction)) +
    ggplot2::geom_line(linewidth = 0.7, alpha = 0.9) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::geom_text(data = ends, ggplot2::aes(label = .data$gene), hjust = -0.15,
                       size = 2.6, colour = INK, check_overlap = TRUE) +
    ggplot2::facet_wrap(~ stratum, scales = "free_y") +
    ggplot2::scale_colour_manual(values = dir_cols, name = NULL) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(0.3, 0.9))) +
    ggplot2::labs(x = NULL, y = "Mean normalised expression",
                  title = "The genes that differ between 3m and 18m, followed through every age",
                  subtitle = paste0("Top ", top_n, " per stratum by effect. 9m and 12m were not used to find them: a gene that",
                                    " sits between the endpoints there is the stronger lead.")) +
    theme_bm()
}

#' Top pathways per stratum from the 3m vs 18m GSEA.
fig_endpoint_gsea <- function(gsea, top_n = 6) {
  g <- gsea[which(gsea$padj < 0.05 & gsea$independent %in% TRUE), ]
  if (!nrow(g)) return(NULL)
  g$stratum <- paste(g$stage, g$sex, sep = " / ")
  top <- do.call(rbind, lapply(split(g, g$stratum), function(x)
    utils::head(x[order(-abs(x$NES)), ], top_n)))
  top$label <- pretty_pathway(top$pathway, 46)
  top$label <- factor(top$label, levels = unique(top$label[order(top$NES)]))
  ggplot2::ggplot(top, ggplot2::aes(.data$NES, .data$label, fill = .data$NES)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, colour = MUTED) +
    ggplot2::facet_wrap(~ stratum, scales = "free_y", ncol = 2) +
    ggplot2::scale_fill_gradient2(low = DIVERGING[["low"]], mid = DIVERGING[["mid"]],
                                  high = DIVERGING[["high"]], midpoint = 0, guide = "none") +
    ggplot2::labs(x = "NES (positive: higher at 18m)", y = NULL,
                  title = "Gene sets moving between 3m and 18m, per stratum",
                  subtitle = paste0("Independent sets at padj < 0.05, top ", top_n, " by |NES|")) +
    theme_bm(base_size = 9)
}

#' Genes clearing the four-age permutation null, with the complexity check.
fig_age_trend_excess <- function(excess) {
  excess$stratum <- paste(excess$stage, excess$sex, sep = " / ")
  excess$stratum <- factor(excess$stratum, levels = excess$stratum[order(excess$n_exceeding_null)])
  excess$check <- ifelse(excess$detection_confounded %in% TRUE,
                         "genes detected trends with age too", "no complexity trend")
  ggplot2::ggplot(excess, ggplot2::aes(.data$n_exceeding_null, .data$stratum, fill = .data$check)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = .data$n_exceeding_null), hjust = -0.3, size = 3,
                       colour = INK) +
    ggplot2::scale_fill_manual(values = c("no complexity trend" = "#2a78d6",
                                          "genes detected trends with age too" = "#eda100"),
                               name = NULL) +
    ggplot2::labs(x = "Genes whose age trend beats every permuted trend", y = NULL,
                  title = "Monotonic age trends above the permutation null, per stratum",
                  subtitle = if ("ages" %in% names(excess)) "Ages per stratum are in age_trend_excess_over_null.csv" else NULL) +
    theme_bm()
}

#' Surface proteins that change with age, per stage and sex.
fig_adt_age_changes <- function(adt, top_n = 25) {
  adt$stratum <- paste(adt$stage, adt$sex, sep = " / ")
  keep <- names(sort(tapply(abs(adt$net), adt$feature, max, na.rm = TRUE), decreasing = TRUE))
  keep <- utils::head(keep, top_n)
  d <- adt[adt$feature %in% keep, ]
  d$feature <- factor(sub("^(Ms|HuMs)\\.", "", d$feature),
                      levels = sub("^(Ms|HuMs)\\.", "", rev(keep)))
  lim <- max(abs(d$net), na.rm = TRUE)
  ggplot2::ggplot(d, ggplot2::aes(.data$stratum, .data$feature, fill = .data$net)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.1f", .data$net)), size = 2.5, colour = INK) +
    ggplot2::scale_fill_gradient2(low = DIVERGING[["low"]], mid = DIVERGING[["mid"]],
                                  high = DIVERGING[["high"]], midpoint = 0,
                                  limits = c(-lim, lim), name = "log2FC\nfirst -> last age") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Surface proteins changing with age (ADT, permutation-calibrated)",
                  subtitle = "Blank: did not pass in that stratum. Protein is the depth-robust readout.") +
    theme_bm(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
                   panel.grid = ggplot2::element_blank())
}

#' Fitted expression along pseudotime for the top trajectory-model genes.
fig_trajectory_fits <- function(hits, fitted, top_n = 12) {
  if (!nrow(hits)) return(NULL)
  top <- utils::head(hits$name, top_n)
  rows <- fitted[match(top, fitted[[1]]), , drop = FALSE]
  long <- do.call(rbind, lapply(names(fitted)[-1], function(col) {
    parts <- regmatches(col, regexec("^pt([0-9.eE+-]+)_([^_]+)_(.+)$", col))[[1]]
    data.frame(gene = rows[[1]], pt = as.numeric(parts[2]), age = parts[3], sex = parts[4],
               log_mu = rows[[col]])
  }))
  long$age <- factor(long$age, levels = intersect(AGE_ORDER, long$age))
  long$gene <- factor(long$gene, levels = top)
  long$driver <- hits$driver[match(long$gene, hits$name)]
  long$label <- paste0(long$gene, " (", long$driver, ")")
  long$label <- factor(long$label, levels = unique(long$label[order(long$gene)]))
  ggplot2::ggplot(long, ggplot2::aes(.data$pt, .data$log_mu, colour = .data$age,
                                     linetype = .data$sex,
                                     group = interaction(.data$age, .data$sex))) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_wrap(~ label, scales = "free_y") +
    scale_age() +
    ggplot2::scale_linetype_manual(values = c(female = "solid", male = "22"), name = "Sex") +
    ggplot2::labs(x = "Pseudotime", y = "Fitted log expression",
                  title = "Developmental profiles that change with age or sex (trajectory model)",
                  subtitle = "Model fits on the shared pseudotime span; driver in brackets") +
    theme_bm(base_size = 9)
}

#' The interaction GSEA, per stage: which pathway trends differ between sexes.
fig_interaction_gsea <- function(gsea, top_n = 12) {
  g <- gsea[which(gsea$padj < 0.05 & gsea$independent %in% TRUE), ]
  if (!nrow(g)) return(NULL)
  top <- do.call(rbind, lapply(split(g, g$stage), function(x)
    utils::head(x[order(-abs(x$NES)), ], top_n)))
  top$label <- pretty_pathway(top$pathway, 46)
  top$label <- factor(top$label, levels = unique(top$label[order(top$NES)]))
  ggplot2::ggplot(top, ggplot2::aes(.data$NES, .data$label, fill = .data$trajectory_more_positive_in)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, colour = MUTED) +
    ggplot2::facet_wrap(~ stage, scales = "free_y") +
    scale_sex("fill", labels = c(female = "trend more positive in females",
                                 male = "trend more positive in males")) +
    ggplot2::labs(x = "NES of the sex x age interaction", y = NULL,
                  title = "Pathways whose age trend differs between the libraries",
                  subtitle = "Within stage. Cross-library: a hypothesis about slopes, not a sex difference.") +
    theme_bm(base_size = 9)
}

# --- 6. The sex contrast ----------------------------------------------------

#' The candidates the ranking says to believe, with intervals.
fig_sex_candidates <- function(cand, top_n = 25) {
  c <- cand[!cand$below_null_floor %in% TRUE & cand$stratum_depth_ok %in% c(TRUE, NA), ]
  if (!nrow(c)) return(NULL)
  c <- utils::head(c[order(-c$rank_score, -abs(c$effect_size)), ], top_n)
  c$label <- paste0(pretty_pathway(c$candidate, 50), "  [", c$stage, "]")
  c$label <- factor(c$label, levels = rev(c$label))
  c$downsample <- ifelse(c$survives_downsampling %in% c(TRUE, "TRUE"), "survives depth matching",
                         "sign flips or shrinks when matched")
  ggplot2::ggplot(c, ggplot2::aes(.data$effect_size, .data$label, colour = .data$direction)) +
    ggplot2::geom_vline(xintercept = 0, colour = MUTED) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper), linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$downsample), size = 2.3) +
    ggplot2::scale_shape_manual(values = c("survives depth matching" = 16,
                                           "sign flips or shrinks when matched" = 1), name = NULL) +
    scale_sex(name = "Higher in") +
    ggplot2::labs(x = "Standardised difference, female - male (module score)", y = NULL,
                  title = "Sex-difference candidates for orthogonal validation",
                  subtitle = "Ranked by what to believe: against the depth bias (male-high) first. Not a test of sex.",
                  caption = "Sex is confounded with library. Intervals are on cells.") +
    theme_bm(base_size = 9)
}

#' Strongest sex separation in each stage against a random split of one library.
fig_within_library_null <- function(null) {
  null$stage <- factor(null$comparison, levels = rev(intersect(STAGE_ORDER, null$comparison)))
  ggplot2::ggplot(null, ggplot2::aes(y = .data$stage)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$null_median, xend = .data$null_floor,
                                       yend = .data$stage), colour = "#86b6ef", linewidth = 3) +
    ggplot2::geom_point(ggplot2::aes(x = .data$n_observed, shape = .data$exceeds_floor),
                        colour = INK, size = 3) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                                labels = c(`TRUE` = "beyond the floor", `FALSE` = "within the floor"),
                                name = NULL) +
    ggplot2::labs(x = "Largest module separation (standardised)", y = NULL,
                  title = "Does the sex contrast beat a random split of one library?",
                  subtitle = "Bar: null median to 95th percentile from random splits of the female library. Point: male vs female.") +
    theme_bm()
}

#' Fitted z-scored age profiles of the per-sex omnibus hits.
fig_omnibus_heatmap <- function(omni, sex, top_n = 40) {
  z_cols <- intersect(AGE_ORDER, names(omni))
  if (!nrow(omni) || length(z_cols) < 2) return(NULL)
  o <- utils::head(omni, top_n)
  long <- do.call(rbind, lapply(z_cols, function(a)
    data.frame(gene = o$name, age = a, z = o[[a]])))
  long$age <- factor(long$age, levels = z_cols)
  long$gene <- factor(long$gene, levels = rev(o$name))
  ggplot2::ggplot(long, ggplot2::aes(.data$age, .data$gene, fill = .data$z)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_gradient2(low = DIVERGING[["low"]], mid = DIVERGING[["mid"]],
                                  high = DIVERGING[["high"]], midpoint = 0, name = "z") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0("Genes that change with age, ", sex, " (NB-GLM, stage-adjusted)"),
                  subtitle = "Fitted means z-scored across ages; permutation-calibrated selection") +
    theme_bm(base_size = 9) + ggplot2::theme(panel.grid = ggplot2::element_blank())
}
