# ---------------------------------------------------------------------------
# Plots. Every function returns a ggplot/patchwork object so scripts can
# choose whether to display or save it.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# One visual system for every figure.
#
# Sex is identity, so it gets two categorical hues (validated for colour-vision
# deficiency: adjacent CVD delta-E 24.7). Age and maturation stage are ordered,
# so each gets a single-hue ramp, light to dark, never a rainbow. Effects that
# run both ways (up / down with age) use a blue-red diverging scale with a grey
# midpoint. Status colours are not used for data.
# ---------------------------------------------------------------------------
SEX_COLOURS   <- c(female = "#2a78d6", male = "#eb6834")
AGE_COLOURS   <- c("3m" = "#86b6ef", "9m" = "#3987e5", "12m" = "#1c5cab", "18m" = "#0d366b")
STAGE_COLOURS <- c(GMPs = "#86b6ef", proNeu = "#5598e7", preNeu = "#2a78d6",
                   immature = "#1c5cab", mature = "#104281")
DIVERGING     <- c(low = "#2a78d6", mid = "#f0efec", high = "#e34948")
INK <- "#1f2328"; MUTED <- "#6b7280"; GRID <- "#e7e9ee"

#' The shared theme: recessive grid, no chart junk, text in ink not data colour.
theme_bm <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = GRID, linewidth = 0.4),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0, colour = INK),
      plot.title = ggplot2::element_text(face = "bold", colour = INK),
      plot.title.position = "plot",
      plot.subtitle = ggplot2::element_text(colour = MUTED),
      plot.caption = ggplot2::element_text(colour = MUTED, hjust = 0),
      plot.caption.position = "plot",
      axis.text = ggplot2::element_text(colour = MUTED),
      legend.position = "top", legend.justification = "left")
}

#' Colour scales keyed by name, so a filtered plot never repaints survivors.
scale_sex <- function(aesthetic = "colour", name = "Sex", ...)
  ggplot2::scale_discrete_manual(aesthetic, values = SEX_COLOURS, name = name, ...)
scale_age <- function(aesthetic = "colour", name = "Age", ...)
  ggplot2::scale_discrete_manual(aesthetic, values = AGE_COLOURS, name = name, ...)
scale_stage <- function(aesthetic = "fill", name = "Stage", ...)
  ggplot2::scale_discrete_manual(aesthetic, values = STAGE_COLOURS, name = name, ...)

#' Readable GO:BP names: no prefix, sentence case, wrapped.
pretty_pathway <- function(x, width = 48) {
  x <- tolower(gsub("_", " ", sub("^GOBP_", "", x)))
  x <- paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  vapply(x, function(s) paste(strwrap(s, width), collapse = "\n"), character(1),
         USE.NAMES = FALSE)
}

#' The four QC distributions, split by age.
plot_qc_metrics <- function(obj, group.by = "age") {
  metrics <- c("log10GenesPerUMI", "mitoRatio", "riboRatio", "percent.hb")
  plots <- lapply(metrics, function(m)
    Seurat::VlnPlot(obj, features = m, group.by = group.by) + ggplot2::labs(x = NULL))
  patchwork::wrap_plots(plots, ncol = 2)
}

#' Doublet score distribution, coloured by the singlet/doublet call.
plot_doublet_scores <- function(obj) {
  ggplot2::ggplot(obj@meta.data,
                  ggplot2::aes(x = .data$scDblFinder.score,
                               fill = .data$scDblFinder.class)) +
    ggplot2::geom_histogram(bins = 80, alpha = 0.8, position = "identity") +
    ggplot2::scale_fill_manual(values = c(singlet = "#4393c3", doublet = "#d73027"),
                               name = "Call") +
    ggplot2::labs(title = "scDblFinder score distribution",
                  x = "Doublet score", y = "Cell count") +
    ggplot2::theme_classic(base_size = 12)
}

#' Elbow plots for whichever reductions an object carries.
plot_elbows <- function(obj, reductions = c("pca", "sct_pca", "adt_pca"), ndims = 30) {
  present <- intersect(reductions, names(obj@reductions))
  patchwork::wrap_plots(lapply(present, function(r)
    Seurat::ElbowPlot(obj, reduction = r, ndims = ndims) + ggplot2::labs(title = r)))
}

#' Pseudobulk PCA of one cell stratum: mean expression per age_sex group,
#' top variable genes, PC1 vs PC2.
#'
#' Returns NULL (with a message) when a stratum has too few groups for a PCA
#' to mean anything, so it can be mapped over every stratum without erroring.
plot_pseudobulk_pca <- function(obj, cfg, label = NULL,
                                stratum_col = "stage",
                                group_col = "age_sex",
                                n_top_genes = 2000,
                                age_levels = cfg$analysis$age_levels) {
  if (!is.null(label)) {
    cells <- colnames(obj)[as.character(obj[[stratum_col]][, 1]) == label &
                             as.character(obj$age) %in% age_levels]
    if (length(cells) < 10) {
      message("skipping '", label, "': ", length(cells), " cells")
      return(NULL)
    }
    obj <- subset(drop_graphs(obj), cells = cells)
  }

  # A cell with no hashtag call has no age_sex, and AverageExpression would
  # turn those into their own "NA" column.
  labelled <- colnames(obj)[!is.na(obj[[group_col]][, 1])]
  obj <- subset(drop_graphs(obj), cells = labelled)

  avg <- Seurat::AverageExpression(obj, group.by = group_col,
                                   assays = "RNA", layer = "data")$RNA
  if (ncol(avg) < 3) {
    message("skipping '", label, "': only ", ncol(avg), " ", group_col, " groups present")
    return(NULL)
  }

  gene_var <- matrixStats::rowVars(as.matrix(avg))
  n_use <- min(n_top_genes, sum(gene_var > 0))
  top_genes <- names(sort(stats::setNames(gene_var, rownames(avg)), decreasing = TRUE))[seq_len(n_use)]

  pca <- stats::prcomp(t(avg[top_genes, , drop = FALSE]), center = TRUE, scale. = TRUE)
  var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  pca_df <- as.data.frame(pca$x[, 1:2])
  pca_df$sample <- rownames(pca_df)
  # AverageExpression turns "3m_M" into "g3m-M": strip the "g" and split.
  parts <- strsplit(sub("^g", "", pca_df$sample), "[-_]")
  pca_df$age <- factor(vapply(parts, `[`, character(1), 1), levels = age_levels)
  pca_df$sex <- factor(vapply(parts, `[`, character(1), 2),
                       levels = c("F", "M"), labels = c("female", "male"))

  ggplot2::ggplot(pca_df, ggplot2::aes(x = .data$PC1, y = .data$PC2,
                                       colour = .data$sex, group = .data$sex)) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$age), vjust = -1, size = 3.5,
                       show.legend = FALSE) +
    ggplot2::labs(title = label,
                  x = paste0("PC1 (", var_explained[1], "%)"),
                  y = paste0("PC2 (", var_explained[2], "%)"),
                  colour = "Sex") +
    ggplot2::theme_minimal(base_size = 13)
}

#' One pseudobulk PCA per developmental stage, laid out together.
plot_pseudobulk_pca_grid <- function(obj, cfg, stratum_col = "stage") {
  labels <- cfg$analysis$stage_levels
  plots <- lapply(labels, function(l) plot_pseudobulk_pca(obj, cfg, l, stratum_col))
  names(plots) <- labels
  plots <- Filter(Negate(is.null), plots)
  patchwork::wrap_plots(plots)
}

#' |rho| distributions for the two sexes, with the selection threshold marked.
plot_rho_distributions <- function(trend_male, trend_female, rho_cutoff = 0.1) {
  df <- rbind(data.frame(rho = abs(trend_male$rho), sex = "male"),
              data.frame(rho = abs(trend_female$rho), sex = "female"))
  ggplot2::ggplot(df, ggplot2::aes(.data$rho, fill = .data$sex)) +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::geom_vline(xintercept = rho_cutoff, linetype = "dashed") +
    ggplot2::labs(x = "|Spearman rho| (age trend)",
                  title = "Age-trend effect sizes by sex") +
    ggplot2::theme_minimal()
}

#' Bar chart of the strongest GSEA hits.
plot_gsea_bars <- function(comparison, n = 25) {
  top <- utils::head(comparison, n)
  top$label <- tolower(gsub("_", " ", sub("^GOBP_", "", top$pathway)))

  ggplot2::ggplot(top, ggplot2::aes(x = stats::reorder(.data$label, .data$NES_full),
                                    y = .data$NES_full, fill = .data$NES_full > 0)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#2C6E9B", `FALSE` = "#B5482A"),
                               labels = c(`TRUE` = "more positive in males",
                                          `FALSE` = "more positive in females"),
                               name = "Ageing trajectory") +
    ggplot2::labs(x = NULL, y = "NES (sex x age interaction)",
                  title = "GO:BP pathways with sex-differential ageing") +
    ggplot2::theme_minimal(base_size = 11)
}

#' NES per age for the pathways the interaction test flagged.
plot_nes_trajectories <- function(classified, cfg, n = 30,
                                  age_levels = cfg$analysis$age_levels) {
  age_levels <- age_levels[paste0("NES_", age_levels) %in% names(classified)]
  top <- utils::head(classified[which(classified$class == "sex-differential ageing (supported)"), ], n)
  if (!nrow(top)) { message("no supported hits to plot"); return(NULL) }

  plot_df <- tidyr::pivot_longer(
    top[, c("pathway", paste0("NES_", age_levels))],
    -"pathway", names_to = "age", values_to = "NES")
  plot_df$age <- factor(sub("^NES_", "", plot_df$age), levels = age_levels)
  plot_df$label <- tolower(gsub("_", " ", sub("^GOBP_", "", plot_df$pathway)))

  ggplot2::ggplot(plot_df, ggplot2::aes(.data$age,
                                        stats::reorder(.data$label, .data$NES),
                                        fill = .data$NES)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_gradient2(low = "#B5482A", mid = "grey95", high = "#2C6E9B",
                                  midpoint = 0, name = "NES\n(male vs female)") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Sex difference by age, for pathways flagged by the interaction test") +
    ggplot2::theme_minimal(base_size = 10)
}

#' Moran's I at the first age against every later age, per gene.
#'
#' Points on the diagonal are genes as structured along the trajectory at both
#' ages; off-diagonal genes gained or lost structure. This replaces a Venn of
#' per-age significance lists, which compared power rather than biology.
plot_moran_by_age <- function(by_age, age_levels, title = "", n_label = 12) {
  i_cols <- paste0("I_", age_levels)
  i_cols <- i_cols[i_cols %in% names(by_age)]
  if (length(i_cols) < 2) return(NULL)
  ref <- i_cols[1]
  long <- do.call(rbind, lapply(i_cols[-1], function(col) {
    keep <- is.finite(by_age[[ref]]) & is.finite(by_age[[col]])
    data.frame(gene = by_age$gene[keep], x = by_age[[ref]][keep],
               y = by_age[[col]][keep], age = sub("^I_", "", col))
  }))
  long$age <- factor(long$age, levels = sub("^I_", "", i_cols[-1]))
  long$delta <- long$y - long$x
  labels <- do.call(rbind, lapply(split(long, long$age), function(d)
    utils::head(d[order(-abs(d$delta)), ], n_label)))
  lim <- range(c(long$x, long$y), na.rm = TRUE)

  ggplot2::ggplot(long, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    ggplot2::geom_point(alpha = 0.25, size = 0.6, colour = "grey30") +
    ggplot2::geom_point(data = labels, colour = "#2a78d6", size = 1.6) +
    ggplot2::geom_text(data = labels, ggplot2::aes(label = .data$gene),
                       size = 2.6, vjust = -0.7, colour = "grey15", check_overlap = TRUE) +
    ggplot2::facet_wrap(~ age, labeller = ggplot2::labeller(age = function(a)
      paste(sub("^I_", "", ref), "vs", a))) +
    ggplot2::coord_equal(xlim = lim, ylim = lim) +
    ggplot2::labs(x = paste("Moran's I at", sub("^I_", "", ref)),
                  y = "Moran's I at the later age", title = title,
                  subtitle = "Structure along each age's own trajectory; labelled: largest changes") +
    theme_bm()
}

#' Stacked stage composition per age, one panel per sex.
plot_stage_composition <- function(composition) {
  ggplot2::ggplot(composition,
                  ggplot2::aes(x = .data$age, y = .data$proportion, fill = .data$stage)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_wrap(~ sex) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = NULL, y = "Share of granulocytic cells", fill = "Stage",
                  title = "Developmental stage composition by age") +
    ggplot2::theme_minimal(base_size = 12)
}

#' Each stage's proportion against age, with Wilson intervals.
#'
#' The interval is the sampling error on the cells actually captured. It says
#' nothing about how much this would vary between mice, which is the
#' uncertainty that matters and which this design cannot estimate.
plot_stage_trends <- function(composition) {
  ggplot2::ggplot(composition,
                  ggplot2::aes(x = .data$age, y = .data$proportion,
                               colour = .data$sex, group = .data$sex)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high),
                           width = 0.15) +
    ggplot2::facet_wrap(~ stage, scales = "free_y") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = NULL, y = "Share of granulocytic cells", colour = "Sex",
                  title = "Stage proportions across age",
                  caption = "Intervals are Wilson binomial CIs on captured cells, not between-animal variation") +
    ggplot2::theme_minimal(base_size = 11)
}

#' ECDF of a per-cell quantity by age, one panel per sex.
#'
#' The ECDF is the right display for these comparisons: the KS statistic is the
#' largest vertical gap between two of these curves, and a rightward shift is
#' cells sitting further along the trajectory.
plot_distribution_by_age <- function(values, meta, cfg, xlab,
                                     age_levels = cfg$analysis$age_levels) {
  df <- data.frame(value = values,
                   age = factor(as.character(meta$age), levels = age_levels),
                   sex = meta$sex)
  df <- df[is.finite(df$value) & !is.na(df$age), ]

  ecdf_plot <- ggplot2::ggplot(df, ggplot2::aes(.data$value, colour = .data$age)) +
    ggplot2::stat_ecdf(linewidth = 0.7) +
    ggplot2::facet_wrap(~ sex) +
    ggplot2::labs(x = xlab, y = "Cumulative share of cells", colour = "Age") +
    ggplot2::theme_minimal(base_size = 12)

  density_plot <- ggplot2::ggplot(df, ggplot2::aes(.data$value, fill = .data$age)) +
    ggplot2::geom_density(alpha = 0.35, colour = NA) +
    ggplot2::facet_wrap(~ sex) +
    ggplot2::labs(x = xlab, y = "Density", fill = "Age") +
    ggplot2::theme_minimal(base_size = 12)

  patchwork::wrap_plots(ecdf_plot, density_plot, ncol = 1)
}

# ---------------------------------------------------------------------------
# Confound diagnostics.
#
# Three figures that between them say whether any cross-sex number in this
# project is interpretable. They belong at the front of a report, not in an
# appendix: everything downstream is conditional on them.
# ---------------------------------------------------------------------------

#' Per-library UMI and gene-count distributions.
#'
#' The depth confound, shown rather than quoted. If the two libraries do not
#' overlap, no amount of normalisation makes a cross-sex comparison safe.
plot_depth_distributions <- function(obj, cfg, group_col = "sex") {
  df <- data.frame(
    group = as.character(obj[[group_col]][, 1]),
    umis = obj$nCount_RNA,
    genes = obj$nFeature_RNA
  )
  df <- df[stats::complete.cases(df), ]

  medians <- stats::aggregate(cbind(umis, genes) ~ group, df, stats::median)
  ratio <- max(medians$umis) / min(medians$umis)

  umi_plot <- ggplot2::ggplot(df, ggplot2::aes(.data$umis, fill = .data$group)) +
    ggplot2::geom_density(alpha = 0.4, colour = NA) +
    ggplot2::geom_vline(data = medians,
                        ggplot2::aes(xintercept = .data$umis, colour = .data$group),
                        linetype = "dashed", show.legend = FALSE) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "UMIs per cell (log scale)", y = "Density", fill = NULL,
                  title = sprintf("Median depth differs %.2fx between libraries", ratio),
                  subtitle = "Sex is confounded with library, so this is the sex contrast's floor") +
    ggplot2::theme_minimal(base_size = 12)

  gene_plot <- ggplot2::ggplot(df, ggplot2::aes(.data$genes, fill = .data$group)) +
    ggplot2::geom_density(alpha = 0.4, colour = NA) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Genes detected per cell (log scale)", y = "Density", fill = NULL) +
    ggplot2::theme_minimal(base_size = 12)

  patchwork::wrap_plots(umi_plot, gene_plot, ncol = 1)
}

#' Detection rate in one group against the other, with the depth expectation.
#'
#' The most persuasive single exhibit available here. Under a pure depth
#' difference the points follow a curve set by the depth ratio and nothing
#' else; genes departing from it are the only ones a difference can be claimed
#' for. Points below the diagonal that were called "up" in the shallower group
#' are the incoherence made visible.
plot_detection_scatter <- function(detection, groups, highlight = character()) {
  cols <- paste0("pct_", groups)
  df <- data.frame(x = detection[[cols[1]]], y = detection[[cols[2]]],
                   gene = detection$gene)
  df <- df[stats::complete.cases(df), ]

  # Expected relationship if the only difference is capture probability:
  # 1 - (1 - p)^r for depth ratio r, fitted on the observed pairs.
  fit_ratio <- tryCatch(stats::optimise(function(r)
    sum((1 - (1 - df$x)^r - df$y)^2), c(0.1, 10))$minimum, error = function(e) NA_real_)
  curve <- data.frame(x = seq(0, 1, length.out = 200))
  curve$y <- 1 - (1 - curve$x)^fit_ratio

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_point(alpha = 0.15, size = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    ggplot2::geom_line(data = curve, colour = "#B5482A", linewidth = 0.8) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = paste("Fraction of", groups[1], "cells detecting the gene"),
                  y = paste("Fraction of", groups[2], "cells detecting the gene"),
                  title = "Detection rates track depth, not biology",
                  subtitle = sprintf(paste("Red: expected under a pure capture difference",
                                           "(fitted ratio %.2f). Dotted: equality."),
                                     fit_ratio)) +
    ggplot2::theme_minimal(base_size = 12)

  if (length(highlight)) {
    marked <- df[df$gene %in% highlight, ]
    p <- p + ggplot2::geom_point(data = marked, colour = "#2C6E9B", size = 1.8) +
      ggplot2::geom_text(data = marked, ggplot2::aes(label = .data$gene),
                         size = 3, vjust = -0.8, colour = "#2C6E9B")
  }
  p
}

#' Retained-cell fraction through QC, by sex and stage.
plot_retention <- function(retention, sex_col = "sex", stage_col = "stage") {
  ggplot2::ggplot(retention, ggplot2::aes(x = .data[[stage_col]],
                                          y = .data$retained_fraction,
                                          fill = .data[[sex_col]])) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = .data$n_after),
                       position = ggplot2::position_dodge(width = 0.8),
                       vjust = -0.4, size = 3) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                limits = c(0, 1.05)) +
    ggplot2::labs(x = NULL, y = "Cells retained through QC", fill = NULL,
                  title = "QC retention by stage and library",
                  subtitle = "Unequal retention converts a depth difference into a composition one",
                  caption = "Labels are cell counts after filtering") +
    ggplot2::theme_minimal(base_size = 12)
}

#' All three confound diagnostics, written together.
write_confound_diagnostics <- function(obj, cfg, retention = NULL,
                                       detection = NULL, prefix = "confound") {
  save_figure(plot_depth_distributions(obj, cfg), cfg,
              paste0(prefix, "_depth_distributions.pdf"), width = 8, height = 7)

  if (!is.null(detection))
    save_figure(plot_detection_scatter(detection, cfg$analysis$sex_levels), cfg,
                paste0(prefix, "_detection_scatter.pdf"), width = 7, height = 7)

  if (!is.null(retention))
    save_figure(plot_retention(retention), cfg,
                paste0(prefix, "_qc_retention.pdf"), width = 9, height = 5)

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Protein gate (step 12)
# ---------------------------------------------------------------------------

#' The gate itself, drawn: Ly6G against CD11b, then CXCR2 against Ly6G within
#' myeloid cells, one row per library, with the thresholds used.
plot_protein_gate <- function(cells, thresholds) {
  cells <- cells[!is.na(cells$protein_class), ]
  thr <- function(sex, m) {
    t <- thresholds[thresholds$sex == sex & thresholds$marker == m, "threshold"]
    if (length(t)) t[1] else NA_real_
  }
  lines <- do.call(rbind, lapply(unique(cells$sex), function(s)
    data.frame(sex = s, cd11b = thr(s, "CD11b"), ly6g = thr(s, "Ly6G"),
               cxcr2 = thr(s, "CXCR2"))))
  a <- ggplot2::ggplot(cells, ggplot2::aes(.data$CD11b, .data$Ly6G)) +
    ggplot2::geom_bin2d(bins = 90) +
    ggplot2::scale_fill_gradient(low = "#cde2fb", high = "#0d366b", trans = "log10",
                                 name = "Cells") +
    ggplot2::geom_vline(data = lines, ggplot2::aes(xintercept = .data$cd11b), colour = INK) +
    ggplot2::geom_hline(data = lines, ggplot2::aes(yintercept = .data$ly6g), colour = INK) +
    ggplot2::facet_wrap(~ sex) +
    ggplot2::labs(x = "CD11b (log, isotype-centred)", y = "Ly6G",
                  title = "Protein gate: myeloid cells, then Ly6G") +
    theme_bm()
  my <- cells[cells$CD11b > lines$cd11b[match(cells$sex, lines$sex)] &
                is.finite(cells$CXCR2), ]
  if (!nrow(my)) return(a)
  b <- ggplot2::ggplot(my, ggplot2::aes(.data$Ly6G, .data$CXCR2)) +
    ggplot2::geom_bin2d(bins = 90) +
    ggplot2::scale_fill_gradient(low = "#cde2fb", high = "#0d366b", trans = "log10",
                                 name = "Cells") +
    ggplot2::geom_vline(data = lines, ggplot2::aes(xintercept = .data$ly6g), colour = INK) +
    ggplot2::geom_hline(data = lines[is.finite(lines$cxcr2), ],
                        ggplot2::aes(yintercept = .data$cxcr2), colour = INK) +
    ggplot2::facet_wrap(~ sex) +
    ggplot2::labs(x = "Ly6G", y = "CXCR2 (CD182)",
                  title = "Within CD11b+ cells: CXCR2 separates mature-like from immature-like") +
    theme_bm()
  patchwork::wrap_plots(a, b, ncol = 1)
}

#' Protein-stage composition by age, before and after RNA QC, side by side.
#'
#' If the two rows of panels differ for one library, RNA QC changed its stage
#' mix -- and the RNA-staged composition inherited that change.
plot_protein_composition <- function(comp, stage_levels) {
  comp$stage <- factor(as.character(comp$stage), levels = stage_levels)
  comp$age <- factor(as.character(comp$age), levels = names(AGE_COLOURS))
  cols <- stats::setNames(c("#86b6ef", "#2a78d6", "#104281"), stage_levels)
  ggplot2::ggplot(comp, ggplot2::aes(.data$age, .data$proportion, fill = .data$stage)) +
    ggplot2::geom_col(width = 0.7, colour = "white", linewidth = 0.4) +
    ggplot2::facet_grid(cell_set ~ sex) +
    ggplot2::scale_fill_manual(values = cols, name = "Protein stage") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = NULL, y = "Share of protein-gated granulocytes",
                  title = "Stage mix from surface protein, before and after RNA QC",
                  subtitle = "Top: every singlet. Bottom: the cells the RNA analysis kept.") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    theme_bm()
}

#' Fraction of each protein-defined class that survived RNA QC, by age.
plot_protein_retention <- function(retention, stage_levels) {
  r <- retention[retention$protein_class %in% stage_levels, ]
  r$protein_class <- factor(r$protein_class, levels = stage_levels)
  r$age <- factor(r$age, levels = names(AGE_COLOURS))
  ggplot2::ggplot(r, ggplot2::aes(.data$age, .data$retained_fraction,
                                  colour = .data$sex, group = .data$sex)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high),
                             size = 0.35) +
    ggplot2::facet_wrap(~ protein_class) +
    scale_sex() +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                limits = c(0, 1.08)) +
    ggplot2::labs(x = NULL, y = "Kept by RNA QC",
                  title = "How much of each protein-defined stage RNA QC kept",
                  subtitle = "A fall with age in one library means QC, not the marrow, removed those cells",
                  caption = "Intervals: Wilson, on cells. Cell counts: protein_gate_retention.csv.") +
    theme_bm()
}
