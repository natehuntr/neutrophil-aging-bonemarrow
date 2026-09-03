# ---------------------------------------------------------------------------
# Plots. Every function returns a ggplot/patchwork object so scripts can
# choose whether to display or save it.
# ---------------------------------------------------------------------------

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
                                stratum_col = "fine_neu_labels",
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
    obj <- subset(obj, cells = cells)
  }

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
plot_pseudobulk_pca_grid <- function(obj, cfg, stratum_col = "fine_neu_labels") {
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
plot_nes_trajectories <- function(classified, cfg, n = 30) {
  age_levels <- cfg$analysis$age_levels
  top <- utils::head(classified[classified$class == "sex-differential ageing (supported)", ], n)
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

#' Venn of the significant trajectory genes from each age.
plot_gene_venn <- function(gene_lists, title = "") {
  ggVennDiagram::ggVennDiagram(gene_lists) +
    ggplot2::scale_fill_gradient(low = "white", high = "lightblue") +
    ggplot2::labs(title = title)
}
