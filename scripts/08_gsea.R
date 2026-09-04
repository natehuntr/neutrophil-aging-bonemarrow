#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 8 - sex x age interaction GSEA on GO:BP, within a maturation stage.
#
#   Rscript scripts/08_gsea.R [stage ...]
#
# Ranks genes by the Fisher-z contrast of their per-sex age trends (see the
# header of R/gsea.R for what the sign of that statistic does and does not
# mean), runs preranked GSEA, then does two checks on the result:
#   - a sensitivity run with the abundant granule transcripts removed;
#   - per-age male-vs-female GSEA, to see whether the per-age picture agrees
#     with the interaction statistic.
#
# STRATIFICATION: each stage in gsea.stages is analysed on its own. Pooling
# stages would let a composition shift masquerade as a within-cell age trend --
# if the immature share rises with age, every gene that is higher in immature
# cells acquires an age trend without changing in either stage. Stratifying is
# what makes the ranking a statement about cells rather than about the mix.
#
# Reads:  results/objects/gmp_neutrophils.rds
# Writes: results/tables/gsea_sex_by_age_GOBP_<stage>.csv
#         results/tables/gsea_sensitivity_comparison_<stage>.csv
#         results/tables/gsea_perage_vs_interaction_<stage>.csv
#         results/tables/sex_by_age_gene_ranking_<stage>.csv
#         results/figures/gsea_sex_by_age_GOBP_<stage>.pdf,
#                         nes_trajectories_<stage>.pdf
# ---------------------------------------------------------------------------

# Run from the project root ("Rscript scripts/08_gsea.R") or from inside
# scripts/; setup.R locates the project root from either.
source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()
require_packages("fgsea", "msigdbr")

age_levels <- cfg$analysis$age_levels
args <- commandArgs(trailingOnly = TRUE)
stages <- if (length(args)) args else cfg$gsea$stages

gmp_neu <- read_object(cfg, "gmp_neutrophils.rds")
require_metadata(gmp_neu, c("fine_neu_labels", "age", "sex"), context = "step 8")

for (stage in stages) {
  log_step("=================== ", stage, " ===================")

  neus <- select_cells(gmp_neu, list(
    "fine_neu_labels is this stage"     = as.character(gmp_neu$fine_neu_labels) == stage,
    "age is one of analysis.age_levels" = gmp_neu$age %in% age_levels
  ), context = paste(stage, "neutrophils"))

  group_sizes <- table(neus$sex, neus$age)
  print(group_sizes)
  if (any(dim(group_sizes) < c(2, length(age_levels))) || any(group_sizes < 20)) {
    log_step("skipping ", stage, ": every sex x age group needs at least 20 cells")
    next
  }

  neus <- join_layers(neus)
  Seurat::DefaultAssay(neus) <- "RNA"
  neus <- Seurat::NormalizeData(neus, verbose = FALSE)

  # --- 1. Per-sex age trends on one shared gene universe ------------------
  rho_male <- compute_rho_by_sex(neus, "male", cfg)
  rho_female <- compute_rho_by_sex(neus, "female", cfg)

  rank_tbl <- fisher_z_contrast(rho_male, rho_female, cfg$gsea$rho_flat)
  log_step("genes in shared universe: ", nrow(rank_tbl))
  log_step("z_diff range: ", paste(round(range(rank_tbl$z_diff), 2), collapse = " to "))
  print(sort(table(rank_tbl$pattern), decreasing = TRUE))
  write_table(rank_tbl, cfg, sprintf("sex_by_age_gene_ranking_%s.csv", stage))

  # --- 2. Preranked GSEA --------------------------------------------------
  pathways <- gobp_pathways(rank_tbl$gene, cfg)
  gsea <- run_interaction_gsea(rank_tbl, pathways, cfg)
  write_table(gsea$result, cfg, sprintf("gsea_sex_by_age_GOBP_%s.csv", stage))

  reportable <- gsea$result[which(gsea$result$padj < 0.05 & gsea$result$independent), ]
  log_step(nrow(reportable), " independent pathways at padj < 0.05")
  if (nrow(reportable))
    print(utils::head(reportable[, c("pathway", "NES", "padj", "signal_type",
                                     "dominant_pattern", "dominant_frac")], 30))

  # --- 3. Granule sensitivity ---------------------------------------------
  comparison <- granule_sensitivity(gsea, pathways, cfg)
  write_table(comparison, cfg, sprintf("gsea_sensitivity_comparison_%s.csv", stage))
  log_step(sum(comparison$robust, na.rm = TRUE), " pathways survive both runs")
  save_figure(plot_gsea_bars(comparison) +
                ggplot2::labs(subtitle = paste(stage, "neutrophils")),
              cfg, sprintf("gsea_sex_by_age_GOBP_%s.pdf", stage))

  # --- 4. Per-age contrasts vs the interaction ----------------------------
  combined <- compare_perage_to_interaction(neus, rank_tbl, pathways, cfg)
  write_table(combined, cfg, sprintf("gsea_perage_vs_interaction_%s.csv", stage))
  print(table(combined$class))

  overlap <- significant_set_overlap(combined, age_levels)
  print(overlap$sizes)
  print(overlap$jaccard)

  nes_plot <- plot_nes_trajectories(combined, cfg)
  if (!is.null(nes_plot))
    save_figure(nes_plot + ggplot2::labs(subtitle = paste(stage, "neutrophils")),
                cfg, sprintf("nes_trajectories_%s.pdf", stage), height = 8)

  rm(neus, rank_tbl, pathways, gsea, comparison, combined)
  gc(verbose = FALSE)
}

log_step("step 8 complete")
