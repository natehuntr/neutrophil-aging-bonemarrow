#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 8 - sex x age interaction GSEA on GO:BP.
#
#   Rscript scripts/08_gsea.R
#
# Ranks genes by the Fisher-z contrast of their per-sex age trends (see the
# header of R/gsea.R for what the sign of that statistic does and does not
# mean), runs preranked GSEA, then does two checks on the result:
#   - a sensitivity run with the abundant granule transcripts removed;
#   - per-age male-vs-female GSEA, to see whether the per-age picture agrees
#     with the interaction statistic.
#
# Reads:  results/objects/gmp_neutrophils.rds
# Writes: results/tables/gsea_sex_by_age_GOBP.csv
#         results/tables/gsea_sensitivity_comparison.csv
#         results/tables/gsea_perage_vs_interaction.csv
#         results/figures/gsea_sex_by_age_GOBP.pdf, nes_trajectories.pdf
# ---------------------------------------------------------------------------

# Run from the project root ("Rscript scripts/08_gsea.R") or from inside
# scripts/; setup.R locates the project root from either.
source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()
require_packages("fgsea", "msigdbr")

age_levels <- cfg$analysis$age_levels
gmp_neu <- read_object(cfg, "gmp_neutrophils.rds")

# Mature and immature neutrophils: the compartments with enough cells in both
# sexes at every age for a per-sex trend to be estimable.
keep <- as.character(gmp_neu$fine_neu_labels) %in% c("mature", "immature") &
  gmp_neu$age %in% age_levels
neus <- subset(gmp_neu, cells = colnames(gmp_neu)[keep])
neus <- Seurat::JoinLayers(neus)
Seurat::DefaultAssay(neus) <- "RNA"
neus <- Seurat::NormalizeData(neus, verbose = FALSE)
print(table(neus$sex, neus$age))

# --- 1. Per-sex age trends on one shared gene universe --------------------
rho_male <- compute_rho_by_sex(neus, "male", cfg)
rho_female <- compute_rho_by_sex(neus, "female", cfg)

rank_tbl <- fisher_z_contrast(rho_male, rho_female, cfg$gsea$rho_flat)
log_step("genes in shared universe: ", nrow(rank_tbl))
log_step("z_diff range: ", paste(round(range(rank_tbl$z_diff), 2), collapse = " to "))
print(sort(table(rank_tbl$pattern), decreasing = TRUE))
write_table(rank_tbl, cfg, "sex_by_age_gene_ranking.csv")

# --- 2. Preranked GSEA ----------------------------------------------------
pathways <- gobp_pathways(rank_tbl$gene, cfg)
gsea <- run_interaction_gsea(rank_tbl, pathways, cfg)
write_table(gsea$result, cfg, "gsea_sex_by_age_GOBP.csv")

reportable <- gsea$result[which(gsea$result$padj < 0.05 & gsea$result$independent), ]
log_step(nrow(reportable), " independent pathways at padj < 0.05")
print(utils::head(reportable[, c("pathway", "NES", "padj", "signal_type",
                                 "dominant_pattern", "dominant_frac")], 30))

# --- 3. Granule sensitivity ----------------------------------------------
comparison <- granule_sensitivity(gsea, pathways, cfg)
write_table(comparison, cfg, "gsea_sensitivity_comparison.csv")
log_step(sum(comparison$robust, na.rm = TRUE), " pathways survive both runs")
save_figure(plot_gsea_bars(comparison), cfg, "gsea_sex_by_age_GOBP.pdf")

# --- 4. Per-age contrasts vs the interaction ------------------------------
combined <- compare_perage_to_interaction(neus, rank_tbl, pathways, cfg)
write_table(combined, cfg, "gsea_perage_vs_interaction.csv")
print(table(combined$class))

overlap <- significant_set_overlap(combined, age_levels)
print(overlap$sizes)
print(overlap$jaccard)

nes_plot <- plot_nes_trajectories(combined, cfg)
if (!is.null(nes_plot)) save_figure(nes_plot, cfg, "nes_trajectories.pdf", height = 8)

log_step("step 8 complete")
