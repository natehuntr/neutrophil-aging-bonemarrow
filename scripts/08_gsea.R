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
# AGES: per stage, the ages from gsea.age_levels at which both sexes clear
# gates.min_cells_per_stratum. The table rows say which.
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

# gsea.age_levels is the CANDIDATE set; each stage then keeps the ages where
# BOTH sexes clear the cell gate, since the interaction compares the two
# sexes' trends over the same ages.
candidate_ages <- cfg$gsea$age_levels %||% cfg$analysis$age_levels
args <- commandArgs(trailingOnly = TRUE)
stages <- if (length(args)) args else cfg$gsea$stages

gmp_neu <- read_object(cfg, "gmp_neutrophils.rds")
gsea_assay <- matched_assay_for(gmp_neu, cfg, step = 8)
require_metadata(gmp_neu, c("stage", "age", "sex"), context = "step 8")

for (stage in stages) {
  log_step("=================== ", stage, " ===================")

  in_stage <- !is.na(gmp_neu$stage) & as.character(gmp_neu$stage) == stage &
    gmp_neu$age %in% candidate_ages
  # Requiring every sex x age cell to clear the gate retired both default
  # stages (male 12m holds 24 mature and 34 immature cells), so this step had
  # produced nothing since the gate was introduced. Thin ages are dropped
  # instead; a trend needs at least two.
  age_levels <- shared_usable_ages(gmp_neu$age[in_stage], gmp_neu$sex[in_stage],
                                   candidate_ages, cfg$gates$min_cells_per_stratum)
  if (length(age_levels) < 2) {
    log_step("skipping ", stage, ": fewer than two ages clear the ",
             cfg$gates$min_cells_per_stratum, "-cell gate in both sexes")
    next
  }
  log_step(stage, ": ages used -- ", paste(age_levels, collapse = ", "))

  neus <- select_cells(gmp_neu, list(
    "stage is this stage"     = as.character(gmp_neu$stage) == stage,
    "age clears the gate in both sexes" = gmp_neu$age %in% age_levels
  ), context = paste(stage, "neutrophils"))
  print(table(neus$sex, neus$age))

  neus <- join_layers(neus)
  Seurat::DefaultAssay(neus) <- gsea_assay
  neus <- Seurat::NormalizeData(neus, assay = gsea_assay, verbose = FALSE)

  # --- 1. Per-sex age trends on one shared gene universe ------------------
  # On the step's own assay. These defaulted to "RNA", so the trends were
  # computed on unmatched counts while every log line said matched.
  rho_male <- compute_rho_by_sex(neus, "male", cfg, assay = gsea_assay)
  rho_female <- compute_rho_by_sex(neus, "female", cfg, assay = gsea_assay)

  rank_tbl <- fisher_z_contrast(rho_male, rho_female, cfg$gsea$rho_flat)
  # Sex-chromosome genes are guaranteed extremes of anything that contrasts
  # the two libraries and say nothing about granulopoiesis. They are reported
  # on their own as a check that the sex assignment is right, then removed
  # from the ranking.
  sex_chr <- rank_tbl[rank_tbl$gene %in% SEX_CHR_GENES, ]
  if (nrow(sex_chr)) {
    log_step("sex-chromosome genes, removed from the ranking:")
    print(as.data.frame(sex_chr[, c("gene", "rho_M", "rho_F", "z_diff")]), row.names = FALSE)
  }
  rank_tbl <- rank_tbl[!rank_tbl$gene %in% SEX_CHR_GENES, ]
  log_step("genes in shared universe: ", nrow(rank_tbl))
  log_step("z_diff range: ", paste(round(range(rank_tbl$z_diff), 2), collapse = " to "))
  print(sort(table(rank_tbl$pattern), decreasing = TRUE))
  write_table(rank_tbl, cfg, sprintf("sex_by_age_gene_ranking_%s.csv", stage))

  # --- 2. Preranked GSEA --------------------------------------------------
  pathways <- gobp_pathways(rank_tbl$gene, cfg)
  gsea <- run_interaction_gsea(rank_tbl, pathways, cfg)
  # The confound has a direction, so a pathway pointing with it and one
  # pointing against it are not equally believable.
  annotated <- annotate_bias(gsea$result, gsea$result$trajectory_more_positive_in, cfg)
  annotated$orthogonal_assay <- suggest_orthogonal_assay(annotated$pathway)
  write_table(strip_inferential_columns(annotated, keep = "padj"), cfg,
              sprintf("gsea_sex_by_age_GOBP_%s.csv", stage))

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
  combined <- compare_perage_to_interaction(neus, rank_tbl, pathways, cfg,
                                            age_levels = age_levels,
                                            assay = gsea_assay)
  write_table(combined, cfg, sprintf("gsea_perage_vs_interaction_%s.csv", stage))
  print(table(combined$class))

  overlap <- significant_set_overlap(combined, age_levels)
  print(overlap$sizes)
  print(overlap$jaccard)

  nes_plot <- plot_nes_trajectories(combined, cfg, age_levels = age_levels)
  if (!is.null(nes_plot))
    save_figure(nes_plot + ggplot2::labs(subtitle = paste(stage, "neutrophils")),
                cfg, sprintf("nes_trajectories_%s.pdf", stage), height = 8)

  rm(neus, rank_tbl, pathways, gsea, comparison, combined)
  gc(verbose = FALSE)
}

log_step("step 8 complete")
