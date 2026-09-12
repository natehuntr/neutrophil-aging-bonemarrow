#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 10 - what differs between 3m and 18m, and do 9m and 12m follow?
#
#   Rscript scripts/10_endpoint_contrast.R
#
# Step 6 asks whether expression trends monotonically across all four ages,
# which spends the stratum's power on the shape of the trend and requires
# every timepoint to clear the cell gate. With 24 male 12m mature cells that
# retires eight of ten stage x sex strata before anything is computed.
#
# This step splits the question in two. Genes are FOUND at the extremes, 3m
# vs 18m, where the difference is largest and the cells are most numerous --
# 8 of 10 strata clear the gate on those two ages. Then 9m and 12m are read
# off for exactly those genes, in the same stage, as a shape check: a gene
# that interpolates across the timepoints the contrast never used is a better
# aging hypothesis than one that differs at the extremes and wanders between.
#
# The shape check reuses the discovery cells and several intermediate strata
# are small, so it REORDERS candidates. It is not independent confirmation,
# and the tables say so.
#
# Reads:  results/objects/gmp_neutrophils.rds
# Writes: results/tables/endpoint_<stage>_<sex>.csv
#         results/tables/endpoint_gsea_<stage>_<sex>.csv
#         results/tables/endpoint_summary.csv
# ---------------------------------------------------------------------------

source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()

gmp_neu <- read_object(cfg, "gmp_neutrophils.rds")
require_metadata(gmp_neu, c("stage", "age", "sex"), context = "step 10")
assay <- matched_assay_for(gmp_neu, cfg, step = 10)

reference <- cfg$endpoint$reference
endpoint <- cfg$endpoint$endpoint
log_step("discovery contrast: ", reference, " vs ", endpoint,
         "; shape read off the remaining ages in analysis.age_levels")

summary_rows <- list()

for (stage in cfg$analysis$stage_levels) {
  for (sex in cfg$analysis$sex_levels) {
    label <- paste(stage, sex, sep = "/")
    cells <- colnames(gmp_neu)[!is.na(gmp_neu$stage) &
                                 as.character(gmp_neu$stage) == stage &
                                 as.character(gmp_neu$sex) == sex]
    if (length(cells) < cfg$gates$min_cells_per_stratum) {
      log_step("  skipping ", label, ": ", length(cells), " cells in total")
      next
    }
    log_step("=== ", label, " ===")
    sub <- subset(drop_graphs(gmp_neu), cells = cells)

    res <- endpoint_contrast(sub, cfg, label = label, assay = assay)
    if (is.null(res)) next

    n_hits <- if (is.null(res$genes)) 0L else nrow(res$genes)
    shapes <- if (n_hits) table(res$genes$shape) else integer()

    if (n_hits) {
      out <- res$genes
      out$stage <- stage
      out$sex <- sex
      write_table(strip_inferential_columns(out), cfg,
                  sprintf("endpoint_%s_%s.csv", stage, sex))
    }

    # Set-level enrichment over the FULL ranking, not just the genes above the
    # threshold. Step 8's interaction GSEA needs a per-sex trend across ages
    # and cannot run in these strata at all; this one only needs the two
    # endpoints, so it covers every stratum the discovery contrast covers.
    gsea <- endpoint_gsea(res$effect, cfg, label = label)
    n_paths <- 0L
    if (!is.null(gsea) && nrow(gsea)) {
      gsea$stage <- stage
      gsea$sex <- sex
      write_table(strip_inferential_columns(gsea, keep = "padj"), cfg,
                  sprintf("endpoint_gsea_%s_%s.csv", stage, sex))
      top <- gsea[which(gsea$padj < 0.05 & gsea$independent), ]
      n_paths <- nrow(top)
      if (n_paths)
        print(utils::head(as.data.frame(
          top[, c("pathway", "NES", "padj", "direction", "leadingEdge_n",
                  "orthogonal_assay")]), 10), row.names = FALSE)
    }

    summary_rows[[label]] <- data.frame(
      stage = stage, sex = sex,
      n_reference = as.integer(res$n_by_age[reference] %||% 0L),
      n_endpoint = as.integer(res$n_by_age[endpoint] %||% 0L),
      permutation_threshold = res$threshold,
      n_genes = n_hits,
      n_monotonic = as.integer(shapes["monotonic"] %||% 0L),
      n_pathways = n_paths,
      row.names = NULL)
  }
}

if (!length(summary_rows)) {
  log_step("no stratum had enough cells at both ", reference, " and ", endpoint)
} else {
  summary_tbl <- do.call(rbind, summary_rows)
  summary_tbl$frac_monotonic <- ifelse(summary_tbl$n_genes > 0,
                                       summary_tbl$n_monotonic / summary_tbl$n_genes,
                                       NA_real_)
  write_table(summary_tbl, cfg, "endpoint_summary.csv")
  log_step(reference, " vs ", endpoint, " by stratum:")
  print(summary_tbl, row.names = FALSE)
  log_step("n_monotonic counts genes whose 9m and 12m means sit between the ",
           "endpoints in age order. Those cells were used to find the genes, ",
           "so this ranks candidates -- it does not validate them.")
}

log_step("step 10 complete")
