#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 12 - stage composition from surface protein, before and after RNA QC.
#
#   Rscript scripts/12_protein_gate.R
#
# Answers one question the RNA analysis cannot: does the male library lose
# maturing neutrophils with age, or does RNA QC remove them? See the header of
# R/protein_gate.R for the reasoning and the limits.
#
# Reads:  results/objects/<sex>_protein_preqc.rds   (written by step 1)
#         results/objects/gmp_neutrophils.rds       (optional: RNA stages, to
#                                                    cross the two calls)
# Writes: results/tables/protein_gate_thresholds.csv
#         results/tables/protein_class_counts.csv
#         results/tables/protein_stage_composition.csv
#         results/tables/protein_stage_age_trends.csv
#         results/tables/protein_gate_retention.csv     <- the decisive table
#         results/tables/protein_vs_rna_stage.csv
#         results/figures/protein_gate.pdf, protein_stage_composition.pdf,
#                         protein_gate_retention.pdf
# ---------------------------------------------------------------------------

source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()

age_levels <- cfg$analysis$age_levels
keys <- names(cfg$samples)
paths <- vapply(keys, function(k) object_path(cfg, paste0(k, "_protein_preqc.rds")),
                character(1))
if (!all(file.exists(paths)))
  stop("no pre-QC protein record for: ", paste(keys[!file.exists(paths)], collapse = ", "),
       "\nStep 1 writes <sex>_protein_preqc.rds from this version on; re-run step 1 ",
       "(and the steps after it) to produce it.")

# ===========================================================================
# 1. Gate each library separately
# ===========================================================================
gated <- lapply(stats::setNames(keys, keys), function(k) {
  record <- readRDS(paths[[k]])
  log_step("=== ", k, ": ", ncol(record$adt), " hashtagged barcodes ===")
  print(table(record$meta$qc_status))
  g <- gate_library(record, cfg)
  log_step("  markers: ", paste(names(g$markers), g$markers, sep = "=", collapse = ", "))
  g$thresholds$sex <- record$sex
  g
})

has_cxcr2 <- all(vapply(gated, `[[`, logical(1), "has_cxcr2"))
stage_levels <- protein_stage_levels(has_cxcr2)
if (!has_cxcr2)
  log_step("CXCR2 is not on the panel in every library: Ly6G+ cells are split ",
           "into Ly6G-int / Ly6G-hi by the Ly6G density valley instead")

cells <- do.call(rbind, lapply(gated, `[[`, "cells"))
thresholds <- do.call(rbind, lapply(gated, `[[`, "thresholds"))
write_table(thresholds, cfg, "protein_gate_thresholds.csv")
log_step("thresholds (per library -- each was stained separately):")
print(thresholds, row.names = FALSE)

class_levels <- unique(c(stage_levels, PROTEIN_CLASSES[-(1:3)]))
counts <- as.data.frame(table(sex = cells$sex,
                              age = factor(cells$age, levels = age_levels),
                              protein_class = factor(cells$protein_class, levels = class_levels),
                              passed_rna_qc = cells$passed_rna_qc),
                        responseName = "n", stringsAsFactors = FALSE)
write_table(counts, cfg, "protein_class_counts.csv")

# ===========================================================================
# 2. Retention: how much of each protein-defined class RNA QC kept
# ===========================================================================
retention <- protein_retention(cells, age_levels, class_levels)
write_table(retention, cfg, "protein_gate_retention.csv")
log_step("fraction of each protein class kept by RNA QC:")
print(retention[retention$protein_class %in% stage_levels,
                c("sex", "protein_class", "age", "n_before_qc", "n_passed",
                  "retained_fraction", "median_umi_failed")], row.names = FALSE)

# ===========================================================================
# 3. Composition, before and after QC, and its trend with age
# ===========================================================================
composition <- protein_composition(cells, cfg, stage_levels)
write_table(composition, cfg, "protein_stage_composition.csv")

trends <- do.call(rbind, lapply(split(composition, composition$cell_set), function(comp) {
  t <- composition_age_test(comp, cfg)$trends
  if (!is.null(t)) t$cell_set <- comp$cell_set[1]
  t
}))
if (!is.null(trends)) {
  write_table(trends, cfg, "protein_stage_age_trends.csv")
  log_step("protein-stage share per age step (log-odds), before vs after RNA QC:")
  print(trends[, c("cell_set", "sex", "stage", "prop_first", "prop_last",
                   "log_odds_per_step", "padj")], row.names = FALSE)

  # The verdict on the male trend, said in words, from the numbers above.
  mature <- stage_levels[3]
  for (this_sex in unique(trends$sex)) {
    before <- trends[trends$sex == this_sex & trends$stage == mature &
                       grepl("^all", trends$cell_set), ]
    after <- trends[trends$sex == this_sex & trends$stage == mature &
                      trends$cell_set == "passed RNA QC", ]
    if (!nrow(before) || !nrow(after)) next
    log_step(sprintf(
      "%s %s: %.1f%% -> %.1f%% before QC (slope %+.2f), %.1f%% -> %.1f%% after (slope %+.2f)",
      this_sex, mature, 100 * before$prop_first, 100 * before$prop_last,
      before$log_odds_per_step, 100 * after$prop_first, 100 * after$prop_last,
      after$log_odds_per_step))
  }
}

# ===========================================================================
# 4. Protein class against RNA stage, for the cells that have both
# ===========================================================================
if (file.exists(object_path(cfg, "gmp_neutrophils.rds"))) {
  gmp_neu <- read_object(cfg, "gmp_neutrophils.rds")
  # Seurat's merge renames barcodes that collide between libraries by adding
  # _1 / _2, so cells are matched on sex plus the barcode without that suffix.
  rna_key <- paste(gmp_neu$sex, sub("_[0-9]+$", "", colnames(gmp_neu)))
  cell_key <- paste(cells$sex, cells$barcode)
  rna_stage <- as.character(gmp_neu$stage)[match(cell_key, rna_key)]
  log_step(sum(!is.na(rna_stage)), " protein-gated cells matched to an RNA stage")
  cross <- protein_vs_rna(cells, rna_stage)
  if (!is.null(cross)) write_table(cross, cfg, "protein_vs_rna_stage.csv")
  rm(gmp_neu)
} else {
  log_step("gmp_neutrophils.rds not found; skipping the protein vs RNA comparison")
}

# ===========================================================================
# 5. Figures
# ===========================================================================
save_figure(plot_protein_gate(cells, thresholds), cfg, "protein_gate.pdf",
            width = 9, height = 9)
save_figure(plot_protein_composition(composition, stage_levels), cfg,
            "protein_stage_composition.pdf", width = 9, height = 7)
save_figure(plot_protein_retention(retention, stage_levels), cfg,
            "protein_gate_retention.pdf", width = 10, height = 5)

log_step("step 12 complete")
