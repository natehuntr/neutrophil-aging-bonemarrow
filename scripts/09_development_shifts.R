#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 9 - does development itself shift with age?
#
#   Rscript scripts/09_development_shifts.R
#
# The other steps ask how gene expression changes. This one asks whether the
# developmental process moves, in the two ways that expression analyses cannot
# see:
#
#   A. Composition - does the mix of maturation stages shift with age?
#      A rise in the immature fraction is the classic aged-marrow signature,
#      and it is invisible to any within-stage comparison.
#
#   B. Position - do cells sit further along, or pile up before, a maturation
#      step? That is a shift in the distribution of cells over pseudotime,
#      which the gene-level trajectory models do not test.
#
# Both comparisons are made across ages *within* a sex. That is the
# design-clean contrast: all ages share a library, so an age shift cannot be a
# batch shift. The cross-sex comparison is reported once, flagged, because with
# one 10x run per sex it is also a library contrast.
#
# Reads:  results/objects/gmp_neutrophils.rds (part A)
#         results/objects/combined_cds.rds    (part B, from step 5)
# Writes: results/tables/stage_composition.csv, stage_age_trends.csv,
#         stage_composition_clr.csv, stage_sex_differences.csv,
#         pseudotime_shifts.csv, potency_shifts.csv (+ their summaries)
#         results/figures/stage_composition.pdf, stage_trends.pdf,
#         pseudotime_by_age.pdf, potency_by_age.pdf
# ---------------------------------------------------------------------------

# Run from the project root ("Rscript scripts/09_development_shifts.R") or from
# inside scripts/; setup.R locates the project root from either.
source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()

age_levels <- cfg$analysis$age_levels
gmp_neu <- read_object(cfg, "gmp_neutrophils.rds")

# ===========================================================================
# A. Stage composition
# ===========================================================================
log_step("=== stage composition ===")

composition <- stage_composition(gmp_neu, cfg)
write_table(composition, cfg, "stage_composition.csv")

print(stats::xtabs(n ~ stage + age + sex, composition))

save_figure(plot_stage_composition(composition), cfg,
            "stage_composition.pdf", width = 9, height = 5)
save_figure(plot_stage_trends(composition), cfg,
            "stage_trends.pdf", width = 10, height = 6)

age_test <- composition_age_test(composition, cfg)
if (!is.null(age_test$trends)) {
  write_table(age_test$trends, cfg, "stage_age_trends.csv")

  moved <- age_test$trends[order(-abs(age_test$trends$log_odds_per_step)), ]
  log_step("stages ordered by how far their share moves per timepoint:")
  print(moved[, c("sex", "stage", "prop_first", "prop_last",
                  "log_odds_per_step", "padj")], row.names = FALSE)
}

# CLR removes the sum-to-one constraint, so a stage rising is not confused
# with the others being pushed down.
clr <- composition_clr_matrix(composition)
write_table(as.data.frame(clr), cfg, "stage_composition_clr.csv", row.names = TRUE)

sex_test <- composition_sex_test(composition, cfg)
if (!is.null(sex_test)) {
  write_table(sex_test, cfg, "stage_sex_differences.csv")
  log_step("stage composition by sex, within age ",
           "(CONFOUNDED WITH LIBRARY -- one 10x run per sex):")
  print(sex_test, row.names = FALSE)
}

# ===========================================================================
# B. Position along the trajectory
# ===========================================================================
if (!file.exists(object_path(cfg, "combined_cds.rds"))) {
  log_step("combined_cds.rds not found -- run step 5 for the pseudotime half")
} else {
  log_step("=== position along the trajectory ===")
  require_packages("monocle3", "SummarizedExperiment")

  combined_cds <- read_object(cfg, "combined_cds.rds")
  cd <- SummarizedExperiment::colData(combined_cds)

  meta <- data.frame(age = as.character(cd$age), sex = as.character(cd$sex),
                     stringsAsFactors = FALSE)
  pt <- monocle3::pseudotime(combined_cds)

  pt_result <- compare_across_ages_by_sex(pt, meta, cfg)
  if (!is.null(pt_result$summary)) write_table(pt_result$summary, cfg, "pseudotime_summary.csv")
  if (!is.null(pt_result$shifts)) {
    write_table(pt_result$shifts, cfg, "pseudotime_shifts.csv")
    log_step("pseudotime shift relative to ", age_levels[1], ":")
    print(pt_result$shifts[, c("sex", "group", "n", "median_shift", "wasserstein",
                               "frac_beyond_reference_median", "ks_padj")],
          row.names = FALSE)
  }

  save_figure(plot_distribution_by_age(pt, meta, cfg, xlab = "Pseudotime"),
              cfg, "pseudotime_by_age.pdf", width = 9, height = 8)

  # The same question asked of CytoTRACE2 potency, which is independent of the
  # trajectory fit: if both move the same way, the shift is not an artefact of
  # where the trajectory root landed.
  potency <- SummarizedExperiment::colData(combined_cds)$CytoTRACE2_Score
  if (!is.null(potency) && any(is.finite(potency))) {
    log_step("=== potency (CytoTRACE2) ===")
    potency_result <- compare_across_ages_by_sex(as.numeric(potency), meta, cfg)
    if (!is.null(potency_result$summary))
      write_table(potency_result$summary, cfg, "potency_summary.csv")
    if (!is.null(potency_result$shifts)) {
      write_table(potency_result$shifts, cfg, "potency_shifts.csv")
      print(potency_result$shifts[, c("sex", "group", "median_shift",
                                      "wasserstein", "ks_padj")], row.names = FALSE)
    }
    save_figure(plot_distribution_by_age(as.numeric(potency), meta, cfg,
                                         xlab = "CytoTRACE2 potency score"),
                cfg, "potency_by_age.pdf", width = 9, height = 8)
  } else {
    log_step("no CytoTRACE2_Score in the cds -- skipping the potency comparison")
  }
}

log_step("step 9 complete")
