#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 7 - negative-binomial GLM models with permutation-calibrated FDR.
#
#   Rscript scripts/07_glm_models.R
#
# Part A, cell-level omnibus test in differentiated neutrophils:
#   - per sex:  ~ age      vs ~ 1          ("does this gene change with age")
#   - joint:    ~ age*sex  vs ~ age + sex  ("does it change differently")
#
# Part B, the same question along the trajectory: does a gene's developmental
# profile differ by age, by sex, or in shape.
#
# BH adjustment is reported but not used for selection: the cutoff comes from
# refitting each model on shuffled labels, which is the only null that accounts
# for cells being pseudoreplicates of pooled animals.
#
# Reads:  results/objects/gmp_neutrophils.rds, combined_cds.rds
# Writes: results/tables/omnibus_<sex>.csv, interaction_age_sex.csv,
#         pt_age_sex_hits.csv, pt_age_sex_fitted.csv
#         results/objects/pt_age_sex_fits.rds
#
# This step refits every model NPERM times. It is by far the slowest part of
# the pipeline; lower glm_de.n_perm in the config for a quick pass.
# ---------------------------------------------------------------------------

# Run from the project root ("Rscript scripts/07_glm_models.R") or from inside
# scripts/; setup.R locates the project root from either.
source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()
require_packages("glmGamPoi", "monocle3", "splines", "SingleCellExperiment")

age_levels <- cfg$analysis$age_levels
target_fdr <- cfg$glm_de$target_fdr
n_perm <- cfg$glm_de$n_perm

# ===========================================================================
# Part A - cell-level omnibus test
# ===========================================================================
gmp_neu <- read_object(cfg, "gmp_neutrophils.rds")

keep <- !is.na(gmp_neu$CytoTRACE2_Potency) &
  gmp_neu$CytoTRACE2_Potency == "Differentiated" &
  gmp_neu$age %in% age_levels
neus <- subset(gmp_neu, cells = colnames(gmp_neu)[keep])
neus <- Seurat::JoinLayers(neus)
Seurat::DefaultAssay(neus) <- "RNA"

meta <- neus@meta.data[, c("age", "sex")]
meta$age <- factor(as.character(meta$age), levels = age_levels)
meta$sex <- factor(as.character(meta$sex), levels = cfg$analysis$sex_levels)

counts <- Seurat::GetAssayData(neus, assay = "RNA", layer = "counts")
ok <- !is.na(meta$age) & !is.na(meta$sex)
counts <- counts[, ok]
meta <- droplevels(meta[ok, ])
print(table(meta$sex, meta$age))

# Detection filter defined once, on all cells, so every model sees the same
# genes and the permuted fits stay comparable.
counts <- counts[Matrix::rowSums(counts > 0) >= cfg$glm_de$min_cells_detected, ]
log_step(sprintf("%d genes x %d cells", nrow(counts), ncol(counts)))

is_female <- meta$sex == "female"
is_male <- meta$sex == "male"

log_step("fitting observed models")
fits <- list(
  female = glm_lrt(counts[, is_female], droplevels(meta[is_female, ]), ~ age, ~ 1),
  male   = glm_lrt(counts[, is_male],   droplevels(meta[is_male, ]),   ~ age, ~ 1),
  interaction = glm_lrt(counts, meta, ~ age * sex, ~ age + sex)
)
# All three spend the same degrees of freedom, which is what makes their
# thresholds comparable.
expected_df <- length(age_levels) - 1
stopifnot(vapply(fits, function(f) all(f$res$df1 == expected_df), logical(1)))

log_step("permuting (", n_perm, " refits per model)")
nulls <- list(
  female = permutation_null(counts[, is_female], droplevels(meta[is_female, ]),
                            ~ age, ~ 1, shuffle_age, n_perm),
  male   = permutation_null(counts[, is_male], droplevels(meta[is_male, ]),
                            ~ age, ~ 1, shuffle_age, n_perm),
  interaction = permutation_null(counts, meta, ~ age * sex, ~ age + sex,
                                 shuffle_sex_within_age, n_perm)
)

calibrations <- Map(function(f, n) calibrate_threshold(f$res$pval, n, target_fdr),
                    fits, nulls)

z_scores <- lapply(fits[c("female", "male")], function(f)
  t(scale(t(glm_group_means(f$fit, age_levels)))))

tables <- list(
  female      = glm_hits_table(fits$female, calibrations$female, z_scores$female),
  male        = glm_hits_table(fits$male, calibrations$male, z_scores$male),
  interaction = glm_hits_table(fits$interaction, calibrations$interaction)
)

# Sensitivity: the middle timepoint is the smallest group in one arm, so check
# how much of the interaction signal survives without it.
middle <- age_levels[2]
sub <- meta$age != middle
no_middle <- glm_lrt(counts[, sub], droplevels(meta[sub, ]), ~ age * sex, ~ age + sex)
robust_interaction <- intersect(tables$interaction$name,
                                no_middle$res$name[no_middle$res$pval < 0.05])

for (nm in names(tables))
  write_table(tables[[nm]], cfg,
              if (nm == "interaction") "interaction_age_sex.csv"
              else paste0("omnibus_", nm, ".csv"))

for (nm in c("female", "male"))
  log_step(sprintf("%-6s: thresh %.2e -> %d genes (BH<0.05: %d)", nm,
                   calibrations[[nm]]$threshold, nrow(tables[[nm]]),
                   sum(fits[[nm]]$res$adj_pval < 0.05)))
log_step(sprintf("age*sex: thresh %.2e -> %d genes; %d survive dropping %s",
                 calibrations$interaction$threshold, nrow(tables$interaction),
                 length(robust_interaction), middle))

# ===========================================================================
# Part B - the same question along the trajectory
# ===========================================================================
combined_cds <- read_object(cfg, "combined_cds.rds")

cd <- SummarizedExperiment::colData(combined_cds)
cd$pt <- monocle3::pseudotime(combined_cds)
cd$age <- factor(as.character(cd$age), levels = age_levels)
cd$sex <- factor(as.character(cd$sex), levels = cfg$analysis$sex_levels)

# The middle timepoint is dropped here: the model needs every group populated
# across the whole pseudotime span, and it is the sparsest.
keep_cells <- is.finite(cd$pt) & !is.na(cd$age) & !is.na(cd$sex) & cd$age != middle
combined_cds <- combined_cds[, keep_cells]
cd <- droplevels(cd[keep_cells, ])
cd$grp <- droplevels(interaction(cd$age, cd$sex, sep = "_"))

# Group-specific splines extrapolate freely outside their own cell range, and
# that divergence reads as a real difference. Trim to the span all groups cover.
ranges <- tapply(cd$pt, cd$grp, range)
lo <- max(vapply(ranges, `[`, numeric(1), 1))
hi <- min(vapply(ranges, `[`, numeric(1), 2))
log_step(sprintf("shared pseudotime span: %.3f - %.3f (full: %.3f - %.3f)",
                 lo, hi, min(cd$pt), max(cd$pt)))

quintiles <- stats::quantile(cd$pt, 0:5 / 5)
print(table(cut(cd$pt, quintiles, include.lowest = TRUE), cd$grp))

in_span <- cd$pt >= lo & cd$pt <= hi
log_step(sprintf("trimming removes %d of %d cells (%.1f%%)",
                 sum(!in_span), length(in_span), 100 * mean(!in_span)))
combined_cds <- combined_cds[, in_span]
cd <- droplevels(cd[in_span, ])
print(table(cd$sex, cd$age))

cts <- SingleCellExperiment::counts(combined_cds)
cts <- cts[Matrix::rowSums(cts > 0) >= cfg$glm_de$min_cells_detected, ]
log_step(sprintf("%d genes x %d cells", nrow(cts), ncol(cts)))

PT_DF <- cfg$glm_de$pt_spline_df

# LEAN is the primary model: it lets both the age effect and the sex effect
# vary along pseudotime, but omits the three-way term, which costs 3 df and is
# the worst-estimated piece. FULL adds it, to see whether it earns them.
f_lean   <- ~ splines::ns(pt, PT_DF) + age * sex +
             splines::ns(pt, PT_DF):age + splines::ns(pt, PT_DF):sex
f_full   <- ~ splines::ns(pt, PT_DF) * age * sex
f_null   <- ~ splines::ns(pt, PT_DF)          # developmental profile only
f_no_age <- ~ splines::ns(pt, PT_DF) * sex    # -> any age effect
f_no_sex <- ~ splines::ns(pt, PT_DF) * age    # -> any sex effect

stopifnot(!anyNA(stats::model.frame(f_lean, as.data.frame(cd), na.action = stats::na.pass)))

traj <- list(
  any  = glm_lrt(cts, cd, f_lean, f_null),
  full = glm_lrt(cts, cd, f_full, f_null),
  age  = glm_lrt(cts, cd, f_full, f_no_age),
  sex  = glm_lrt(cts, cd, f_full, f_no_sex)
)
log_step("df spent -- ", paste(sprintf("%s: %d", names(traj),
         vapply(traj, function(t) t$res$df1[1], numeric(1))), collapse = "  "))

# Null: no dependence on age or sex, so the two are shuffled together and
# pseudotime stays fixed. That keeps each cell's developmental position and
# library depth while destroying the group labels.
traj_nulls <- list(
  any = permutation_null(cts, cd, f_lean, f_null, shuffle_age_sex, n_perm),
  age = permutation_null(cts, cd, f_full, f_no_age, shuffle_age_sex, n_perm),
  sex = permutation_null(cts, cd, f_full, f_no_sex, shuffle_age_sex, n_perm)
)
traj_cal <- Map(function(t, n) calibrate_threshold(t$res$pval, n, target_fdr),
                traj[names(traj_nulls)], traj_nulls)

hits <- Map(function(t, cal) t$res$name[t$res$pval < cal$threshold],
            traj[names(traj_cal)], traj_cal)

tab <- traj$any$res[traj$any$res$name %in% hits$any,
                    c("name", "pval", "adj_pval", "f_statistic")]
tab <- tab[order(tab$pval), ]
tab$age_effect <- tab$name %in% hits$age
tab$sex_effect <- tab$name %in% hits$sex
tab$driver <- with(tab, ifelse(age_effect & sex_effect, "both",
                        ifelse(age_effect, "age",
                        ifelse(sex_effect, "sex", "joint_only"))))

# Where along the trajectory each hit differs: fitted values per group at five
# points across the span, for characterisation rather than testing.
#
# The prediction grid must reuse the spline basis the model was fitted with.
# Calling model.matrix() on a fresh grid re-derives the knots from that grid,
# which silently produces a different basis and therefore wrong fitted values.
# model.frame() records the fitted basis in the terms "predvars" attribute,
# and model.matrix() reuses it when handed those terms.
pt_grid <- seq(min(cd$pt), max(cd$pt), length.out = 5)
newdat <- expand.grid(pt = pt_grid, age = levels(cd$age), sex = levels(cd$sex))
fitted_terms <- stats::terms(stats::model.frame(f_lean, as.data.frame(cd)))
mm_new <- stats::model.matrix(fitted_terms, newdat)
stopifnot(identical(colnames(mm_new), colnames(traj$any$fit$Beta)))

pred <- traj$any$fit$Beta[tab$name, , drop = FALSE] %*% t(mm_new)
colnames(pred) <- with(newdat, paste0("pt", round(pt, 2), "_", age, "_", sex))

log_step(sprintf("any age/sex effect: thresh %.2e -> %d genes (BH<0.05: %d)",
                 traj_cal$any$threshold, length(hits$any),
                 sum(traj$any$res$adj_pval < 0.05)))
log_step(sprintf("  age-driven %d | sex-driven %d | both %d | joint only %d",
                 sum(tab$driver == "age"), sum(tab$driver == "sex"),
                 sum(tab$driver == "both"), sum(tab$driver == "joint_only")))
log_step(sprintf("three-way model returns %d genes at BH<0.05 vs %d for lean",
                 sum(traj$full$res$adj_pval < 0.05), sum(traj$any$res$adj_pval < 0.05)))

log_step("overlap with the cell-level female hits: ",
         length(intersect(hits$any, tables$female$name)),
         " | trajectory-only: ", length(setdiff(hits$any, tables$female$name)))

write_table(tab, cfg, "pt_age_sex_hits.csv")
write_table(as.data.frame(pred), cfg, "pt_age_sex_fitted.csv", row.names = TRUE)
save_object(list(fits = traj, calibrations = traj_cal, col_data = cd),
            cfg, "pt_age_sex_fits.rds")

log_step("step 7 complete")
