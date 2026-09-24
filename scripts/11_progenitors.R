#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 11 - the progenitor compartment upstream of GMP
#
#   Rscript scripts/11_progenitors.R
#
# Steps 3-10 analyse granulopoiesis from the GMP onward. That begins one step
# too late to answer where an age-related shift originates: if aged marrow
# holds proportionally more GMPs, the question is whether GMPs are being made
# faster, or whether the HSC and MPP pool upstream has itself shifted toward
# the myeloid branch.
#
# The cells are already labelled -- SingleR against ImmGen resolves LT-HSC,
# ST-HSC, MPP, CMP, GMP, MEP, CLP, MDP and CDP, renamed by IMMGEN_LABEL_MAP --
# so this step needs no new annotation, only a different subset of the merged
# object and the same imbalance machinery step 9 applies to stages.
#
# The sibling branches are kept deliberately. Composition is a closed system:
# whether GMP output rises at the expense of the erythroid or the lymphoid
# branch is part of the answer, and restricting to the myeloid axis would
# hide it.
#
# Reads:  results/objects/bm_merged.rds
# Writes: results/tables/progenitor_composition.csv
#         results/tables/progenitor_imbalance.csv
#         results/tables/progenitor_composition_clr.csv
#         results/tables/endpoint_progenitor_<population>_<sex>.csv
#         results/figures/progenitor_composition.pdf, progenitor_trends.pdf
# ---------------------------------------------------------------------------

source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()

levels_wanted <- unlist(cfg$analysis$progenitor_levels)
if (!length(levels_wanted))
  stop("analysis.progenitor_levels is empty; nothing to analyse.")

bm <- read_object(cfg, "bm_merged.rds")
require_metadata(bm, c("fine_label_readable", "age", "sex"), context = "step 11")

label <- as.character(bm$fine_label_readable)
present <- intersect(levels_wanted, unique(label))
absent <- setdiff(levels_wanted, present)
if (length(absent))
  log_step("progenitor labels configured but absent from the data: ",
           paste(absent, collapse = ", "))
if (!length(present))
  stop("none of analysis.progenitor_levels appear in fine_label_readable.\n",
       "The object has: ", paste(utils::head(sort(unique(label)), 30),
                                 collapse = ", "))

# The reverse check. A label the data carries but the config does not list is
# silently outside the analysis, and a typo in progenitor_levels looks exactly
# like a population that was never there. Only progenitor-like labels are
# reported, since the object also holds every mature lineage.
unlisted <- setdiff(grep("HSC|Progenitor|MEP|CMP|GMP|CLP|MDP|CDP|Stem",
                         unique(label), value = TRUE), levels_wanted)
if (length(unlisted))
  log_step("progenitor-like labels in the data but NOT in ",
           "analysis.progenitor_levels (excluded from this analysis): ",
           paste(sort(unlisted), collapse = ", "))

prog <- select_cells(bm, list(
  "fine_label_readable is a configured progenitor" = label %in% present,
  "age is one of analysis.age_levels" = bm$age %in% cfg$analysis$age_levels
), context = "progenitor compartment")
prog$progenitor <- factor(as.character(prog$fine_label_readable), levels = present)

log_step("progenitor populations by age and sex:")
print(table(prog$progenitor, prog$age, prog$sex))

# ===========================================================================
# 1. Composition
# ===========================================================================
# stage_composition() names its output column `stage` and the plotting and CLR
# helpers read it, so the object keeps that name and only the written copy is
# relabelled. Renaming in place and back again would leave the plots one edit
# away from breaking silently.
composition <- stage_composition(prog, cfg, stage_col = "progenitor",
                                 levels = present)
written <- composition
names(written)[names(written) == "stage"] <- "population"
write_table(written, cfg, "progenitor_composition.csv")

save_figure(plot_stage_composition(composition), cfg,
            "progenitor_composition.pdf", width = 10, height = 5)
save_figure(plot_stage_trends(composition), cfg,
            "progenitor_trends.pdf", width = 11, height = 7)

clr <- composition_clr_matrix(composition)
write_table(as.data.frame(clr), cfg, "progenitor_composition_clr.csv",
            row.names = TRUE)

# ===========================================================================
# 2. Imbalance
# ===========================================================================
counts <- population_counts(prog, cfg, population_col = "progenitor",
                            levels = present)

# A log ratio built on a handful of cells is a statement about a handful of
# cells. Populations that never reach the floor in any group of a contrast are
# dropped before scoring rather than reported with a wide interval, which
# reads as a measured null when it is an absence of data.
floor_n <- cfg$imbalance$min_cells_per_group %||% 20
per_population <- stats::aggregate(n ~ population, counts, max)
keep <- as.character(per_population$population[per_population$n >= floor_n])
dropped <- setdiff(as.character(unique(counts$population)), keep)
if (length(dropped))
  log_step("below ", floor_n, " cells in every group, dropped from the ",
           "imbalance scoring: ", paste(dropped, collapse = ", "))

scored <- counts[as.character(counts$population) %in% keep, ]
scored$population <- factor(as.character(scored$population),
                            levels = intersect(present, keep))

imbalance <- imbalance_report(scored, cfg)
if (!is.null(imbalance)) {
  write_table(imbalance, cfg, "progenitor_imbalance.csv")
  report_imbalance(imbalance)
}

# ===========================================================================
# 3. Expression change within each progenitor population
# ===========================================================================
# The same 3m vs 18m contrast step 10 runs on maturation stages. A shift in
# how many GMPs there are and a shift in what GMPs express are different
# claims, and the composition tables above cannot separate them.
# The matched assay built in step 3 covers the GMP/neutrophil subset only, so
# it is not on this object. Matching is rebuilt here, over the progenitor
# compartment's own cells: thinning toward a target computed from a different
# population would be the wrong target for this one.
assay <- "RNA"
if (isTRUE(cfg$depth$match)) {
  match_on <- cfg$depth$match_on %||% "sex"
  require_metadata(prog, match_on, context = "progenitor depth matching")
  prog <- add_matched_assay(prog, cfg, group_col = match_on)
  assay <- "RNAmatched"
}
log_step("progenitor expression contrasts read the ", assay, " assay")

summary_rows <- list()

for (population in present) {
  for (sex in cfg$analysis$sex_levels) {
    label_txt <- paste(population, sex, sep = "/")
    cells <- colnames(prog)[!is.na(prog$progenitor) &
                              as.character(prog$progenitor) == population &
                              as.character(prog$sex) == sex]
    if (length(cells) < cfg$gates$min_cells_per_stratum) {
      log_step("  skipping ", label_txt, ": ", length(cells), " cells in total")
      next
    }
    log_step("=== ", label_txt, " ===")
    sub <- subset(drop_graphs(prog), cells = cells)

    res <- endpoint_contrast(sub, cfg, label = label_txt, assay = assay)
    if (is.null(res)) next

    n_hits <- if (is.null(res$genes)) 0L else nrow(res$genes)
    if (n_hits) {
      out <- res$genes
      out$population <- population
      out$sex <- sex
      write_table(strip_inferential_columns(out), cfg,
                  sprintf("endpoint_progenitor_%s_%s.csv",
                          gsub("[^A-Za-z0-9]+", "_", population), sex))
    }
    summary_rows[[label_txt]] <- data.frame(
      population = population, sex = sex,
      n_reference = as.integer(res$n_by_age[cfg$endpoint$reference] %||% 0L),
      n_endpoint = as.integer(res$n_by_age[cfg$endpoint$endpoint] %||% 0L),
      fdr_threshold = res$fdr_threshold,
      fwer_threshold = res$fwer_threshold,
      n_genes = n_hits,
      n_monotonic = if (n_hits) sum(res$genes$shape == "monotonic") else 0L,
      row.names = NULL)
  }
}

if (length(summary_rows)) {
  summary_tbl <- do.call(rbind, summary_rows)
  write_table(summary_tbl, cfg, "endpoint_progenitor_summary.csv")
  log_step(cfg$endpoint$reference, " vs ", cfg$endpoint$endpoint,
           " within each progenitor population:")
  print(summary_tbl, row.names = FALSE)
} else {
  log_step("no progenitor population had enough cells at both endpoints")
}

log_step("step 11 complete")
