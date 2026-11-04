#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 13 - summary figures from the results tables.
#
#   Rscript scripts/13_figures.R                    # reads results/tables
#   Rscript scripts/13_figures.R results/summary    # or any folder of tables
#
# Draws one figure per finding, from tables only: no Seurat object, no
# analysis package, nothing but ggplot2, patchwork and scales. It runs in
# seconds on a laptop, can be pointed at a results bundle, and redraws every
# figure after any step re-runs.
#
# A figure whose table is missing is skipped and listed as such in the index
# -- the index is the checklist of what the pipeline has and has not produced.
#
# Writes: results/figures/summary/figNN_*.pdf (+ .png)
#         results/figures/summary/FIGURES.md
# ---------------------------------------------------------------------------

source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
# Only the config and the plotting code: this step must run without Seurat.
cfg <- load_config()
for (pkg in c("ggplot2", "patchwork", "scales"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("step 13 needs ", pkg, "; install.packages('", pkg, "')")
for (m in c("plots", "figures", "protein_gate", "composition"))
  source(project_path("R", paste0(m, ".R")))

args <- commandArgs(trailingOnly = TRUE)
table_dir <- if (length(args)) normalizePath(args[1]) else cfg$paths$tables
out_dir <- file.path(cfg$paths$figures, "summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cfg$paths$figures <- out_dir
log_step("reading tables from ", table_dir)

read_tbl <- function(name) {
  path <- file.path(table_dir, name)
  if (!file.exists(path)) return(NULL)
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(df)) NULL else df
}

#' Read every table matching a pattern and bind them, adding columns parsed
#' from the file name where the table does not carry them.
read_many <- function(pattern, parse = NULL, exclude = NULL) {
  files <- list.files(table_dir, pattern)
  if (!is.null(exclude)) files <- files[!grepl(exclude, files)]
  tabs <- lapply(files, function(f) {
    df <- read_tbl(f)
    if (is.null(df)) return(NULL)
    if (!is.null(parse)) {
      extra <- parse(f)
      for (nm in names(extra)) if (!nm %in% names(df)) df[[nm]] <- extra[[nm]]
    }
    df
  })
  tabs <- Filter(Negate(is.null), tabs)
  if (!length(tabs)) return(NULL)
  common <- Reduce(intersect, lapply(tabs, names))
  do.call(rbind, lapply(tabs, `[`, common))
}

index <- list()
figure <- function(id, file, what, source, build, width = 10, height = 6) {
  plot <- tryCatch(build(), error = function(e) {
    log_step("  ", id, ": could not draw -- ", conditionMessage(e))
    structure(list(), class = "draw_error", message = conditionMessage(e))
  })
  status <- if (inherits(plot, "draw_error")) paste("NOT DRAWN:", attr(plot, "message"))
            else if (is.null(plot)) "NOT DRAWN: table missing, empty, or nothing passes its filter"
            else {
              save_figure(plot, cfg, paste0(id, "_", file, ".pdf"), width = width, height = height)
              "drawn"
            }
  index[[length(index) + 1]] <<- data.frame(figure = paste0(id, "_", file), shows = what,
                                            from = source, status = status,
                                            stringsAsFactors = FALSE)
}
need <- function(...) { x <- list(...); if (any(vapply(x, is.null, logical(1)))) NULL else TRUE }

# --- what the data can be trusted for --------------------------------------
depth <- read_tbl("depth_by_age_within_sex.csv")
figure("fig01", "depth_by_stage",
       "Median UMIs per stage, age and library -- which comparisons are depth-controlled",
       "depth_by_age_within_sex.csv",
       function() if (is.null(depth)) NULL else fig_depth_by_stage(depth), 12, 4.5)

targets <- read_tbl("depth_matching_targets.csv")
figure("fig02", "depth_matching",
       "Depth ratio inside each matching stratum, before and after thinning",
       "depth_matching_targets.csv",
       function() if (is.null(targets)) NULL else fig_matching_targets(targets), 10, 6)

confusion <- read_tbl("stage_assignment_confusion.csv")
figure("fig03", "stage_rna_vs_protein",
       "RNA-signature stage against surface-protein stage",
       "stage_assignment_confusion.csv",
       function() if (is.null(confusion)) NULL else fig_stage_confusion(confusion), 8, 6)

# --- composition -----------------------------------------------------------
comp <- read_tbl("stage_composition.csv")
trends <- read_tbl("stage_age_trends.csv")
figure("fig04", "stage_composition",
       "Maturation stage mix by age, per library",
       "stage_composition.csv",
       function() if (is.null(comp)) NULL else fig_stage_composition(comp), 9, 5.5)
figure("fig05", "stage_trends",
       "Each stage's share across age, with the fitted per-step slope",
       "stage_composition.csv, stage_age_trends.csv",
       function() if (is.null(comp)) NULL else fig_stage_trends(comp, trends), 13, 4.5)

stage_imb <- read_tbl("stage_imbalance.csv")
figure("fig06", "stage_imbalance",
       "Which stages are over- or under-represented at each age",
       "stage_imbalance.csv",
       function() if (is.null(stage_imb)) NULL else
         fig_imbalance(stage_imb, "Maturation stages over- and under-represented, by age",
                       order = STAGE_ORDER), 12, 4.5)

prog_imb <- read_tbl("progenitor_imbalance.csv")
figure("fig07", "progenitor_imbalance",
       "Progenitor populations over- and under-represented by age (HSC expansion, MDP loss)",
       "progenitor_imbalance.csv",
       function() if (is.null(prog_imb)) NULL else
         fig_imbalance(prog_imb, "Progenitor compartment: over- and under-representation by age",
                       short_progenitor), 13, 6)

prog_comp <- read_tbl("progenitor_composition.csv")
figure("fig08", "progenitor_trends",
       "Progenitor shares across age",
       "progenitor_composition.csv",
       function() if (is.null(prog_comp)) NULL else fig_progenitor_trends(prog_comp), 12, 8)

# --- protein gate ----------------------------------------------------------
p_comp <- read_tbl("protein_stage_composition.csv")
p_ret <- read_tbl("protein_gate_retention.csv")
p_levels <- protein_stage_levels(!is.null(p_comp) && any(grepl("CXCR2", p_comp$stage)))
figure("fig09", "protein_stage_composition",
       "Stage mix from surface protein, before and after RNA QC",
       "protein_stage_composition.csv",
       function() if (is.null(p_comp)) NULL else plot_protein_composition(p_comp, p_levels), 9, 7)
figure("fig10", "protein_qc_retention",
       "Share of each protein-defined stage that RNA QC kept, by age -- settles the male mature-cell trend",
       "protein_gate_retention.csv",
       function() if (is.null(p_ret)) NULL else plot_protein_retention(p_ret, p_levels), 11, 4.8)

# --- position along the trajectory -----------------------------------------
pt <- read_tbl("pseudotime_shifts.csv")
pt_w <- read_tbl("pseudotime_shifts_within_stage.csv")
figure("fig11", "pseudotime_shifts",
       "Where cells sit along the trajectory relative to 3m, pooled and within stage",
       "pseudotime_shifts.csv, pseudotime_shifts_within_stage.csv",
       function() if (is.null(pt)) NULL else fig_distribution_shifts(pt, pt_w, "Pseudotime"), 13, 4)
pot <- read_tbl("potency_shifts.csv")
pot_w <- read_tbl("potency_shifts_within_stage.csv")
figure("fig12", "potency_shifts",
       "CytoTRACE2 potency relative to 3m, pooled and within stage",
       "potency_shifts.csv, potency_shifts_within_stage.csv",
       function() if (is.null(pot)) NULL else fig_distribution_shifts(pot, pot_w, "Potency"), 13, 4)

# --- expression --------------------------------------------------------------
ep_sum <- read_tbl("endpoint_summary.csv")
ep_prog <- read_tbl("endpoint_progenitor_summary.csv")
figure("fig13", "endpoint_gene_counts",
       "Genes changing 3m -> 18m per stratum, with the detection check",
       "endpoint_summary.csv, endpoint_progenitor_summary.csv",
       function() if (is.null(ep_sum)) NULL else fig_endpoint_summary(ep_sum, ep_prog), 9, 7)

ep_genes <- read_many("^endpoint_(GMPs|proNeu|preNeu|immature|mature)_(female|male)\\.csv$")
figure("fig14", "endpoint_gene_profiles",
       "The 3m vs 18m genes followed through all four ages (Lcn2, Retnlg, Fcer1a, ...)",
       "endpoint_<stage>_<sex>.csv",
       function() if (is.null(ep_genes)) NULL else fig_endpoint_profiles(ep_genes), 12, 8)

ep_gsea <- read_many("^endpoint_gsea_(GMPs|proNeu|preNeu|immature|mature)_(female|male)\\.csv$")
figure("fig15", "endpoint_gsea",
       "Gene sets moving 3m -> 18m per stratum (OXPHOS, MHC-I, interferon, ribosome biogenesis)",
       "endpoint_gsea_<stage>_<sex>.csv",
       function() if (is.null(ep_gsea)) NULL else fig_endpoint_gsea(ep_gsea), 13, 10)

excess <- read_tbl("age_trend_excess_over_null.csv")
figure("fig16", "age_trends_over_null",
       "Monotonic age trends beating the permutation null, per stratum",
       "age_trend_excess_over_null.csv",
       function() if (is.null(excess)) NULL else fig_age_trend_excess(excess), 8, 4.5)

parse_changing <- function(f) {
  m <- regmatches(f, regexec("^age_changing_([^_]+)_(female|male)_", f))[[1]]
  list(stage = m[2], sex = m[3])
}
adt <- read_many("^age_changing_.*_(female|male)_ADT\\.csv$", parse_changing)
figure("fig17", "adt_age_changes",
       "Surface proteins changing with age per stage and sex (FceRIa, Ly6G, ...)",
       "age_changing_<stage>_<sex>_ADT.csv",
       function() if (is.null(adt)) NULL else fig_adt_age_changes(adt), 10, 8)

for (sex in cfg$analysis$sex_levels) {
  omni <- read_tbl(paste0("omnibus_", sex, ".csv"))
  figure(if (sex == "female") "fig18" else "fig19", paste0("omnibus_", sex),
         paste("Genes changing with age in", sex, "differentiated cells (NB-GLM, stage-adjusted)"),
         paste0("omnibus_", sex, ".csv"),
         function() if (is.null(omni)) NULL else fig_omnibus_heatmap(omni, sex), 6, 9)
}

pt_hits <- read_tbl("pt_age_sex_hits.csv")
pt_fit <- read_tbl("pt_age_sex_fitted.csv")
figure("fig20", "trajectory_model_fits",
       "Developmental profiles that change with age or sex (step 7B)",
       "pt_age_sex_hits.csv, pt_age_sex_fitted.csv",
       function() if (is.null(need(pt_hits, pt_fit))) NULL else fig_trajectory_fits(pt_hits, pt_fit), 13, 9)

int_gsea <- read_many("^gsea_sex_by_age_GOBP_.*\\.csv$",
                      function(f) list(stage = sub("^gsea_sex_by_age_GOBP_(.*)\\.csv$", "\\1", f)))
figure("fig21", "interaction_gsea",
       "Pathways whose age trend differs between the two libraries, per stage",
       "gsea_sex_by_age_GOBP_<stage>.csv",
       function() if (is.null(int_gsea)) NULL else fig_interaction_gsea(int_gsea), 14, 9)

# --- the sex contrast --------------------------------------------------------
cand <- read_tbl("sex_candidates.csv")
figure("fig22", "sex_candidates",
       "Sex-difference candidates that clear the null and depth checks, with intervals",
       "sex_candidates.csv",
       function() if (is.null(cand)) NULL else fig_sex_candidates(cand), 11, 8)
null <- read_tbl("within_library_null.csv")
figure("fig23", "within_library_null",
       "Male-vs-female separation per stage against random splits of one library",
       "within_library_null.csv",
       function() if (is.null(null)) NULL else fig_within_library_null(null), 8, 4)

# --- index -------------------------------------------------------------------
idx <- do.call(rbind, index)
lines <- c("# Summary figures", "",
           paste0("Drawn ", format(Sys.time(), "%Y-%m-%d %H:%M"), " from `", table_dir, "`."), "",
           "Each figure is drawn from the table named beside it, so it can be checked against",
           "the numbers. A figure marked NOT DRAWN names the table it is waiting for.", "",
           "| figure | shows | from | status |", "|---|---|---|---|",
           sprintf("| `%s` | %s | `%s` | %s |", idx$figure, idx$shows, idx$from, idx$status))
writeLines(lines, file.path(out_dir, "FIGURES.md"))
log_step(sum(idx$status == "drawn"), " of ", nrow(idx), " figures drawn; index in ",
         file.path(out_dir, "FIGURES.md"))
not_drawn <- idx[idx$status != "drawn", ]
if (nrow(not_drawn))
  for (i in seq_len(nrow(not_drawn)))
    log_step("  ", not_drawn$figure[i], ": ", not_drawn$status[i])
log_step("step 13 complete")
