#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Collect the results into one small, self-describing bundle.
#
#   Rscript scripts/collect_results.R
#
# results/tables holds a few hundred megabytes, most of it per-gene tables
# with 32344 rows that nobody reads whole. This copies the tables that carry
# conclusions verbatim, filters the large ones down to the rows that could
# change a conclusion, and writes a manifest saying what each file is and how
# many rows it lost.
#
# The point is a directory small enough to commit or attach, and complete
# enough that someone reading only the bundle can cross-reference tables
# without going back to the cluster.
#
# Writes: results/summary/*.csv and results/summary/MANIFEST.md
# ---------------------------------------------------------------------------

source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
# Only the config, not init_project(): this reads and writes CSVs with base R,
# and should run on a login node or a laptop without Seurat installed.
cfg <- load_config()

out_dir <- file.path(dirname(cfg$paths$tables), "summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- list()

record <- function(file, description, kept, total, note = "") {
  manifest[[length(manifest) + 1]] <<- data.frame(
    file = file, rows = kept, of = total, description = description,
    note = note, stringsAsFactors = FALSE)
}

read_if <- function(name) {
  path <- table_path(cfg, name)
  if (!file.exists(path)) return(NULL)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

#' Copy a table whole. These are small and every row can matter.
take_whole <- function(name, description) {
  df <- read_if(name)
  if (is.null(df)) { record(name, description, 0, 0, "NOT PRODUCED"); return(invisible()) }
  utils::write.csv(df, file.path(out_dir, name), row.names = FALSE)
  record(name, description, nrow(df), nrow(df))
}

#' A logical column if the table has it, all-TRUE if it does not.
#'
#' NOT `is.null(df$col) | df$col %in% TRUE`. On a table without the column
#' that reads TRUE | logical(0), which is logical(0), and the filter silently
#' returns zero rows -- the same zero-length collapse select_cells() exists to
#' prevent, and it bites here whenever an older table predates a column.
column_or <- function(df, name, default = TRUE) {
  if (name %in% names(df)) df[[name]] %in% TRUE else rep(default, nrow(df))
}

#' Copy the rows that could change a conclusion, and say how many were left.
take_filtered <- function(name, description, keep, note) {
  df <- read_if(name)
  if (is.null(df)) { record(name, description, 0, 0, "NOT PRODUCED"); return(invisible()) }
  rows <- tryCatch(keep(df), error = function(e) rep(TRUE, nrow(df)))
  rows[is.na(rows)] <- FALSE
  utils::write.csv(df[rows, , drop = FALSE], file.path(out_dir, name),
                   row.names = FALSE)
  record(name, description, sum(rows), nrow(df), note)
}

# --- conclusions, copied whole ---------------------------------------------
take_whole("endpoint_summary.csv",
           "3m vs 18m per stage x sex: thresholds, gene counts, detection ratio")
take_whole("endpoint_progenitor_summary.csv",
           "the same, per progenitor population")
take_whole("stage_imbalance.csv",
           "maturation stage over/under-representation, both contrasts")
take_whole("progenitor_imbalance.csv",
           "progenitor over/under-representation, both contrasts")
take_whole("stage_composition.csv", "stage counts and proportions by age and sex")
take_whole("progenitor_composition.csv", "progenitor counts and proportions")
take_whole("stage_composition_clr.csv", "stage composition, CLR-transformed")
take_whole("progenitor_composition_clr.csv", "progenitor composition, CLR")
take_whole("stage_age_trends.csv", "per-stage share change per timepoint")
take_whole("stage_sex_differences.csv", "stage composition by sex (CONFOUNDED)")
take_whole("pseudotime_shifts.csv", "pseudotime shift vs 3m, pooled over stages")
take_whole("pseudotime_shifts_within_stage.csv",
           "pseudotime shift vs 3m, WITHIN stage -- read against the pooled one")
take_whole("potency_shifts.csv", "CytoTRACE2 shift vs 3m, pooled over stages")
take_whole("potency_shifts_within_stage.csv",
           "CytoTRACE2 shift vs 3m, WITHIN stage")
take_whole("depth_by_age_within_sex.csv",
           "depth and complexity per age per sex, matched counts, per stage")
take_whole("age_trend_excess_over_null.csv",
           "four-age trend: genes over the permutation null, with detection_rho")
take_whole("stage_assignment_confusion.csv", "RNA vs ADT stage calls")
take_whole("stage_assignment_margins.csv",
           "whether stage disagreements are fuzzy or contradictory, per stage")
take_whole("sex_control_gates.csv", "gate outcomes for the sex contrast")
take_whole("sex_gene_panel_check.csv", "sex-chromosome genes present in the probe set")

# --- per-gene endpoint tables, whole: they are already short ----------------
# list.files() takes no `perl` argument and its regex engine has no lookahead,
# so the exclusions are applied separately.
endpoint_files <- list.files(cfg$paths$tables, "^endpoint_")
endpoint_files <- endpoint_files[!grepl("^endpoint_(gsea|summary|progenitor_summary)",
                                        endpoint_files)]
for (file in endpoint_files) {
  take_whole(file, "genes changing 3m -> 18m, with their 9m/12m shape")
}

# --- sex candidates: the ones the ranking says to believe --------------------
take_filtered("sex_candidates.csv",
              "sex-difference candidates for orthogonal validation",
              function(df) !column_or(df, "below_null_floor", FALSE) &
                column_or(df, "stratum_depth_ok"),
              "dropped rows below the null floor or from unmatched strata")

# --- GSEA: significant and independent --------------------------------------
for (file in list.files(cfg$paths$tables, "^(endpoint_gsea|gsea_sex_by_age)")) {
  take_filtered(file, "gene set enrichment",
                function(df) df$padj < 0.05 & column_or(df, "independent"),
                "kept padj < 0.05 and independent after collapsing")
}

# --- step 7, if it completed -------------------------------------------------
take_filtered("pt_age_sex_hits.csv", "trajectory GLM hits",
              function(df) rep(TRUE, nrow(df)), "")
for (file in list.files(cfg$paths$tables, "^omnibus_")) {
  take_filtered(file, "cell-level NB-GLM omnibus test",
                function(df) seq_len(nrow(df)) %in%
                  utils::head(order(df$p_val_diagnostic %||% df$pval %||%
                                      seq_len(nrow(df))), 500),
                "top 500 by the diagnostic statistic")
}

# --- manifest ----------------------------------------------------------------
tbl <- do.call(rbind, manifest)
tbl <- tbl[order(tbl$file), ]
utils::write.csv(tbl, file.path(out_dir, "MANIFEST.csv"), row.names = FALSE)

lines <- c(
  "# Results bundle",
  "",
  paste0("Collected ", format(Sys.time(), "%Y-%m-%d %H:%M"), " from ",
         normalizePath(cfg$paths$tables)),
  "",
  "Tables that carry conclusions are copied whole. Large tables are filtered",
  "to the rows that could change a conclusion, and `of` says how many rows the",
  "original held. A row count of 0 with note NOT PRODUCED means that step did",
  "not run or produced nothing -- which is itself informative.",
  "",
  "Read `docs/methods.md` alongside this: several columns (evidence_tier,",
  "below_null_floor, detection_skewed, shape) only mean something with the",
  "method behind them.",
  "",
  "| file | rows | of | description | note |",
  "|---|---|---|---|---|",
  sprintf("| `%s` | %d | %d | %s | %s |", tbl$file, tbl$rows, tbl$of,
          tbl$description, tbl$note))
writeLines(lines, file.path(out_dir, "MANIFEST.md"))

size <- sum(file.info(list.files(out_dir, full.names = TRUE))$size, na.rm = TRUE)
log_step("wrote ", nrow(tbl), " files to ", out_dir,
         " (", round(size / 1024), " KB)")
missing <- tbl$file[tbl$note == "NOT PRODUCED"]
if (length(missing))
  log_step("not produced by this run: ", paste(missing, collapse = ", "))
