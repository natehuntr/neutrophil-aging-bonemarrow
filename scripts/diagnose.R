#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Inspect the saved objects: what is present, and what any later step is
# about to fail on.
#
#   Rscript scripts/diagnose.R
#
# Reads whatever exists under results/objects and reports nothing but facts:
# cell counts, the metadata columns each step needs, and how many cells carry
# a value for each. Run it when a step reports an empty subset, or before a
# long run, to find out which upstream step did not write what.
# ---------------------------------------------------------------------------

# Run from the project root ("Rscript scripts/diagnose.R") or from inside
# scripts/; setup.R locates the project root from either.
source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project(quiet = TRUE)
load_modules()

# Columns each step selects on, and so cannot run without.
REQUIREMENTS <- list(
  "<sex>_filtered.rds"   = c("age", "sex", "MULTI_ID", "scDblFinder.class",
                             "log10GenesPerUMI", "mitoRatio", "percent.hb"),
  "<sex>_annotated.rds"  = c("age", "sex", "singleR_main_label",
                             "singleR_fine_label", "fine_label_readable"),
  "<sex>_neustem.rds"    = c("age", "sex", "singleR_fine_label",
                             "fine_label_readable", "clusters",
                             "CytoTRACE2_Potency", "CytoTRACE2_Score"),
  "bm_merged.rds"        = c("age", "sex", "singleR_main_label",
                             "CytoTRACE2_Potency"),
  "gmp_neutrophils.rds"  = c("age", "sex", "fine_neu_labels",
                             "CytoTRACE2_Potency", "seurat_clusters")
)

describe <- function(name, required) {
  path <- object_path(cfg, name)
  if (!file.exists(path)) {
    cat("\n", name, ": not present\n", sep = "")
    return(invisible(NULL))
  }

  obj <- readRDS(path)
  cat("\n=== ", name, " ===\n", sep = "")
  cat(ncol(obj), "cells x", nrow(obj), "features | assays:",
      paste(assay_names(obj), collapse = ", "), "\n")
  if (length(obj@reductions))
    cat("reductions:", paste(names(obj@reductions), collapse = ", "), "\n")

  coverage <- metadata_coverage(obj, required)
  coverage$pct <- round(100 * coverage$non_na / ncol(obj), 1)
  print(coverage, row.names = FALSE)

  broken <- coverage[!coverage$present | coverage$non_na == 0, ]
  if (nrow(broken))
    cat("!! empty or absent: ", paste(broken$column, collapse = ", "),
        "\n   Any step selecting on these will find no cells.\n", sep = "")

  for (col in c("age", "sex", "fine_neu_labels", "CytoTRACE2_Potency",
                "singleR_main_label")) {
    if (!col %in% colnames(obj@meta.data)) next
    values <- obj@meta.data[[col]]
    if (length(unique(values)) > 15) next
    cat("\n", col, ":\n", sep = "")
    print(table(values, useNA = "ifany"))
  }
  invisible(NULL)
}

for (key in names(cfg$samples)) {
  for (template in names(REQUIREMENTS)[1:3])
    describe(sub("<sex>", key, template), REQUIREMENTS[[template]])
}
for (name in names(REQUIREMENTS)[4:5])
  describe(name, REQUIREMENTS[[name]])

cat("\nDone. A column marked absent or with non_na = 0 is what an empty",
    "\nsubset downstream is complaining about; re-run the step that writes it.\n")
