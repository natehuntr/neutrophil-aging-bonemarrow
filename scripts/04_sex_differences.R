#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Step 4 - male vs female differences within each age, in differentiated
# neutrophils, on both RNA and ADT, followed by GO over-representation.
#
#   Rscript scripts/04_sex_differences.R
#
# Reads:  results/objects/bm_merged.rds
# Writes: results/tables/degs_<age>_<assay>.csv
#         results/tables/go_ora_<age>_<direction>.csv
#         results/objects/merged_neutrophils.rds
#
# ident.1 is male throughout, so a positive avg_log2FC means higher in males.
# ---------------------------------------------------------------------------

# Run from the project root ("Rscript scripts/04_sex_differences.R") or from
# inside scripts/; setup.R locates the project root from either.
source(if (file.exists("R/setup.R")) "R/setup.R" else "../R/setup.R")
cfg <- init_project()
load_modules()
require_packages("clusterProfiler", "org.Mm.eg.db")

age_levels <- cfg$analysis$age_levels
bm <- read_object(cfg, "bm_merged.rds")

# Differentiated neutrophils only. Cells with no potency call sit outside the
# stem/neutrophil subset that step 2 scored, and are excluded by the NA rule
# inside select_cells().
require_metadata(bm, c("singleR_main_label", "CytoTRACE2_Potency", "age", "sex"),
                 context = "step 4")

neutrophils <- select_cells(bm, list(
  "singleR_main_label is Neutrophils"  = bm$singleR_main_label == "Neutrophils",
  "CytoTRACE2_Potency is Differentiated" = bm$CytoTRACE2_Potency == "Differentiated",
  "age is one of analysis.age_levels"  = bm$age %in% age_levels
), context = "differentiated neutrophils")

print(table(neutrophils$sex, neutrophils$age))

neutrophils <- join_layers(neutrophils)
Seurat::DefaultAssay(neutrophils) <- "RNA"
neutrophils <- Seurat::NormalizeData(neutrophils, verbose = FALSE)
Seurat::Idents(neutrophils) <- neutrophils$age_sex
save_object(neutrophils, cfg, "merged_neutrophils.rds")

# --- 1. Sex contrast within each age, on each modality --------------------
degs <- list()
for (assay in c("RNA", "ADT")) {
  Seurat::DefaultAssay(neutrophils) <- assay
  degs[[assay]] <- list()

  group_sizes <- table(neutrophils$age_sex)
  for (age in age_levels) {
    idents <- paste0(age, c("_M", "_F"))
    present <- idents[idents %in% names(group_sizes)]
    if (length(present) < 2 || any(group_sizes[present] < 3)) {
      log_step("skipping ", age, " (", assay, "): ",
               paste(idents, ifelse(idents %in% present, group_sizes[idents], 0),
                     sep = "=", collapse = ", "))
      next
    }

    log_step("DE: ", age, " male vs female (", assay, ")")
    res <- Seurat::FindMarkers(neutrophils,
                               ident.1 = paste0(age, "_M"),
                               ident.2 = paste0(age, "_F"),
                               min.pct = 0.05,
                               logfc.threshold = 0,
                               min.cells.feature = 3,
                               min.cells.group = 3,
                               verbose = FALSE)
    res$gene <- rownames(res)     # keep names as a column, not just rownames
    degs[[assay]][[age]] <- res
    write_table(res, cfg, sprintf("degs_%s_%s.csv", age, assay))
  }
}
Seurat::DefaultAssay(neutrophils) <- "RNA"

# --- 2. GO over-representation, per age and direction ---------------------
# One gene per symbol (the largest absolute fold change), and the same
# universe for both directions: every gene that was tested.
ora_for_age <- function(res) {
  dedup <- res |>
    dplyr::group_by(.data$gene) |>
    dplyr::slice_max(abs(.data$avg_log2FC), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  universe <- dedup$gene
  directions <- list(
    male   = dedup$gene[dedup$p_val_adj < 0.05 & dedup$avg_log2FC > 0],
    female = dedup$gene[dedup$p_val_adj < 0.05 & dedup$avg_log2FC < 0]
  )

  lapply(directions, function(genes) {
    if (length(genes) < 5) {
      message("  only ", length(genes), " significant genes -- skipping ORA")
      return(NULL)
    }
    ego <- clusterProfiler::enrichGO(gene = genes, universe = universe,
                                     OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                                     keyType = "SYMBOL", ont = "BP",
                                     pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    if (is.null(ego)) return(NULL)
    ego@result[ego@result$pvalue < 0.05, ]
  })
}

for (age in names(degs[["RNA"]])) {
  log_step("GO ORA: ", age)
  terms <- ora_for_age(degs[["RNA"]][[age]])
  for (direction in names(terms)) {
    if (is.null(terms[[direction]])) next
    write_table(terms[[direction]], cfg,
                sprintf("go_ora_%s_up_in_%s.csv", age, direction))
  }
}

log_step("step 4 complete")
