# ---------------------------------------------------------------------------
# Cell identity: SingleR reference labelling, readable label names,
# CytoTRACE2 potency, and neutrophil maturation module scores.
# ---------------------------------------------------------------------------

#' Label cells against the ImmGen reference at both main and fine resolution.
#'
#' Cells SingleR cannot label confidently are pruned to NA; they are relabelled
#' "NA" so they show up explicitly in plots and tables rather than vanishing.
annotate_singler <- function(obj, cfg) {
  ref <- celldex::fetchReference(cfg$annotation$singler_ref,
                                 cfg$annotation$singler_ref_version)
  sce <- Seurat::as.SingleCellExperiment(obj, assay = "RNA")

  main <- SingleR::SingleR(test = sce, ref = ref, labels = ref$label.main)
  fine <- SingleR::SingleR(test = sce, ref = ref, labels = ref$label.fine)

  obj$singleR_main_label <- main$pruned.labels
  obj$singleR_fine_label <- fine$pruned.labels
  # `scores` is a cell x label matrix; only the winning score is a per-cell
  # value, so that is what goes into metadata. A cell with no finite score at
  # all gets NA rather than -Inf.
  obj$main_label_score <- apply(main$scores, 1, function(x) {
    x <- x[is.finite(x)]
    if (length(x)) max(x) else NA_real_
  })

  obj$singleR_main_label[is.na(obj$singleR_main_label)] <- "NA"
  obj$singleR_fine_label[is.na(obj$singleR_fine_label)] <- "NA"
  obj
}

#' ImmGen fine labels -> names a reader can interpret.
IMMGEN_LABEL_MAP <- c(
  "Stem cells (SC.MEP)"         = "Megakaryocyte-Erythroid Progenitor (MEP)",
  "Stem cells (MLP)"            = "Multipotent Lymphoid Progenitor",
  "Stem cells (SC.MDP)"         = "Macrophage-DC Progenitor (MDP)",
  "Stem cells (SC.CMP.DR)"      = "Common Myeloid Progenitor (CMP)",
  "Stem cells (GMP)"            = "Granulocyte-Monocyte Progenitor (GMP)",
  "Stem cells (SC.CDP)"         = "Common Dendritic Cell Progenitor (CDP)",
  "Stem cells (SC.ST34F)"       = "Short-Term HSC (CD34+ Flt3-)",
  "Stem cells (SC.CD150-CD48-)" = "HSC subset (CD150- CD48-, SLAM MPP)",
  "Stem cells (SC.LT34F)"       = "Long-Term HSC (CD34- Flt3-)",
  "Stem cells (SC.MPP34F)"      = "Multipotent Progenitor (CD34+ Flt3-)",
  "Stem cells (SC.STSL)"        = "Short-Term HSC (Sca1+ Lin-, SLAM)",
  "Stem cells (LTHSC)"          = "Long-Term HSC (LT-HSC)",
  "Stem cells (proB.CLP)"       = "Common Lymphoid Progenitor/pro-B cell",
  "Neutrophils (GN)"            = "Neutrophils",
  "Neutrophils (GN.ARTH)"       = "Neutrophils",
  "Neutrophils (GN.URAC)"       = "Neutrophils",
  "Neutrophils (GN.Thio)"       = "Neutrophils"
)

#' Add the readable fine label. Labels with no mapping keep their ImmGen name
#' rather than becoming NA.
add_readable_labels <- function(obj, map = IMMGEN_LABEL_MAP,
                                from = "singleR_fine_label",
                                to = "fine_label_readable") {
  raw <- as.character(obj[[from]][, 1])
  mapped <- unname(map[raw])
  obj[[to]] <- ifelse(is.na(mapped), raw, mapped)
  obj
}

#' Columns CytoTRACE2 is expected to contribute. Steps 4 and 7 select
#' differentiated cells on CytoTRACE2_Potency, so its absence is fatal there.
CYTOTRACE_COLUMNS <- c("CytoTRACE2_Score", "CytoTRACE2_Potency",
                       "CytoTRACE2_Relative", "preKNN_CytoTRACE2_Score",
                       "preKNN_CytoTRACE2_Potency")

#' Developmental potency per cell (CytoTRACE2).
#'
#' The result is verified rather than trusted: if CytoTRACE2 returns columns
#' under different names, or covers only some cells, AddMetaData succeeds
#' quietly and the problem only surfaces two steps later as an empty subset.
run_cytotrace2 <- function(obj, cfg, assay = "RNA") {
  expr <- as.matrix(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
  result <- CytoTRACE2::cytotrace2(expr, is_seurat = FALSE,
                                   species = cfg$annotation$cytotrace_species)

  if (!"CytoTRACE2_Potency" %in% colnames(result))
    stop("CytoTRACE2 returned no CytoTRACE2_Potency column. It returned: ",
         paste(colnames(result), collapse = ", "),
         "
Steps 4 and 7 select cells on that column, so the pipeline cannot ",
         "continue without it. Check the CytoTRACE2 version.")

  overlap <- length(intersect(rownames(result), colnames(obj)))
  if (overlap < ncol(obj))
    warning("CytoTRACE2 returned values for ", overlap, " of ", ncol(obj),
            " cells; the rest will carry NA potency")

  obj <- Seurat::AddMetaData(obj, metadata = result)

  scored <- sum(!is.na(obj$CytoTRACE2_Potency))
  log_step(sprintf("CytoTRACE2: %d of %d cells scored", scored, ncol(obj)))
  if (scored == 0)
    stop("CytoTRACE2 produced no usable potency calls after AddMetaData. ",
         "This is usually a cell-name mismatch between the result and the object.")
  print(table(obj$CytoTRACE2_Potency, useNA = "ifany"))
  obj
}

#' Marker modules for the granulocyte maturation series.
NEUTROPHIL_MODULES <- list(
  proNeu  = c("Elane", "Mpo", "Prtn3", "Ctsg", "Ms4a3", "Cebpe", "Gfi1",
              "Fcnb", "Rab44", "Nkg7", "Plac8", "Cd34", "Kit", "Srgn"),
  proNeu2 = c("Il5ra", "Cebpe", "Fcnb", "Ltf", "Camp", "Mki67"),
  preNeu  = c("Ltf", "Camp", "Ngp", "Lcn2", "Chil3", "Cebpe", "Fcnb",
              "Anxa1", "Hp", "Ifitm6", "Mki67"),
  immNeu  = c("Ltf", "Camp", "Lcn2", "Mmp8", "Cd177", "Ifitm1",
              "Chil3", "Cd101"),
  matNeu  = c("Cxcr2", "Sell", "Il1b", "Csf3r", "Mmp9", "S100a8", "S100a9",
              "Retnlg", "Ifitm1", "Msrb1", "Slpi", "Fpr1"),
  cycling = c("Mki67", "Top2a", "Ccnb1", "Ccna2", "Birc5", "Ube2c", "Cdk1")
)

#' Score the maturation modules and give the columns their module names.
#'
#' AddModuleScore appends a numeric suffix (MOD_1..MOD_n) in list order, which
#' is renamed here so downstream code refers to `preNeu_score` and friends.
add_module_scores <- function(obj, modules = NEUTROPHIL_MODULES, ctrl = 50, seed = 42) {
  obj <- Seurat::AddModuleScore(obj, features = modules, name = "MOD_",
                                ctrl = ctrl, seed = seed)
  generated <- paste0("MOD_", seq_along(modules))
  idx <- match(generated, colnames(obj@meta.data))
  stopifnot(!anyNA(idx))
  colnames(obj@meta.data)[idx] <- paste0(names(modules), "_score")
  obj
}

#' Map merged-object cluster ids to developmental stage names from the config.
#'
#' Cluster numbering is not stable across Seurat versions: if the merge is
#' re-run, check the cluster/module plots before trusting this mapping.
add_stage_labels <- function(obj, cfg, cluster_col = "seurat_clusters",
                             to = "fine_neu_labels") {
  map <- unlist(cfg$analysis$neutrophil_cluster_labels)
  clusters <- as.character(obj[[cluster_col]][, 1])
  unmapped <- setdiff(unique(clusters), names(map))
  if (length(unmapped))
    warning("clusters with no stage label in config: ", paste(unmapped, collapse = ", "))
  obj[[to]] <- factor(unname(map[clusters]), levels = cfg$analysis$stage_levels)
  obj
}
