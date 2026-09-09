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

#' The potency labels CytoTRACE2 assigns, least to most potent. Used to check
#' that the column still holds labels rather than factor codes.
CYTOTRACE_POTENCY_LEVELS <- c("Differentiated", "Unipotent", "Oligopotent",
                              "Multipotent", "Pluripotent", "Totipotent")

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

# ---------------------------------------------------------------------------
# Stage assignment from surface protein.
#
# Assigning maturation stage by clustering RNA and then testing RNA
# differences within those stages puts the confound inside the stratification:
# depth drives the label and the result. The ADT layer is the way out. It is
# CLR-normalised per cell, its isotype controls are non-significant, and it is
# the modality this experiment already paid for.
#
# It also gives a stage assignment that does not move when Seurat's cluster
# numbering changes between versions, which the cluster-id mapping did.
# ---------------------------------------------------------------------------

#' Map canonical marker names onto whatever the ADT panel actually calls them.
#'
#' Vendors name the same antibody differently (CXCR2 / CD182, Ly-6G / Ly6G), so
#' the config carries patterns rather than exact names and this resolves them
#' against the assay. Unresolved markers are reported, not silently dropped:
#' a panel scored on half its markers is worse than one that refuses to score.
resolve_adt_markers <- function(obj, cfg, assay = "ADT") {
  available <- rownames(obj[[assay]])
  aliases <- cfg$stage_assignment$adt_marker_aliases

  resolved <- vapply(names(aliases), function(canonical) {
    patterns <- unlist(aliases[[canonical]])
    hits <- unlist(lapply(patterns, function(p)
      grep(p, available, value = TRUE, ignore.case = TRUE)))
    if (length(hits)) hits[1] else NA_character_
  }, character(1))

  missing <- names(resolved)[is.na(resolved)]
  if (length(missing))
    warning("ADT markers not found in the panel: ", paste(missing, collapse = ", "),
            "\nPanel contains: ", paste(available, collapse = ", "),
            "\nAdjust stage_assignment.adt_marker_aliases in the config.")

  log_step("resolved ADT markers:")
  for (nm in names(resolved))
    log_step(sprintf("  %-10s -> %s", nm, resolved[[nm]] %||% "NOT FOUND"))
  resolved
}

#' Score each stage panel per cell from the ADT data.
#'
#' A panel's score is the mean of its scaled "high" markers minus the mean of
#' its scaled "low" markers. Markers are z-scored across cells first so that
#' antibodies with different dynamic ranges contribute comparably; the ADT data
#' layer is already isotype-centred per cell, so this scaling is across cells,
#' not within them.
score_adt_panels <- function(obj, cfg, assay = "ADT") {
  resolved <- resolve_adt_markers(obj, cfg, assay)
  panels <- cfg$stage_assignment$adt_panels

  adt <- as.matrix(Seurat::GetAssayData(obj, assay = assay, layer = "data"))
  scaled <- t(scale(t(adt)))
  scaled[!is.finite(scaled)] <- 0

  panel_score <- function(panel) {
    take <- function(side) {
      markers <- stats::na.omit(unname(resolved[unlist(panel[[side]])]))
      markers <- intersect(markers, rownames(scaled))
      if (!length(markers)) return(rep(0, ncol(scaled)))
      colMeans(scaled[markers, , drop = FALSE])
    }
    take("high") - take("low")
  }

  scores <- vapply(panels, panel_score, numeric(ncol(obj)))
  rownames(scores) <- colnames(obj)
  scores
}

#' Assign each cell the stage whose panel scores highest.
#'
#' `margin` is the gap between the best and second-best score. A cell whose top
#' two stages are indistinguishable is labelled NA rather than assigned by a
#' coin flip, and the fraction of such cells is reported: it says how well the
#' panel actually separates the stages.
assign_stage_adt <- function(obj, cfg, min_margin = 0.1, to = "stage") {
  scores <- score_adt_panels(obj, cfg)
  stage_levels <- cfg$analysis$stage_levels
  scores <- scores[, intersect(stage_levels, colnames(scores)), drop = FALSE]

  ordered <- t(apply(scores, 1, function(x) sort(x, decreasing = TRUE)))
  best <- colnames(scores)[apply(scores, 1, which.max)]
  margin <- ordered[, 1] - ordered[, 2]

  best[margin < min_margin] <- NA_character_
  obj[[to]] <- factor(best, levels = stage_levels)
  obj[[paste0(to, "_margin")]] <- margin

  log_step(sprintf("ADT stage assignment: %d of %d cells assigned (%.1f%% ambiguous below margin %.2f)",
                   sum(!is.na(best)), length(best), 100 * mean(is.na(best)), min_margin))
  print(table(obj[[to]][, 1], useNA = "ifany"))
  obj
}

#' Assign stage by rank-based scoring of RNA signatures.
#'
#' UCell ranks genes within each cell before scoring, so the score depends on
#' the ordering of genes rather than their absolute counts. That makes it far
#' less depth-sensitive than a mean-expression module score -- though it is
#' still RNA, so it does not break the circularity the way the ADT route does.
assign_stage_markers <- function(obj, cfg, assay = "RNA", to = "stage_rna") {
  require_packages("UCell")
  signatures <- lapply(cfg$stage_assignment$rna_signatures, unlist)

  obj <- UCell::AddModuleScore_UCell(obj, features = signatures, assay = assay,
                                     name = "_UCell")
  score_cols <- paste0(names(signatures), "_UCell")
  present <- intersect(score_cols, colnames(obj@meta.data))
  scores <- as.matrix(obj@meta.data[, present, drop = FALSE])

  best <- sub("_UCell$", "", present)[apply(scores, 1, which.max)]
  obj[[to]] <- factor(best, levels = cfg$analysis$stage_levels)
  log_step("UCell stage assignment:")
  print(table(obj[[to]][, 1], useNA = "ifany"))
  obj
}

#' Stage labels from the old cluster-id mapping, for comparison only.
assign_stage_clusters <- function(obj, cfg, cluster_col = "seurat_clusters",
                                  to = "stage_clusters") {
  map <- unlist(cfg$analysis$neutrophil_cluster_labels)
  if (is.null(map)) {
    log_step("no cluster->stage map in the config; skipping the cluster assignment")
    return(obj)
  }
  clusters <- as.character(obj[[cluster_col]][, 1])
  obj[[to]] <- factor(unname(map[clusters]), levels = cfg$analysis$stage_levels)
  obj
}

#' Confusion matrix between two stage assignments.
#'
#' Run as a formal comparison rather than a spot check: if protein-defined and
#' RNA-cluster-defined stages disagree substantially, that disagreement is
#' itself a result, and it decides which stratification the analyses should use.
stage_confusion <- function(obj, a = "stage", b = "stage_clusters") {
  missing <- setdiff(c(a, b), colnames(obj@meta.data))
  if (length(missing)) {
    log_step("cannot compare stage assignments; missing: ", paste(missing, collapse = ", "))
    return(NULL)
  }

  tab <- table(as.character(obj[[a]][, 1]), as.character(obj[[b]][, 1]),
               dnn = c(a, b), useNA = "ifany")
  both <- !is.na(obj[[a]][, 1]) & !is.na(obj[[b]][, 1])
  agreement <- mean(as.character(obj[[a]][, 1])[both] ==
                      as.character(obj[[b]][, 1])[both])

  log_step(sprintf("stage assignments agree on %.1f%% of the %d cells labelled by both",
                   100 * agreement, sum(both)))
  print(tab)
  structure(list(table = tab, agreement = agreement, n_compared = sum(both)),
            class = "stage_confusion")
}

#' Dispatch on stage_assignment.method, and always compute the comparison.
add_stage_labels <- function(obj, cfg, cluster_col = "seurat_clusters", to = "stage") {
  method <- cfg$stage_assignment$method %||% "clusters"
  log_step("stage assignment method: ", method)

  obj <- switch(method,
    adt      = assign_stage_adt(obj, cfg, to = to),
    markers  = assign_stage_markers(obj, cfg, to = to),
    clusters = assign_stage_clusters(obj, cfg, cluster_col, to = to),
    stop("unknown stage_assignment.method: ", method))

  comparison <- cfg$stage_assignment$compare_against
  if (!is.null(comparison) && !identical(comparison, method)) {
    obj <- switch(comparison,
      adt      = assign_stage_adt(obj, cfg, to = "stage_adt"),
      markers  = assign_stage_markers(obj, cfg, to = "stage_rna"),
      clusters = assign_stage_clusters(obj, cfg, cluster_col, to = "stage_clusters"),
      obj)
  }
  obj
}
