# ---------------------------------------------------------------------------
# Protein-gated maturation stage, including the cells RNA QC removed.
#
# The question this exists for. The male library loses immature and mature
# neutrophils with age (22% -> 7% immature, 14% -> 8% mature, 3m -> 18m) --
# statistically the strongest age trend in the data. But those are the cells
# with the least RNA (male mature: 593 UMIs at 3m, 328 at 18m), RNA QC removes
# the least-RNA cells, and the stage labels themselves come from RNA. So the
# trend could be biology, or it could be aged mature neutrophils failing QC.
# The RNA data cannot tell those apart.
#
# Surface protein can. Antibody capture does not depend on how much mRNA a
# cell yields, and step 1 keeps the protein record of every hashtagged cell,
# QC failures included. Gating those cells on Ly6G / CD11b / CXCR2 gives a
# stage composition that does not depend on RNA at all, and the retention
# table says directly how many protein-defined mature cells RNA QC took at
# each age. If protein-gated mature cells are there and fail QC, the male
# trend is QC loss; if they are absent before QC, it is in the marrow.
#
# Limits, stated where the numbers are. Cells Cell Ranger never called are
# invisible here too. The panel carries no CD101 or CD117, so protein stage is
# coarser than RNA stage: Ly6G-negative precursors, Ly6G+ CXCR2- (immature-
# like) and Ly6G+ CXCR2+ (mature-like). Thresholds are set per library, since
# each library was stained separately; comparisons across ages within a
# library share a threshold, which is what makes them clean.
# ---------------------------------------------------------------------------

#' Canonical marker -> patterns to look for in the ADT rownames.
PROTEIN_GATE_MARKERS <- list(
  Ly6G   = c("Ly.?6G", "Ly6G"),
  CD11b  = c("CD11b", "ITGAM"),
  CXCR2  = c("CXCR2", "CD182"),
  CD117  = c("CD117", "c.?Kit"),
  SiglecF = c("Siglec.?F", "CD170"),
  FceRIa = c("FceRI", "FcεRI", "FcERI"),
  CD115  = c("CD115", "CSF1R")
)

#' Granulocytic protein stages, in maturation order.
PROTEIN_STAGES <- c("Ly6G- precursor", "Ly6G+ CXCR2- (immature-like)",
                    "Ly6G+ CXCR2+ (mature-like)")

#' Every protein class a cell can receive, granulocytic stages first.
PROTEIN_CLASSES <- c(PROTEIN_STAGES, "eosinophil-like (Siglec-F+)",
                     "basophil-like (FceRIa+)", "monocyte-like (CD115+)",
                     "CD11b- (non-myeloid)")

#' Resolve each canonical marker against the panel. Absent markers are NA.
resolve_gate_markers <- function(features, markers = PROTEIN_GATE_MARKERS) {
  vapply(names(markers), function(m) {
    hits <- unlist(lapply(markers[[m]], function(p)
      grep(p, features, value = TRUE, ignore.case = TRUE)))
    hits <- setdiff(hits, grep("isotype", hits, value = TRUE, ignore.case = TRUE))
    if (length(hits)) hits[1] else NA_character_
  }, character(1))
}

#' The background level: a high quantile of the isotype controls.
#'
#' The ADT data are log1p counts minus each cell's median isotype signal, so
#' isotypes sit near zero. A marker above the 99th percentile of every isotype
#' value in the library is above anything non-specific binding produces.
isotype_background <- function(adt, isotype_pattern, q = 0.99) {
  iso <- grep(isotype_pattern, rownames(adt), value = TRUE)
  if (!length(iso)) return(NA_real_)
  stats::quantile(as.numeric(adt[iso, , drop = FALSE]), q, names = FALSE, na.rm = TRUE)
}

#' Density valley between the two dominant modes of x, if there is one.
#'
#' The standard way to split a bimodal flow channel. Returns NA when the
#' distribution has one mode, so the caller falls back to background rather
#' than inventing a split.
valley_threshold <- function(x, min_separation = 0.5, adjust = 1) {
  x <- x[is.finite(x)]
  if (length(x) < 50) return(NA_real_)
  d <- stats::density(x, adjust = adjust, n = 512)
  y <- d$y
  peaks <- which(diff(sign(diff(y))) == -2) + 1
  if (length(peaks) < 2) return(NA_real_)
  # The two tallest peaks far enough apart to be different populations.
  peaks <- peaks[order(-y[peaks])]
  first <- peaks[1]
  second <- peaks[-1][abs(d$x[peaks[-1]] - d$x[first]) >= min_separation][1]
  if (is.na(second)) return(NA_real_)
  lo <- min(first, second); hi <- max(first, second)
  # A shallow dip is not a valley: require the minimum to sit well below the
  # smaller peak.
  v <- lo - 1 + which.min(y[lo:hi])
  if (y[v] > 0.8 * min(y[first], y[second])) return(NA_real_)
  d$x[v]
}

#' A positive threshold for one marker: the density valley when the channel
#' is bimodal and the valley clears background, otherwise background.
marker_threshold <- function(x, background) {
  v <- valley_threshold(x)
  if (is.finite(v) && (!is.finite(background) || v > background))
    return(list(threshold = v, method = "density valley"))
  list(threshold = background, method = "isotype background (99th pct)")
}

#' Gate one library's cells into protein classes.
#'
#' @param record the list step 1 saves: `adt` (antibodies x cells, isotype-
#'   centred) and `meta` (one row per cell, with `qc_status`).
#' @return list(cells = per-cell table, thresholds = per-marker table)
gate_library <- function(record, cfg) {
  adt <- record$adt
  meta <- record$meta
  stopifnot(identical(colnames(adt), meta$barcode))

  markers <- resolve_gate_markers(rownames(adt))
  if (is.na(markers[["Ly6G"]]) || is.na(markers[["CD11b"]]))
    stop("the protein gate needs Ly6G and CD11b; the panel has: ",
         paste(rownames(adt), collapse = ", "))

  # Doublets are not cells, so they set no threshold and get no class.
  usable <- meta$qc_status != "doublet"
  bg <- isotype_background(adt[, usable, drop = FALSE],
                           cfg$adt$isotype_pattern %||% "^Isotype")
  value <- function(m) if (is.na(markers[[m]])) NULL else adt[markers[[m]], ]

  thresholds <- list()
  set_threshold <- function(m, x) {
    t <- marker_threshold(x, bg)
    feature <- markers[[if (m %in% names(markers)) m else sub("-hi$", "", m)]]
    thresholds[[m]] <<- data.frame(marker = m, feature = feature,
                                   threshold = t$threshold, method = t$method,
                                   background = bg, stringsAsFactors = FALSE)
    t$threshold
  }

  cd11b <- value("CD11b")
  ly6g <- value("Ly6G")
  t_cd11b <- set_threshold("CD11b", cd11b[usable])
  myeloid <- cd11b > t_cd11b
  # Ly6G is split within myeloid cells, where it is bimodal (neutrophils vs
  # everything else); across all cells the lymphoid mass hides the valley.
  t_ly6g <- set_threshold("Ly6G", ly6g[usable & myeloid])
  ly6g_pos <- ly6g > t_ly6g

  cxcr2 <- value("CXCR2")
  cxcr2_pos <- if (is.null(cxcr2)) NULL else
    cxcr2 > set_threshold("CXCR2", cxcr2[usable & myeloid & ly6g_pos])

  flag <- function(m, among) {
    x <- value(m)
    if (is.null(x)) return(rep(FALSE, ncol(adt)))
    x > set_threshold(m, x[usable & among])
  }
  siglecf <- flag("SiglecF", myeloid & !ly6g_pos)
  fceri <- flag("FceRIa", myeloid & !ly6g_pos)
  cd115 <- flag("CD115", myeloid & !ly6g_pos)

  ly6g_split <- if (!is.null(cxcr2_pos)) cxcr2_pos else {
    # No CXCR2 on the panel: split Ly6G+ cells into int / hi instead.
    hi <- ly6g > set_threshold("Ly6G-hi", ly6g[usable & myeloid & ly6g_pos])
    hi
  }

  cls <- ifelse(!myeloid, "CD11b- (non-myeloid)",
         ifelse(ly6g_pos & ly6g_split, PROTEIN_STAGES[3],
         ifelse(ly6g_pos, PROTEIN_STAGES[2],
         ifelse(siglecf, "eosinophil-like (Siglec-F+)",
         ifelse(fceri, "basophil-like (FceRIa+)",
         ifelse(cd115, "monocyte-like (CD115+)", PROTEIN_STAGES[1]))))))
  cls[!usable] <- NA_character_
  if (is.null(cxcr2_pos)) {
    cls[cls == PROTEIN_STAGES[2]] <- "Ly6G-int (immature-like)"
    cls[cls == PROTEIN_STAGES[3]] <- "Ly6G-hi (mature-like)"
  }

  cells <- data.frame(meta, protein_class = cls,
                      passed_rna_qc = meta$qc_status == "pass",
                      Ly6G = ly6g, CD11b = cd11b,
                      CXCR2 = if (is.null(cxcr2)) NA_real_ else cxcr2,
                      stringsAsFactors = FALSE)
  list(cells = cells, thresholds = do.call(rbind, thresholds),
       markers = markers, has_cxcr2 = !is.null(cxcr2_pos))
}

#' The granulocytic stage levels a gated library actually uses.
protein_stage_levels <- function(has_cxcr2) {
  if (has_cxcr2) PROTEIN_STAGES
  else c(PROTEIN_STAGES[1], "Ly6G-int (immature-like)", "Ly6G-hi (mature-like)")
}

#' How many cells of each protein class RNA QC kept, per sex and age.
#'
#' The decisive table. `retained_fraction` for mature-like cells falling with
#' age, in one library, is QC removing aged mature neutrophils -- and the RNA
#' composition trend inherits it.
protein_retention <- function(cells, age_levels, class_levels) {
  cells <- cells[!is.na(cells$protein_class), ]
  key <- interaction(cells$sex, factor(cells$age, levels = age_levels),
                     factor(cells$protein_class, levels = class_levels),
                     drop = TRUE, sep = "|")
  rows <- lapply(split(seq_len(nrow(cells)), key), function(i) {
    data.frame(sex = cells$sex[i[1]], age = cells$age[i[1]],
               protein_class = cells$protein_class[i[1]],
               n_before_qc = length(i), n_passed = sum(cells$passed_rna_qc[i]),
               median_umi = stats::median(cells$nCount_RNA[i]),
               median_umi_failed = if (any(!cells$passed_rna_qc[i]))
                 stats::median(cells$nCount_RNA[i][!cells$passed_rna_qc[i]]) else NA_real_,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$retained_fraction <- out$n_passed / out$n_before_qc
  ci <- wilson_interval(out$n_passed, out$n_before_qc)
  out$ci_low <- ci$lower
  out$ci_high <- ci$upper
  out$age <- factor(out$age, levels = age_levels)
  out$protein_class <- factor(out$protein_class, levels = class_levels)
  out <- out[order(out$sex, out$protein_class, out$age), ]
  out$age <- as.character(out$age)
  out$protein_class <- as.character(out$protein_class)
  rownames(out) <- NULL
  out
}

#' Protein-stage composition before and after RNA QC, stacked.
protein_composition <- function(cells, cfg, stage_levels) {
  gran <- cells[!is.na(cells$protein_class) & cells$protein_class %in% stage_levels, ]
  sets <- list("all singlets (before RNA QC)" = rep(TRUE, nrow(gran)),
               "passed RNA QC" = gran$passed_rna_qc)
  do.call(rbind, lapply(names(sets), function(nm) {
    keep <- sets[[nm]]
    comp <- composition_from_labels(gran$protein_class[keep], gran$age[keep],
                                    gran$sex[keep], levels = stage_levels,
                                    age_levels = cfg$analysis$age_levels,
                                    sex_levels = cfg$analysis$sex_levels)
    comp$cell_set <- nm
    comp
  }))
}

#' Protein class against the RNA stage label, for cells that have both.
protein_vs_rna <- function(cells, rna_stage) {
  both <- !is.na(cells$protein_class) & !is.na(rna_stage)
  if (!any(both)) return(NULL)
  tab <- as.data.frame(table(sex = cells$sex[both],
                             protein_class = cells$protein_class[both],
                             rna_stage = rna_stage[both]),
                       responseName = "n", stringsAsFactors = FALSE)
  tab <- tab[tab$n > 0, ]
  totals <- stats::aggregate(n ~ sex + rna_stage, tab, sum)
  names(totals)[3] <- "n_rna_stage"
  tab <- merge(tab, totals, by = c("sex", "rna_stage"))
  tab$fraction_of_rna_stage <- tab$n / tab$n_rna_stage
  tab[order(tab$sex, tab$rna_stage, -tab$n), ]
}
