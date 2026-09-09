# ---------------------------------------------------------------------------
# Depth matching.
#
# The two libraries differ several-fold in capture per cell, and sex is
# confounded with library, so any differential result computed on unmatched
# depth is uninterpretable: the deeper library detects more of everything.
# Matching is therefore the primary analysis path, not a sensitivity check.
#
# Downsampling is binomial thinning of the counts, which is the right operation
# because it reproduces what a shallower run would have measured. Normalisation
# does not substitute for it: LogNormalize rescales totals but cannot restore
# a transcript that was never captured, which is why detection rates -- and so
# anything computed from them -- stay biased after normalisation.
# ---------------------------------------------------------------------------

#' Median UMIs per cell, per group.
group_depth <- function(counts, groups) {
  totals <- Matrix::colSums(counts)
  vapply(split(totals, groups), stats::median, numeric(1))
}

#' Ratio of the deepest group's median depth to the shallowest.
depth_ratio <- function(counts, groups) {
  depths <- group_depth(counts, groups)
  max(depths) / min(depths)
}

#' Downsample every cell toward a common depth.
#'
#' @param target one of "median" (the lowest group median) or a number of UMIs.
#' @return the thinned count matrix, with the target recorded as an attribute.
downsample_counts <- function(counts, groups, target = "median", seed = 42) {
  require_packages("DropletUtils")
  set.seed(seed)

  totals <- Matrix::colSums(counts)
  target_depth <- if (identical(target, "median")) min(group_depth(counts, groups)) else as.numeric(target)

  # Cells already at or below the target are left alone: thinning them further
  # would discard real signal to no purpose.
  prop <- pmin(1, target_depth / pmax(totals, 1))
  thinned <- DropletUtils::downsampleMatrix(counts, prop = prop, bycol = TRUE)

  log_step(sprintf("depth matching: target %d UMIs; %d of %d cells thinned",
                   round(target_depth), sum(prop < 1), length(prop)))
  attr(thinned, "target_depth") <- target_depth
  attr(thinned, "depth_ratio_before") <- depth_ratio(counts, groups)
  attr(thinned, "depth_ratio_after") <- depth_ratio(thinned, groups)
  thinned
}

#' Attach a depth-matched counts layer to a Seurat object.
#'
#' Written to a separate assay rather than overwriting RNA, so the unmatched
#' counts stay available for the comparison that shows what matching changed.
add_matched_assay <- function(obj, cfg, group_col = "sex", assay_name = "RNAmatched") {
  counts <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
  groups <- as.character(obj[[group_col]][, 1])

  before <- depth_ratio(counts, groups)
  log_step(sprintf("depth ratio between %s groups before matching: %.2fx",
                   group_col, before))

  thinned <- downsample_counts(counts, groups, cfg$depth$target, cfg$seed)
  obj[[assay_name]] <- Seurat::CreateAssayObject(counts = thinned)
  Seurat::DefaultAssay(obj) <- assay_name
  obj <- Seurat::NormalizeData(obj, assay = assay_name, verbose = FALSE)
  Seurat::DefaultAssay(obj) <- "RNA"

  log_step(sprintf("depth ratio after matching: %.2fx",
                   attr(thinned, "depth_ratio_after")))
  obj
}

#' The three controls that separate a real effect from a depth artefact.
#'
#' Each answers a different question, and a result needs all three to be
#' interpretable:
#'   cross_sex         does the effect survive matching the deeper library down?
#'   within_sex_null   does downsampling alone manufacture the effect, inside
#'                     one library where no sex difference can exist?
#'   cell_number       is it a power artefact rather than a depth one?
#'
#' @param analysis a function of (counts, groups) returning a named numeric
#'   vector of whatever the caller considers the result.
depth_controls <- function(counts, groups, analysis, cfg, seed = 42) {
  out <- list()

  out$observed <- analysis(counts, groups)

  if (isTRUE(cfg$depth$controls$cross_sex_downsample)) {
    log_step("  control: matching the deeper group down")
    out$depth_matched <- analysis(downsample_counts(counts, groups, "median", seed), groups)
  }

  if (isTRUE(cfg$depth$controls$within_sex_null_downsample)) {
    # Inside one library, split at random and thin one half. Any "difference"
    # this produces is manufactured by the thinning itself.
    log_step("  control: downsampling within one library (no real difference exists)")
    deeper <- names(which.max(group_depth(counts, groups)))
    cells <- which(groups == deeper)
    sub <- counts[, cells, drop = FALSE]

    set.seed(seed)
    pseudo <- rep("b", length(cells))
    pseudo[sample(length(cells), length(cells) %/% 2)] <- "a"

    # Thin one pseudo-half by the SAME ratio that separates the real groups, so
    # this control asks exactly the question the real comparison does, in a
    # setting where the true answer is known to be "no difference".
    ratio <- depth_ratio(counts, groups)
    totals <- Matrix::colSums(sub)
    prop <- ifelse(pseudo == "a", pmin(1, (stats::median(totals) / ratio) / pmax(totals, 1)), 1)

    set.seed(seed)
    thinned <- DropletUtils::downsampleMatrix(sub, prop = prop, bycol = TRUE)
    log_step(sprintf("    thinned pseudo-group by %.2fx; ratio now %.2fx",
                     ratio, depth_ratio(thinned, pseudo)))
    out$within_library_null <- analysis(thinned, pseudo)
  }

  if (isTRUE(cfg$depth$controls$match_cell_numbers)) {
    log_step("  control: matching cell numbers instead of depth")
    set.seed(seed)
    smallest <- min(table(groups))
    keep <- unlist(lapply(split(seq_along(groups), groups), function(i)
      sample(i, min(length(i), smallest))))
    out$cell_matched <- analysis(counts[, keep, drop = FALSE], groups[keep])
  }

  out
}

#' Did a result survive depth matching, and did it keep its sign?
#'
#' A result that changes sign or loses coherence under matching is a depth
#' artefact and should be dropped rather than reported with a caveat.
survives_matching <- function(observed, matched, tolerance = 0.5) {
  if (is.na(observed) || is.na(matched)) return(NA)
  if (sign(observed) != sign(matched)) return("sign-flip")
  abs(matched) >= tolerance * abs(observed)
}
