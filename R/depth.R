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

#' Downsample within strata, each toward its own target.
#'
#' One global target leaves every cell already below it untouched, so a
#' population that is shallow by nature -- mature neutrophils, at 330-590 UMIs
#' against a 2225 target in the male library -- is never matched at all, and
#' its age or sex comparisons run on counts several-fold apart while the global
#' depth gate reads 1.00x. Matching inside each stratum (a stage, or a stage
#' within a sex) sets the target from the groups actually being compared.
#'
#' Groups smaller than `min_cells` do not set the target: a median of four
#' cells is not a depth, and letting it set one would thin a whole stratum to
#' noise. They are thinned to it like everyone else.
#'
#' @param strata one label per cell; NA is treated as its own stratum.
downsample_within <- function(counts, groups, strata, min_cells = 1, seed = 42) {
  require_packages("DropletUtils")
  strata <- ifelse(is.na(strata), "<unassigned>", as.character(strata))
  totals <- Matrix::colSums(counts)
  prop <- rep(1, ncol(counts))
  targets <- list()

  for (st in unique(strata)) {
    idx <- which(strata == st)
    g <- groups[idx]
    sizes <- table(g)
    eligible <- names(sizes)[sizes >= min_cells]
    if (!length(eligible)) eligible <- names(sizes)
    med <- vapply(eligible, function(x) stats::median(totals[idx][g == x]), numeric(1))
    target <- min(med)
    prop[idx] <- pmin(1, target / pmax(totals[idx], 1))
    targets[[st]] <- data.frame(stratum = st, target = target,
                                set_by = names(which.min(med)),
                                n_cells = length(idx),
                                n_thinned = sum(prop[idx] < 1),
                                ratio_before = max(med) / min(med),
                                stringsAsFactors = FALSE)
  }

  set.seed(seed)
  thinned <- DropletUtils::downsampleMatrix(counts, prop = prop, bycol = TRUE)
  report <- do.call(rbind, targets)
  report$ratio_after <- vapply(report$stratum, function(st) {
    idx <- which(strata == st)
    g <- groups[idx]
    sizes <- table(g)
    keep <- names(sizes)[sizes >= min_cells]
    if (length(keep) < 2) return(NA_real_)
    sub_idx <- idx[g %in% keep]
    depth_ratio(thinned[, sub_idx, drop = FALSE], groups[sub_idx])
  }, numeric(1))
  attr(thinned, "targets") <- report
  thinned
}

#' Attach a depth-matched counts layer to a Seurat object.
#'
#' Written to a separate assay rather than overwriting RNA, so the unmatched
#' counts stay available for the comparison that shows what matching changed.
#'
#' @param within_cols metadata columns defining the strata matched separately
#'   (see downsample_within()). NULL matches all cells toward one target, the
#'   old behaviour.
add_matched_assay <- function(obj, cfg, group_col = "sex", assay_name = "RNAmatched",
                              within_cols = NULL) {
  counts <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
  groups <- as.character(obj[[group_col]][, 1])

  before <- depth_ratio(counts, groups)
  log_step(sprintf("%s: depth ratio between %s groups before matching: %.2fx",
                   assay_name, group_col, before))

  if (length(within_cols)) {
    strata <- do.call(paste, c(lapply(within_cols, function(col)
      as.character(obj[[col]][, 1])), sep = " | "))
    thinned <- downsample_within(counts, groups, strata,
                                 min_cells = cfg$gates$min_cells_per_stratum,
                                 seed = cfg$seed)
    targets <- attr(thinned, "targets")
    log_step("  matched within ", paste(within_cols, collapse = " x "),
             ", each stratum toward its own shallowest ", group_col, " group:")
    print(targets[order(targets$stratum), ], row.names = FALSE)
    attr(obj, paste0(assay_name, "_targets")) <- targets
  } else {
    thinned <- downsample_counts(counts, groups, cfg$depth$target, cfg$seed)
    log_step(sprintf("  depth ratio after matching: %.2fx",
                     attr(thinned, "depth_ratio_after")))
  }

  obj[[assay_name]] <- Seurat::CreateAssayObject(counts = thinned)
  Seurat::DefaultAssay(obj) <- assay_name
  obj <- Seurat::NormalizeData(obj, assay = assay_name, verbose = FALSE)
  Seurat::DefaultAssay(obj) <- "RNA"
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

#' Depth per age WITHIN each sex, optionally within stage as well.
#'
#' The within-sex age comparisons in step 9 are shifts relative to the first
#' age, so a constant depth difference between the two libraries cancels. What
#' does not cancel is depth varying BETWEEN HASHTAGS INSIDE one library: each
#' age is a separate hashtag in the same run, and both pseudotime position and
#' CytoTRACE2 potency track transcriptional complexity.
#'
#' POOLING STAGES MAKES THIS UNREADABLE. Complexity is a property of the cell
#' type -- GMPs carry far more detectable genes than mature neutrophils -- so
#' when the stage mix shifts with age, pooled complexity shifts with it and
#' nothing technical need be wrong. Males go from 29% GMP / 14% mature at 3m to
#' 43% / 8% at 18m, which on its own raises pooled genes detected. Only the
#' per-stage rows separate a real depth trend from a composition shift, so both
#' are computed and the per-stage rows are the ones to read.
#'
#' `umi_ratio` and `gene_ratio` are each age's median depth and median genes
#' detected over the reference age's. BOTH have to be near 1.0: genes detected
#' is the one that matters for complexity-based measures, and it is the one
#' that moves independently of total counts.
depth_by_age_within_sex <- function(obj, cfg, assay = "RNA",
                                    age_col = "age", sex_col = "sex",
                                    stage_col = "stage",
                                    age_levels = cfg$analysis$age_levels) {
  counts <- Seurat::GetAssayData(obj, assay = assay, layer = "counts")
  totals <- Matrix::colSums(counts)
  detected <- Matrix::colSums(counts > 0)

  age <- as.character(obj[[age_col]][, 1])
  sex <- as.character(obj[[sex_col]][, 1])
  stage <- if (stage_col %in% colnames(obj@meta.data))
    as.character(obj[[stage_col]][, 1]) else rep(NA_character_, ncol(obj))
  usable <- !is.na(age) & !is.na(sex) & age %in% age_levels

  summarise <- function(idx, stage_label) {
    by_age <- split(idx, factor(age[idx], levels = age_levels))
    by_age <- by_age[lengths(by_age) > 0]
    if (length(by_age) < 2) return(NULL)
    tbl <- data.frame(
      stage = stage_label,
      sex = sex[idx][1],
      age = names(by_age),
      n_cells = vapply(by_age, length, integer(1)),
      median_umi = vapply(by_age, function(i) stats::median(totals[i]), numeric(1)),
      median_genes = vapply(by_age, function(i) stats::median(detected[i]), numeric(1)),
      row.names = NULL, stringsAsFactors = FALSE)
    tbl$umi_ratio <- tbl$median_umi / tbl$median_umi[1]
    tbl$gene_ratio <- tbl$median_genes / tbl$median_genes[1]
    tbl$reference_age <- tbl$age[1]
    tbl
  }

  pooled <- lapply(split(which(usable), sex[usable]), summarise,
                   stage_label = "ALL STAGES (composition-confounded)")

  per_stage <- list()
  if (any(!is.na(stage))) {
    keep <- usable & !is.na(stage)
    for (st in sort(unique(stage[keep]))) {
      idx <- which(keep & stage == st)
      per_stage <- c(per_stage, lapply(split(idx, sex[idx]), summarise,
                                       stage_label = st))
    }
  }

  do.call(rbind, c(Filter(Negate(is.null), c(pooled, per_stage)),
                   list(make.row.names = FALSE)))
}

#' Report the depth-by-age table and say plainly what it implies.
report_depth_by_age <- function(tbl, cfg,
                                min_cells = cfg$gates$min_cells_per_stratum) {
  limit <- cfg$gates$max_depth_ratio %||% 1.3
  pooled_label <- "ALL STAGES (composition-confounded)"

  log_step("sequencing depth by age, within each sex ",
           "(each age is a separate hashtag in the same library):")
  print(tbl[, c("stage", "sex", "age", "n_cells", "median_umi", "median_genes",
                "umi_ratio", "gene_ratio")], row.names = FALSE)

  spread <- function(rows) {
    if (!nrow(rows)) return(NULL)
    by_group <- split(rows, paste(rows$stage, rows$sex))
    do.call(rbind, lapply(by_group, function(g) data.frame(
      stage = g$stage[1], sex = g$sex[1], n_min = min(g$n_cells),
      umi_spread = max(g$umi_ratio) / min(g$umi_ratio),
      gene_spread = max(g$gene_ratio) / min(g$gene_ratio),
      row.names = NULL)))
  }

  pooled <- spread(tbl[tbl$stage == pooled_label, ])
  log_step("  pooled over stages -- reported for completeness, NOT a verdict: ",
           "a stage-mix shift moves these on its own.")
  if (!is.null(pooled))
    for (i in seq_len(nrow(pooled)))
      log_step(sprintf("    %-7s UMIs %.2fx   genes detected %.2fx",
                       pooled$sex[i], pooled$umi_spread[i], pooled$gene_spread[i]))

  staged <- spread(tbl[tbl$stage != pooled_label, ])
  if (is.null(staged)) {
    log_step("  no per-stage rows (no stage column): the verdict cannot be given.")
    return(invisible(tbl))
  }

  # Strata too small to have a stable median say nothing either way.
  staged$evaluated <- staged$n_min >= min_cells
  log_step("  per stage, within each library -- THIS is the verdict:")
  for (i in seq_len(nrow(staged)))
    log_step(sprintf("    %-9s %-7s n>=%-4d UMIs %.2fx   genes %.2fx%s",
                     staged$stage[i], staged$sex[i], staged$n_min[i],
                     staged$umi_spread[i], staged$gene_spread[i],
                     if (!staged$evaluated[i]) "   (too few cells to judge)"
                     else if (max(staged$umi_spread[i], staged$gene_spread[i]) > limit)
                       "   <- over limit" else ""))

  judged <- staged[staged$evaluated, ]
  over <- judged[pmax(judged$umi_spread, judged$gene_spread) > limit, ]
  if (!nrow(judged)) {
    log_step("  every stage-by-sex stratum is below ", min_cells,
             " cells: no verdict.")
  } else if (nrow(over)) {
    log_step(sprintf("  WARNING: %d of %d evaluable strata exceed %.2fx: ",
                     nrow(over), nrow(judged), limit),
             paste(paste(over$stage, over$sex), collapse = ", "), ".")
    log_step("  Pseudotime position and CytoTRACE2 potency both track ",
             "transcriptional complexity, so an age trend in either could be ",
             "this trend. They do not corroborate each other while this holds.")
  } else {
    log_step(sprintf("  all %d evaluable strata are within %.2fx on both ",
                     nrow(judged), limit),
             "measures: the within-sex age trends are not explained by a ",
             "depth or complexity trend.")
  }
  invisible(tbl)
}

#' Which assay a step should read: a depth-matched one, or raw RNA.
#'
#' Two matched assays are built in step 3. RNAmatched equalises the AGES
#' within each sex and stage, which is what every within-sex age comparison
#' needs. RNAmatched_sex equalises every age x sex group within each stage,
#' which is what a comparison ACROSS the two libraries needs; the steps in
#' depth.cross_sex_steps read it. A step that expects matched counts and does
#' not find them says so rather than silently analysing unmatched ones.
matched_assay_for <- function(obj, cfg, step, assay_name = NULL) {
  wanted <- isTRUE(cfg$depth$match) &&
    step %in% unlist(cfg$depth$matched_steps %||% list())
  if (!wanted) {
    log_step("step ", step, " reads the RNA assay (not in depth.matched_steps)")
    return("RNA")
  }
  assay_name <- assay_name %||%
    if (step %in% unlist(cfg$depth$cross_sex_steps %||% list())) "RNAmatched_sex"
    else "RNAmatched"
  if (!assay_name %in% assay_names(obj))
    stop("step ", step, " is configured to use depth-matched counts, but '",
         assay_name, "' is not on this object.\n",
         "It is built in step 3; re-run step 3 to create it, or remove ",
         step, " from depth.matched_steps in config/config.yml.")
  log_step("step ", step, " reads the ", assay_name, " assay (depth-matched)")
  assay_name
}
