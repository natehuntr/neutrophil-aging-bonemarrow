# ---------------------------------------------------------------------------
# Per-sample preprocessing: ADT/HTO assays, hashtag demultiplexing, QC
# metrics, doublet calls and filtering.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Ambient RNA correction.
#
# Reinstated after being removed earlier in this project's history. The
# argument for it here is specific rather than general: ambient contamination
# scales with library depth, so it is NOT independent of the sex confound, and
# the transcripts that bleed most between droplets are the highly expressed
# ones -- which in this compartment are exactly the granule genes that carry
# the maturation signal and top every differential list (S100a8/9, Retnlg,
# Mmp8/9, Camp, Ngp, Ltf).
#
# Correct first, then re-check whether the granule signal survives.
# ---------------------------------------------------------------------------

#' Estimate and remove ambient RNA with SoupX.
#'
#' SoupX needs a clustering to estimate contamination, so a throwaway object is
#' clustered here purely to supply one.
run_soupx <- function(mats, cfg) {
  require_packages("SoupX")
  if (is.null(mats$raw_gex))
    stop("ambient correction needs the raw matrix; read the sample with need_raw = TRUE")

  sc <- SoupX::SoupChannel(tod = mats$raw_gex, toc = mats$gex,
                           channelName = mats$sample_id)

  quick <- Seurat::CreateSeuratObject(counts = mats$gex)
  quick <- Seurat::NormalizeData(quick, verbose = FALSE)
  quick <- Seurat::FindVariableFeatures(quick, verbose = FALSE)
  quick <- Seurat::ScaleData(quick, verbose = FALSE)
  quick <- Seurat::RunPCA(quick, npcs = cfg$ambient$quick_cluster_dims, verbose = FALSE)
  quick <- Seurat::FindNeighbors(quick, dims = seq_len(cfg$ambient$quick_cluster_dims),
                                 verbose = FALSE)
  quick <- Seurat::FindClusters(quick, resolution = cfg$ambient$quick_cluster_resolution,
                                verbose = FALSE)

  sc <- SoupX::setClusters(sc, stats::setNames(as.character(quick$seurat_clusters),
                                               colnames(quick)))
  sc <- SoupX::autoEstCont(sc)

  rho <- sc$fit$rhoEst
  log_step(sprintf("  SoupX contamination estimate: %.1f%% of counts are ambient",
                   100 * rho))
  top <- utils::head(sc$soupProfile[order(-sc$soupProfile$est), ], 20)
  log_step("  top ambient genes: ", paste(utils::head(rownames(top), 10), collapse = ", "))

  list(counts = SoupX::adjustCounts(sc, roundToInt = TRUE),
       rho = rho, top_genes = top)
}

#' Build the RNA Seurat object, ambient-corrected unless switched off.
#'
#' The contamination fraction is recorded on the object: it is per library, so
#' a difference between the sexes is another face of the depth confound and
#' belongs in the diagnostics rather than being discarded.
create_rna_object <- function(mats, cfg) {
  method <- cfg$ambient$method %||% "none"

  if (identical(method, "none")) {
    log_step("  ambient correction: off (ambient.method)")
    obj <- Seurat::CreateSeuratObject(counts = mats$gex, project = mats$sample_id)
    obj$ambient_rho <- NA_real_
    return(obj)
  }

  soup <- switch(method,
    soupx = run_soupx(mats, cfg),
    stop("unknown ambient.method: ", method))

  obj <- Seurat::CreateSeuratObject(counts = soup$counts, project = mats$sample_id)
  obj$ambient_rho <- soup$rho
  attr(obj, "ambient_top_genes") <- soup$top_genes
  obj
}

#' Isotype-centred ADT normalisation.
#'
#' log1p of the raw ADT counts, then subtract each cell's median isotype
#' control signal. This is the background-subtraction half of DSB: it removes
#' per-cell background but does not rescale by empty-droplet variance, so
#' values are comparable across cells but are not in DSB units.
#'
#' The isotype rows themselves stay in the matrix and end up centred near
#' zero, which is what makes them useful as a sanity check on the panel.
normalise_adt_isotype <- function(mats, cells, features, cfg) {
  isotype_features <- grep(cfg$adt$isotype_pattern, features, value = TRUE)
  if (length(isotype_features) == 0)
    stop("no isotype controls matched '", cfg$adt$isotype_pattern, "' in the ADT panel")

  adt_log <- log1p(as.matrix(mats$adt[features, cells, drop = FALSE]))
  isotype_median <- matrixStats::colMedians(adt_log[isotype_features, , drop = FALSE])
  sweep(adt_log, 2, isotype_median, "-")
}

#' Antibody rows that are not hashtags.
adt_feature_names <- function(mats, hashtags) {
  setdiff(rownames(mats$adt), hashtags)
}

#' Attach the ADT and HTO assays to an RNA object.
#'
#' The normalised matrix is written straight into the ADT `data` layer: it is
#' already on a log scale, so Seurat's own NormalizeData must not be run on
#' this assay afterwards.
add_protein_assays <- function(obj, mats, cfg, hashtags) {
  cells <- colnames(obj)
  missing_htos <- setdiff(hashtags, rownames(mats$adt))
  if (length(missing_htos))
    stop("hashtags named in config are absent from the ADT panel: ",
         paste(missing_htos, collapse = ", "),
         "\nPanel rows are: ", paste(rownames(mats$adt), collapse = ", "))

  adt_features <- adt_feature_names(mats, hashtags)

  obj[["ADT"]] <- Seurat::CreateAssayObject(
    counts = mats$adt[adt_features, cells, drop = FALSE])
  obj[["ADT"]]$data <- normalise_adt_isotype(mats, cells, adt_features, cfg)

  obj[["HTO"]] <- Seurat::CreateAssayObject(
    counts = mats$adt[hashtags, cells, drop = FALSE])
  obj
}

#' Demultiplex hashtags and translate the calls into age / age_sex metadata.
demultiplex_hashtags <- function(obj, sample_cfg) {
  obj <- Seurat::NormalizeData(obj, assay = "HTO",
                               normalization.method = "CLR", margin = 2)
  obj <- Seurat::MULTIseqDemux(obj, assay = "HTO", autoThresh = TRUE)

  # Factor levels here would silently reorder the age mapping below.
  obj$MULTI_ID <- as.character(obj$MULTI_ID)
  print(table(obj$MULTI_ID))

  hto_to_age <- unlist(sample_cfg$hashtags)
  obj$age <- unname(hto_to_age[match(obj$MULTI_ID, names(hto_to_age))])
  obj$sex <- sample_cfg$sex
  obj$age_sex <- ifelse(is.na(obj$age), NA_character_,
                        paste0(obj$age, "_", toupper(substr(sample_cfg$sex, 1, 1))))
  obj
}

#' Standard QC metrics: complexity, mitochondrial, ribosomal and haemoglobin
#' content. Ratios are 0-1; percent.hb is 0-100, matching the QC thresholds.
add_qc_metrics <- function(obj, cfg) {
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  obj$mitoRatio <- Seurat::PercentageFeatureSet(obj, pattern = cfg$qc$mito_pattern) / 100
  obj$riboRatio <- Seurat::PercentageFeatureSet(obj, pattern = cfg$qc$ribo_pattern) / 100

  hb_genes <- grep(cfg$qc$hb_pattern, rownames(obj), value = TRUE)
  obj[["percent.hb"]] <- if (length(hb_genes)) {
    Seurat::PercentageFeatureSet(obj, features = hb_genes)
  } else {
    0
  }
  obj
}

#' Call doublets per hashtag group with scDblFinder.
#'
#' The SCE is built straight from the counts: as.SingleCellExperiment() wants
#' a populated `data` layer, which does not exist this early in the pipeline.
#'
#' `samples` is the hashtag call, so doublet rates are estimated per
#' multiplexed group. Cells called Doublet or Negative by the demultiplexer
#' form their own pseudo-groups; they are dropped by filter_cells() regardless
#' of what scDblFinder decides about them.
add_doublet_calls <- function(obj) {
  # scDblFinder runs `counts(sce) <- ...` internally, and that method needs
  # `assay<-` resolvable through the search path. Checking here turns a
  # confusing error raised deep inside another package into an actionable one,
  # and costs nothing.
  if (!exists("assay<-", mode = "function"))
    stop("SummarizedExperiment is not attached, so scDblFinder will fail with\n",
         '  Error in assay(object, "counts") <- value : ',
         'could not find function "assay<-"\n',
         "Call require_packages(\"SummarizedExperiment\") before this step.")

  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")))
  sce <- scDblFinder::scDblFinder(sce, samples = obj$MULTI_ID)

  obj$scDblFinder.class <- sce$scDblFinder.class
  obj$scDblFinder.score <- sce$scDblFinder.score

  log_step(sprintf("%d doublets detected (%.1f%% of cells)",
                   sum(obj$scDblFinder.class == "doublet"),
                   mean(obj$scDblFinder.class == "doublet") * 100))
  obj
}

#' Apply QC thresholds, computed per library.
#'
#' A global nFeature/nCount cutoff is the wrong instrument here. Mature
#' neutrophils carry the least RNA in this compartment, so a fixed threshold
#' removes them preferentially -- and it removes more of them from the
#' shallower library. That converts a depth difference between the sexes into
#' an apparent composition difference, which is the most likely explanation for
#' a progenitor-skewed male sample.
#'
#' scuttle's MAD approach sets the threshold from each library's own
#' distribution, so "outlier" means outlying for that library rather than
#' relative to the deeper one.
filter_cells <- function(obj, cfg, hashtags, batch_col = "sex") {
  require_metadata(obj, c("age", "MULTI_ID", "scDblFinder.class",
                          "log10GenesPerUMI", "mitoRatio", "percent.hb"),
                   context = "QC filtering")
  n_before <- ncol(obj)

  # Demultiplexing and doublets first: these are not depth-driven, and the MAD
  # thresholds should be computed on real single cells.
  obj <- select_cells(obj, list(
    "hashtag call is a real sample" = obj$MULTI_ID %in% hashtags,
    "age was assigned"              = !is.na(obj$age),
    "singlet"                       = obj$scDblFinder.class == "singlet"
  ), context = "demultiplexing and doublets")

  discard <- if (identical(cfg$qc$method, "mad")) {
    mad_outliers(obj, cfg, batch_col)
  } else {
    log_step("  QC method 'fixed': global thresholds (see qc.method in the config)")
    !(obj$log10GenesPerUMI > cfg$qc$min_log10_genes_per_umi)
  }

  obj <- select_cells(obj, list(
    "not a per-library count/feature outlier" = !discard,
    "mitoRatio below cutoff"  = obj$mitoRatio < cfg$qc$max_mito_ratio,
    "percent.hb below cutoff" = obj$percent.hb < cfg$qc$max_percent_hb
  ), context = "QC thresholds")

  log_step(sprintf("cell filtering: %d -> %d cells retained (%.1f%% kept)",
                   n_before, ncol(obj), ncol(obj) / n_before * 100))
  obj
}

#' Per-library median-absolute-deviation outlier calls.
#'
#' `batch` is what makes this per-library: scuttle computes the thresholds
#' separately within each batch, so a shallower library is judged against
#' itself.
mad_outliers <- function(obj, cfg, batch_col = "sex") {
  require_packages("scuttle", "SingleCellExperiment")

  counts <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  qc <- scuttle::perCellQCMetrics(sce)

  batch <- if (batch_col %in% colnames(obj@meta.data))
    as.character(obj[[batch_col]][, 1]) else NULL
  if (is.null(batch))
    log_step("  no ", batch_col, " column; MAD thresholds computed on all cells at once")

  filters <- scuttle::perCellQCFilters(qc, batch = batch, nmads = cfg$qc$mad_nmads)
  discard <- filters$discard

  # Report the thresholds actually used, per library. If they differ a lot
  # between the sexes, that difference IS the confound, made visible.
  for (b in unique(batch %||% "all")) {
    idx <- if (is.null(batch)) rep(TRUE, ncol(obj)) else batch == b
    log_step(sprintf("  %s: %d of %d cells discarded (%.1f%%); sum >= %.0f, detected >= %.0f",
                     b, sum(discard[idx]), sum(idx), 100 * mean(discard[idx]),
                     min(qc$sum[idx & !discard]), min(qc$detected[idx & !discard])))
  }
  discard
}

#' Retained-cell fraction by group, before and after filtering.
#'
#' The table that decides whether a composition finding survives interpretation.
#' If the shallower library loses disproportionately many of one stage, the
#' composition difference is a QC artefact and this is the number that shows it.
retention_table <- function(before, after, group_cols = c("sex", "stage")) {
  present <- intersect(group_cols, intersect(colnames(before@meta.data),
                                             colnames(after@meta.data)))
  if (!length(present)) stop("none of ", paste(group_cols, collapse = ", "), " are present")

  tabulate_by <- function(obj) {
    keys <- lapply(present, function(col) as.character(obj[[col]][, 1]))
    names(keys) <- present
    as.data.frame(do.call(table, keys), responseName = "n", stringsAsFactors = FALSE)
  }

  n_before <- tabulate_by(before)
  n_after <- tabulate_by(after)
  names(n_before)[names(n_before) == "n"] <- "n_before"
  names(n_after)[names(n_after) == "n"] <- "n_after"

  out <- merge(n_before, n_after, by = present, all = TRUE)
  out$n_before[is.na(out$n_before)] <- 0
  out$n_after[is.na(out$n_after)] <- 0
  out$retained_fraction <- ifelse(out$n_before > 0, out$n_after / out$n_before, NA_real_)
  out[order(out[[present[1]]], -out$n_before), ]
}

#' Does retention differ enough between libraries to explain a composition
#' difference on its own?
retention_asymmetry <- function(retention, sex_col = "sex", stage_col = "stage") {
  if (!all(c(sex_col, stage_col) %in% names(retention))) return(NULL)
  wide <- stats::reshape(retention[, c(sex_col, stage_col, "retained_fraction")],
                         idvar = stage_col, timevar = sex_col, direction = "wide")
  names(wide) <- sub("^retained_fraction\\.", "retained_", names(wide))
  fractions <- wide[, setdiff(names(wide), stage_col), drop = FALSE]
  wide$retention_ratio <- apply(fractions, 1, function(x)
    if (any(is.na(x)) || min(x) == 0) NA_real_ else max(x) / min(x))
  wide[order(-wide$retention_ratio), ]
}
