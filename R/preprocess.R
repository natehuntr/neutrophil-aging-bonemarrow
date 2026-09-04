# ---------------------------------------------------------------------------
# Per-sample preprocessing: ADT/HTO assays, hashtag demultiplexing, QC
# metrics, doublet calls and filtering.
# ---------------------------------------------------------------------------

#' Build the RNA Seurat object for a sample.
create_rna_object <- function(mats) {
  Seurat::CreateSeuratObject(counts = mats$gex, project = mats$sample_id)
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

#' Apply the QC thresholds from the config.
#'
#' Cells whose hashtag call was Doublet or Negative have no age and are
#' dropped here along with the QC failures.
filter_cells <- function(obj, cfg, hashtags) {
  qc <- cfg$qc
  require_metadata(obj, c("age", "MULTI_ID", "scDblFinder.class",
                          "log10GenesPerUMI", "mitoRatio", "percent.hb"),
                   context = "QC filtering")

  n_before <- ncol(obj)
  # NA counts as a failure inside select_cells(), which also covers the cell
  # with a single UMI whose log10GenesPerUMI is Inf/NaN.
  obj <- select_cells(obj, list(
    "hashtag call is a real sample" = obj$MULTI_ID %in% hashtags,
    "age was assigned"              = !is.na(obj$age),
    "singlet"                       = obj$scDblFinder.class == "singlet",
    "log10GenesPerUMI above cutoff" = obj$log10GenesPerUMI > qc$min_log10_genes_per_umi,
    "mitoRatio below cutoff"        = obj$mitoRatio < qc$max_mito_ratio,
    "percent.hb below cutoff"       = obj$percent.hb < qc$max_percent_hb
  ), context = "QC filtering")

  log_step(sprintf("cell filtering: %d -> %d cells retained (%.1f%% kept)",
                   n_before, ncol(obj), ncol(obj) / n_before * 100))
  obj
}
