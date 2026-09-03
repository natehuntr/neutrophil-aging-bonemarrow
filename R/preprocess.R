# ---------------------------------------------------------------------------
# Per-sample preprocessing: ambient RNA, ADT/HTO assays, demultiplexing,
# QC metrics, doublets, filtering.
# ---------------------------------------------------------------------------

#' Estimate and remove ambient ("soup") RNA with SoupX.
#'
#' SoupX needs a clustering to estimate contamination, so a throwaway Seurat
#' object is clustered here purely to supply one.
run_soupx <- function(mats, cfg) {
  sc <- SoupX::SoupChannel(tod = mats$raw_gex,
                           toc = mats$filtered_gex,
                           channelName = mats$sample_id)

  quick <- Seurat::CreateSeuratObject(counts = mats$filtered_gex)
  quick <- Seurat::NormalizeData(quick, verbose = FALSE)
  quick <- Seurat::FindVariableFeatures(quick, verbose = FALSE)
  quick <- Seurat::ScaleData(quick, verbose = FALSE)
  quick <- Seurat::RunPCA(quick, npcs = cfg$soupx$quick_cluster_dims, verbose = FALSE)
  quick <- Seurat::FindNeighbors(quick, dims = seq_len(cfg$soupx$quick_cluster_dims),
                                 verbose = FALSE)
  quick <- Seurat::FindClusters(quick, resolution = cfg$soupx$quick_cluster_resolution,
                                verbose = FALSE)

  clusters <- stats::setNames(as.character(quick$seurat_clusters), colnames(quick))
  sc <- SoupX::setClusters(sc, clusters)
  sc <- SoupX::autoEstCont(sc)

  log_step(sprintf("SoupX contamination estimate: %.3f", sc$fit$rhoEst))
  top_soup <- utils::head(sc$soupProfile[order(-sc$soupProfile$est), ], 20)

  list(
    channel        = sc,
    rho            = sc$fit$rhoEst,
    top_soup_genes = top_soup,
    counts         = SoupX::adjustCounts(sc, roundToInt = TRUE)
  )
}

#' Build the RNA Seurat object for a sample.
#'
#' `use_corrected` selects between the SoupX-adjusted counts and the raw
#' filtered counts. See the note in config/config.yml: the original notebook
#' used the uncorrected counts.
create_rna_object <- function(mats, soup, cfg) {
  counts <- if (isTRUE(cfg$soupx$use_corrected_counts) && !is.null(soup)) {
    log_step("using SoupX-corrected counts")
    soup$counts
  } else {
    log_step("using uncorrected filtered counts")
    mats$filtered_gex
  }
  Seurat::CreateSeuratObject(counts = counts, project = mats$sample_id)
}

#' Isotype-centred ADT normalisation.
#'
#' log1p of the raw ADT counts, then subtract each cell's median isotype
#' control signal. This is the background-subtraction half of DSB: it removes
#' per-cell background but does not rescale by the empty-droplet variance.
normalise_adt_isotype <- function(mats, cells, cfg) {
  isotype_features <- grep(cfg$adt$isotype_pattern, rownames(mats$raw_adt), value = TRUE)
  if (length(isotype_features) == 0)
    stop("no isotype controls matched '", cfg$adt$isotype_pattern, "' in the ADT panel")

  adt_features <- adt_feature_names(mats, cfg)
  cell_adt <- as.matrix(mats$filtered_adt[adt_features, cells, drop = FALSE])

  adt_log <- log1p(cell_adt)
  isotype_median <- matrixStats::colMedians(adt_log[isotype_features, , drop = FALSE])
  sweep(adt_log, 2, isotype_median, "-")
}

#' Antibody rows that are not hashtags.
adt_feature_names <- function(mats, cfg, hashtags = NULL) {
  hashtags <- hashtags %||% grep("Hashtag", rownames(mats$raw_adt), value = TRUE)
  setdiff(rownames(mats$raw_adt), hashtags)
}

#' Attach the ADT and HTO assays to an RNA object.
add_protein_assays <- function(obj, mats, cfg, hashtags) {
  cells <- colnames(obj)
  missing_htos <- setdiff(hashtags, rownames(mats$filtered_adt))
  if (length(missing_htos))
    stop("hashtags named in config are absent from the ADT panel: ",
         paste(missing_htos, collapse = ", "))

  adt_features <- adt_feature_names(mats, cfg, hashtags)

  obj[["ADT"]] <- Seurat::CreateAssayObject(
    counts = mats$filtered_adt[adt_features, cells, drop = FALSE])
  obj[["ADT"]]$data <- normalise_adt_isotype(mats, cells, cfg)[adt_features, , drop = FALSE]

  obj[["HTO"]] <- Seurat::CreateAssayObject(
    counts = mats$filtered_adt[hashtags, cells, drop = FALSE])
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
add_doublet_calls <- function(obj) {
  sce <- Seurat::as.SingleCellExperiment(obj)
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
  keep <- !is.na(obj$age) &
    obj$MULTI_ID %in% hashtags &
    obj$scDblFinder.class == "singlet" &
    obj$log10GenesPerUMI > qc$min_log10_genes_per_umi &
    obj$mitoRatio < qc$max_mito_ratio &
    obj$percent.hb < qc$max_percent_hb

  n_before <- ncol(obj)
  obj <- subset(obj, cells = colnames(obj)[keep])
  log_step(sprintf("cell filtering: %d -> %d cells retained (%.1f%% kept)",
                   n_before, ncol(obj), ncol(obj) / n_before * 100))
  obj
}
