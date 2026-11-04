# Embryonic neutrophils (mouse prenatal atlas)

`embryonic_neutrophils.qmd` is the reproducible version of the analysis behind
*"Beyond immunity: Mid-gestation neutrophils prioritise homeostasis over immune
function in the mouse embryo"*. It is separate from the bone-marrow ageing
pipeline in the rest of this repository.

## Running it

1. Download the atlas files from the
   [JAX atlas download page](https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax/download/mtx/Example_Script)
   into `embryo/data/`: `df_cell.rds`, `gene_annotation.csv.gz`, and
   `gene_count.run_<N>.mtx.gz` + `cell_annotation.run_<N>.csv.gz` for runs 4,
   13–16, 18–29, 17_sub1 and 17_sub2.
2. Install the packages, then lock their versions:
   ```bash
   Rscript install_dependencies.R
   Rscript -e 'renv::snapshot(type = "all")'   # commit the renv.lock this writes
   ```
3. Render:
   ```bash
   quarto render embryonic_neutrophils.qmd
   ```

Every expensive step is saved in `results/objects/` and reloaded on the next
render; use `-P force:true` after changing an early setting. Figures go to
`results/figures/` (PNG + PDF, named after their chunk) and tables to
`results/tables/`. All analysis settings are in the `settings` chunk at the top.

## What changed from the original notebook, and why

### It did not run top to bottom

- Most chunks were `eval: false` and the rest depended on files made by hand
  or by code that was not in the document (`neustem_subset_v3.rds` with
  `cytotrace_potency`, `e14_neus_allcells.rds`, `p0_neus_allcells.rds`,
  `e14_cellchat.rds`, `e18_cellchat.rds`, `neutrophils_only.rds`, an
  `adult_neu_subset`). Every object is now made in the document. The CytoTRACE2
  step was missing and has been added. The all-cell E14/P0 objects for CellChat
  are now built from the atlas (see "Check these" below).
- QC started with `readRDS("clean_rds/unsubset_object.rds")`, a file that was
  only written *after* QC. The steps now pass objects forward.
- Undefined objects: `neu_subset` was plotted before it was made, and `rats`,
  `gg2`, `cellchat` (as a merged object) and `adult_neu_subset` were never
  defined. A stray `d)` argument in a `DimPlot` call was a syntax error.
- `FindNeighbors()` used the default reduction `"pca"`, but the PCA was named
  `"rna_pca"`, so it would have failed.

### Steps that need a person clicking

- `trace(createCellChat, edit = TRUE)` hand-edited a package function each
  session. CellChat now gets the expression matrix and metadata directly,
  which avoids the Seurat v5 problem the edit was working around.
- Islands were removed by lasso in a Shiny app. There is now an automatic rule
  (DBSCAN on the UMAP) or a saved cell list (`params$retained_cells_file`); the
  Shiny selector is kept, and writes that file.
- monocle3's `choose_graph_segments()` and `order_cells()` were interactive.
  The path start and end, and the pseudotime root, are now chosen from
  CytoTRACE2 potency.

### Errors in the analysis

- **Cell cycle scoring used human genes.** `cc.genes` is uppercase human
  symbols (`MCM5`), which match almost nothing in mouse data (`Mcm5`), so the
  S/G2M scores, and the `CC.Difference` regressed out of every scaling, were
  close to meaningless. They are now converted to mouse symbols.
- **Wrong log base.** scran's `logNormCounts()` gives log2 values; Seurat
  assumes natural log, so `avg_log2FC` from `FindMarkers`/`FindAllMarkers` was
  off. The scran-normalised values are now stored as `log1p`.
- **Assays out of step.** After decontX became the default assay, variable
  features and scaling ran on decontX's *un-normalised* counts and PCA then ran
  on the RNA assay. Everything now runs on normalised decontX data, and the
  pseudobulk uses the raw RNA counts (DESeq2 needs integers).
- **decontX after normalisation.** Ambient correction now comes before
  normalisation, so the normalised values are of corrected counts.
- **Pseudobulk samples by position.** Embryos were removed by column number
  (`[, -c(1, 2, 3, 12, ...)]`) and given a day range by position in a typed
  vector. If anything changed the column order, samples would be silently
  mislabelled. Groups now come from each embryo's own `Embryonic_day`.
- **Clustering on the wrong graph.** Leiden was run on `RNA_nn` (the kNN graph)
  rather than the shared-nearest-neighbour graph it is designed for.
- The neutrophil re-analysis ran `ScaleData` before `FindVariableFeatures`
  (so it scaled the parent object's genes), used different numbers of PCs for
  the neighbours and the UMAP, and plotted stale clusters before re-clustering.
- `FindAllMarkers(group.by = ...)` is ignored; identities are now set
  explicitly.
- Day labels sorted alphabetically (E10 before E8); they are now ordered by
  developmental time.
- The text called the data 10x Chromium; it is sci-RNA-seq3.

### Reproducibility

- Random sampling, scDblFinder, decontX, scran clustering, UMAP, Leiden,
  monocle3 and CellChat all use the one seed in the settings.
- Ensembl ID → symbol conversion used a live biomaRt query, whose answer
  changes with each Ensembl release. It now uses the annotation shipped with
  the atlas.
- `renv::snapshot()` was called from inside the document. Versions are now
  recorded once, in a committed `renv.lock`, and the document ends with
  `sessionInfo()`.
- Paths are no longer scattered (`clean_rds/`, `pseudobulk_res/`, the working
  directory); figures previously pasted in as PNGs are now drawn by the code.
- Packages used but never loaded (scran, celda, DESeq2, CellChat, monocle3,
  SeuratWrappers, pheatmap) are now checked and loaded at the start; unused ones
  (harmony, garnett, GEOquery, AnnotationDbi, org.Mm.eg.db, UCell) are gone.
- The "run through clusterProfiler" step was described but not done; it is now
  an over-representation test of GO biological processes (fgsea + msigdbr),
  against all tested genes as background.

## Check these before trusting a re-run

- **Clustering resolution.** The text said 0.8 was best, but the code used 0.2.
  The settings keep 0.2 (what the code did); change `cfg$resolution` if 0.8 was
  intended.
- **Pseudobulk clusters.** `cfg$pseudobulk_clusters = c("2", "5")`. Cluster
  numbers change between software versions, so check the marker dot plot in
  the Pseudobulking section after every re-run.
- **Doublet sample column.** The original used `mouse.id`, which is not an
  atlas column I could identify. The settings use `RT_group`; change
  `cfg$doublet_sample_col` if you meant something else (e.g. `embryo_id`).
- **CellChat inputs.** `e14_neus_allcells.rds` and `p0_neus_allcells.rds` are
  rebuilt from the atlas metadata. That is my reconstruction of how they were
  made. Also, your plots referred to groups called "CNS neurons",
  "Sensory/PNS neurons" and "Neural progenitors", but the merge map in the
  notebook puts all of them into "Neurons". The plots now use the groups the
  map produces. If your final figures used a finer map, replace
  `celltype_merge_map`.
- **Island removal.** The DBSCAN settings (`cfg$islands`) will need tuning by
  eye the first time. Once you are happy, copy
  `results/tables/retained_cells.txt` into the repository and set
  `params$retained_cells_file` to it, so the selection is fixed.
- **Removed exploratory plots.** The PCA plot of the whole white-blood-cell
  object was dropped. It needed a slow regression on every cell and was not
  used afterwards.
