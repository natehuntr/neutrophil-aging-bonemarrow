# What changed from `docs/original/bone_marrowaging.qmd`

The original was a 5,799-line Quarto notebook. The analysis it describes is
preserved; what changed is structure, and a small number of things that were
actually broken.

## Structure

**Three copies of the per-sample pipeline became one.** The notebook contained
the same ~800-line block three times: "Male bone marrow", "Female bone marrow
(v2)" and "Female bone marrow". Normalising object names away, the male and
female-v2 blocks differed only in the sample ID, the hashtag-to-age suffix, and
nine tuning numbers (PCA dimensions and `use_partition`). Those nine numbers
now live in `config/config.yml` under each sample's `dims`, and the code exists
once in `R/` + `scripts/01_preprocess.R` and `02_annotate.R`.

**The third block was dead.** "Female bone marrow" (original lines 1614-2134)
reprocessed the *same* sample as "Female bone marrow (v2)"
(`MCC8209A4_Female_BM_4`), was shorter, lacked the trajectory analysis, and its
objects were overwritten by the v2 block. It is not carried over. If it was
meant to be a different sample, add it to `config.yml` and it will run through
the same code path.

**Copy-pasted chunks became loops or functions.** Per-age trajectory blocks
(6 copies), the GO ORA blocks (3 identical copies differing only in the age),
the pseudobulk PCA blocks (4 copies, before a function doing the same thing
appeared later in the notebook), and the per-stage `age_trend_clusters` blocks
(6 copies) are now single parameterised calls.

**Paths and outputs.** Every `/Users/costahn/Desktop/...` path and the one
`setwd()` are gone. Inputs come from `paths.cellranger_root`; outputs go to
`results/{objects,tables,figures}` through helpers that log where they wrote.
Output filenames now encode stage, sex and assay consistently, so a table's
name says what is in it.

## Genuine bugs fixed

1. **SoupX was computed and discarded.** `soup_male_bm` was built from the
   corrected counts, but the object carried forward (`male_bm`) was built from
   the uncorrected `filtered_gex`, so the correction reached nothing. Ambient
   correction has since been removed from the pipeline entirely, which makes
   the original behaviour the actual behaviour. The raw (empty-droplet) matrix
   is no longer read at all, since SoupX was its only consumer.

2. **`main_label_score` was a matrix.** `obj$main_label_score <-
   cell_identity$scores` assigns a cells-x-labels matrix into a metadata column.
   It now stores the winning score, `apply(scores, 1, max)`.

3. **The ADT age-trend functions silently overwrote the RNA ones.**
   `run_pairwise()`, `calibrate()` and `classify()` were each defined twice at
   the top level, the second (ADT) definition clobbering the first (RNA). Any
   RNA analysis re-run after that point in the session would have used the ADT
   code. There is now one `age_changing(..., assay = )` covering both.

4. **`age_trend_clusters()` was defined twice**, the second version adding the
   GSEA ranking vector and NA handling. Only the second is kept; it is a strict
   superset.

5. **`bm_merge` used `CytoTRACE2_Potency`, which it did not have.** Potency was
   computed on the `neustem_*` subsets, but the merge was built from the full
   `male_bm`/`female_bm` objects, so the neutrophil subset filter was testing an
   absent column. Step 3 now transfers the potency calls onto the merged object
   by barcode; cells outside the scored subset are `NA` and are excluded
   explicitly rather than by accident.

6. **Undefined variables in the GLM sections.** `cd$pt` was used but never
   assigned; `cds` was used where `combined_cds` was meant; `f_lean` appeared in
   a `stopifnot()` above its own definition; `mat_neus` was used one chunk
   before it was created; `male_GAM` was plotted before it was fitted; a male
   plot coloured points by `neustem_3m_female_cds`. Step 7 defines pseudotime
   from `monocle3::pseudotime()` and orders the definitions correctly.

7. **A stray `Default` on its own line** (original line ~3078) would have thrown
   an object-not-found error mid-chunk. Removed.

8. **`AverageExpression` column names were assumed.** Several chunks hard-coded
   `c("g3m", "g12m", "g18m")` and would silently reorder or fail if Seurat's
   name-mangling changed. Columns are now matched on the trailing age string.

9. **`max(abs(pct_diff))` on an all-NA group returned `-Inf`** with a warning;
   it now returns `NA`.

## Behaviour changes worth knowing about

- **Trajectory roots are chosen programmatically.** `order_cells()` and
  `choose_graph_segments()` without a root open a Shiny window, so the original
  could not run unattended. `learn_trajectory()` now takes `root_group` and
  picks the principal-graph node closest to the GMP cells. Roots picked by hand
  in the original may differ, which shifts pseudotime values; pass
  `root_group = NULL` interactively to reproduce the manual choice.
  `find_traj_genes()`, which called `choose_graph_segments()` interactively to
  trim the graph and then re-learned it, is not carried over: its cell
  selection was not recorded, so it was not reproducible.

- **Per-sample trajectories dropped.** The notebook ran per-age trajectories
  twice: once per sample on the `neustem_*` objects, and again on the merged
  GMP+neutrophil object. Only the merged version is kept — it is the later one
  and the one whose outputs the downstream sections read.

- **9m handling is explicit.** The notebook dropped 9m by writing
  `age != "9m"` in about fifteen places. `analysis.age_levels` now holds the
  three ages used for comparisons, and `all_age_levels` the four that exist.

## Not carried over

- Exploratory `plotSmoothers()` calls for individual genes (roughly 40 of them,
  one gene at a time). Fit the GAMs with step 5 and call `plotSmoothers()` on
  the saved objects for whatever genes you want.
- `print()`/`head()` calls used to eyeball an object mid-chunk.
- The empty chunks and the `## Stem cells and neutrophils` heading with no code.
- The `library()` calls scattered through 30 chunks: dependencies are declared
  once in `R/packages.R` and installed by `scripts/00_install_dependencies.R`.

## Second review pass

A further read of the restructured code found these. The first four were
present in the original notebook too; the rest are robustness gaps that the
original avoided only by being run by hand, one chunk at a time.

- **Spline predictions used the wrong knots.** The `pt_age_sex_fitted.csv`
  grid was built with `model.matrix(f_lean, newdat)`, which re-derives the
  `ns()` knots *from the prediction grid* instead of reusing the ones the model
  was fitted with. The two bases differ (up to ~0.09 in basis value on a
  realistic grid), so every fitted value in that table was wrong. Step 7 now
  builds the grid from the fitted model's `terms()`, whose `predvars` attribute
  carries the original knots, and asserts the design columns match the fitted
  coefficients.

- **`padj < 0.05` on fgsea output.** fgsea leaves `padj` as NA for pathways it
  cannot score, so the filter produced NA rows, which `collapsePathways()`
  rejects. Filters now go through `which()`.

- **Per-gene `cor.test()` in a loop.** The age-trend, pseudotime and per-sex
  rho functions each called `apply(expr, 1, cor.test)`, which densifies the
  whole matrix and runs tens of thousands of tests. `spearman_rows()` ranks
  once and calls `cor()`, then computes the p-value from the same asymptotic t
  approximation `cor.test(method = "spearman", exact = FALSE)` uses — verified
  identical to ten decimal places, including the NA for a zero-variance gene.
  This function is called twelve times per run, so it was a large share of the
  total runtime.

- **`min(x, na.rm = TRUE)` over an all-NA row** returned `Inf` with a warning
  rather than reporting "not significant at any age".

- **PCA and t-SNE on small strata.** `RunPCA` defaults to 50 components, which
  errors on a subset with fewer cells — or, for the ~30-antibody ADT panel,
  fewer features. `safe_npcs()` and `safe_dims()` cap both, and t-SNE
  perplexity is capped at `(n - 2) / 3`.

- **Assay conversion before normalisation.** `as.SingleCellExperiment()` wants
  a populated `data` layer, which does not exist when doublets are called. The
  SCE is now built directly from the counts.

- **Non-finite QC metrics.** A cell with a single UMI gives
  `log10GenesPerUMI = Inf/NaN`, making the filter vector NA and breaking the
  subset. Non-TRUE now counts as a failure, and the count is logged.

- **Split layers after a merge.** `JoinLayers()` only joins the default assay,
  so the ADT assay stayed split and `ScaleData` on it did not do what it
  appeared to. `join_all_layers()` joins every v5 assay.

- **Missing groups aborted whole scripts.** `FindMarkers` on an absent or
  3-cell ident, and `fitGAM` on an empty gene set, both failed with opaque
  errors. Those cases are now checked, skipped or explained.

- **`FindMarkers(layer = "data")`** was removed: Seurat builds disagree on
  whether that argument is `slot` or `layer`, and "data" is the default in both.

- **`glm_group_means()` assumed exactly three ages.** It indexed `Beta[, 2]`
  and `Beta[, 3]` directly; it now derives the means from however many levels
  the design has, and errors if the design is not `~ age`.
