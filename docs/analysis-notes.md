# Reading the results

## What each method is asking

| Method | Question | Selection |
|---|---|---|
| `age_trend_clusters()` (step 6) | does expression change monotonically with age, and in what shape | `\|rho\| > age_trend.rho_cutoff` |
| `age_changing()` (step 6) | does any pairwise age comparison move more than chance | permutation-calibrated `\|log2FC\|` |
| `glm_lrt()` (step 7) | does the fitted NB model need an age term at all | permutation-calibrated p-value |
| `run_interaction_gsea()` (step 8) | do the two sexes' age trajectories differ, pathway-wise | fgsea `padj`, then a granule sensitivity run |
| `FindMarkers()` (step 4) | which genes differ between the sexes at one age | BH-adjusted Wilcoxon p |
| `composition_age_test()` (step 9) | does the *mix* of maturation stages shift with age | log-odds slope per stage, BH-adjusted |
| `compare_across_ages_by_sex()` (step 9) | do cells sit further along the trajectory | Wasserstein distance and median shift; KS p as a ranking |

## Pseudoreplication

Each age/sex group is one pooled library of several animals. Cells within a
group are not independent samples of the biology, so a per-cell p-value
answers "are these cells different" rather than "are these mice different",
and is far too small.

The permutation-calibrated methods (steps 6 and 7) shuffle the labels and
refit, which produces a null with the right correlation structure. Their
thresholds are trustworthy in a way BH-adjusted per-cell p-values are not.
Step 4's DEGs and the fgsea `padj` in step 8 are **not** calibrated this way:
treat them as rankings, not as evidence. `docs/original` sketches a sex-label
permutation for the GSEA that would fix this; it is not implemented here.

## Sex is confounded with batch

There is one 10x run per sex. Every sex difference in this project is a
sex-or-run difference, and nothing in the data can separate them. Age is
multiplexed *within* each run, so age comparisons within a sex are clean; the
sex x age interaction is the more defensible cross-sex statistic, because a
constant per-run offset cancels out of it. A plain male-vs-female contrast at
one age does not.

## Reading the sex x age interaction sign

`z_diff > 0`, and so `NES > 0`, means only: **the age trajectory is more
positive in males than in females**. It does not mean the gene or pathway is
higher in males, and it does not mean males age faster. A pathway can score
positive because it rises in males, or because it falls in females, or both.

`signal_type` and `dominant_pattern` in `gsea_sex_by_age_GOBP.csv` resolve
this by going back to the two per-sex rho values for every leading-edge gene.
Read those columns before writing a sentence about direction.

## Composition versus cell-intrinsic change

These are different findings and they are easy to confuse. "Gene X rises with
age in neutrophils" can mean the gene rises inside each cell, or that the
proportion of cells that express it highly has risen while no cell changed.

Everything in step 6 is stratified by stage, so it measures the first. Step 9
part A measures the second directly. Step 8 is stratified by stage for exactly
this reason: run on pooled stages, every gene that differs between stages
inherits an apparent age trend from a shift in the mix.

When a gene appears in both, the honest description is that both are moving.
When it appears only in the pooled analysis, it is composition.

## Cell numbers

Group sizes are very uneven — one arm has roughly ten times the cells of the
other at the middle timepoint. Consequences:

- Wilcoxon p-values scale with the smaller group's size, which is why the
  per-age GSEA ranks on `avg_log2FC` rather than on p.
- The interaction models are checked with the middle timepoint dropped
  (step 7 reports how many hits survive). A hit that only exists with it in is
  leverage, not biology.
- `scripts/06_age_trends.R` skips any stratum with fewer than 30 cells, and
  `age_changing()` warns when too few features pass the detection filter.

## Sparsity artefacts

A gene detected in 3 of 26 cells can produce a 256-fold "change". Two guards:
`detectable_genes()` requires detection in a minimum number *and* fraction of
cells in **every** age group, and `age_changing()` warns when the median
`|log2FC|` of its own hits exceeds 2, which is what a sparsity-driven result
looks like. If that warning fires, raise `permutation_de.min_cells`.

## Stage labels

`analysis.neutrophil_cluster_labels` maps cluster ids to maturation stages.
The mapping was made by eye from the module scores. Cluster numbering is not
stable across Seurat versions or after any change to step 3, so check
`results/figures/gmp_neutrophil_modules.pdf` against
`gmp_neutrophil_clusters.pdf` before trusting a stage-stratified result.

## ADT normalisation

The ADT `data` layer is `log1p(counts)` with each cell's median isotype-control
signal subtracted. That is the background-subtraction half of DSB: it removes
per-cell background but does not rescale by empty-droplet variance, so the
values are comparable across cells but are not DSB units. The ADT age-trend
analysis therefore filters on median signal rather than on a detection rate.
