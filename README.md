# Neutrophil ageing in mouse bone marrow

CITE-seq analysis of bone marrow granulopoiesis across the mouse lifespan
(3, 9, 12 and 18 months), one hashed 10x run per sex, with a paired antibody
panel (ADT) and hashtag oligos (HTO) multiplexing the ages within each run.

This is a self-contained project: clone it, point `config/config.yml` at the
Cell Ranger output, and run.

## Requirements

R >= 4.3 and **Seurat v5**. The code uses the v5 `layer =` argument throughout
(`GetAssayData`, `AverageExpression`, `FindMarkers`); on Seurat v4 those calls
need `slot =` instead. `scripts/00_install_dependencies.R` installs everything
else, including the three packages that come from GitHub (monocle3,
SeuratWrappers, CytoTRACE2).

## Quick start

```bash
Rscript scripts/00_install_dependencies.R   # once per machine
ln -s /path/to/per_sample_outs data/per_sample_outs
Rscript scripts/run_all.R                   # or: run_all.R 3 4 5
quarto render analysis/report.qmd
```

Scripts are run from the project root. They also work from inside `scripts/`;
nothing else uses the working directory, and no path is hard-coded outside
`config/config.yml`.

## Layout

```
config/config.yml   every path, threshold and tuning parameter
R/                  the analysis functions, one file per concern
scripts/            the pipeline, one numbered step per stage
analysis/report.qmd the write-up; reads results, recomputes nothing
data/               inputs (not tracked)
results/            objects, tables, figures (not tracked)
docs/               provenance and notes
```

| File | What it holds |
|---|---|
| `R/setup.R` | project root, config loading, path and I/O helpers |
| `R/packages.R` | package loading with one actionable error |
| `R/io.R` | reading Cell Ranger matrices |
| `R/preprocess.R` | ADT/HTO assays, demultiplexing, QC, doublets |
| `R/dimred.R` | the RNA / SCT / ADT / WNN reduction recipes |
| `R/annotate.R` | SingleR, ImmGen label map, CytoTRACE2, module scores |
| `R/trajectory.R` | monocle3 trajectories and pseudotime correlations |
| `R/de.R` | age trends, permutation-calibrated DE, NB-GLM tests |
| `R/gsea.R` | sex x age interaction ranking and GSEA |
| `R/plots.R` | every figure |

## Pipeline

| Step | Script | Produces |
|---|---|---|
| 1 | `01_preprocess.R` | `<sex>_filtered.rds` — QC-filtered, demultiplexed, ADT/HTO assays attached |
| 2 | `02_annotate.R` | `<sex>_annotated.rds`, `<sex>_neustem.rds` — reductions, SingleR labels, CytoTRACE2 potency |
| 3 | `03_merge.R` | `bm_merged.rds`, `gmp_neutrophils.rds` — both sexes, granulocytes staged |
| 4 | `04_sex_differences.R` | per-age male-vs-female DEGs (RNA and ADT) plus GO ORA |
| 5 | `05_trajectory.R` | monocle3 trajectories per sex x age, tradeSeq GAMs |
| 6 | `06_age_trends.R` | Spearman age trends and permutation-calibrated pairwise DE, per stage |
| 7 | `07_glm_models.R` | NB-GLM omnibus, age x sex interaction, and the trajectory model |
| 8 | `08_gsea.R` | sex x age interaction GSEA, run within each maturation stage, with sensitivity and per-age checks |
| 9 | `09_development_shifts.R` | stage composition across age, and where cells sit along the trajectory |

Steps 1-3 must run in order. Steps 4-9 depend only on step 3, except that
steps 7 and 9 also use the trajectory built in step 5, so they can be run
individually as parameters are tuned.

Step 9 is the one that asks whether *development* moves rather than whether
expression does — the stage mix and the distribution of cells along
pseudotime. It is the most direct answer to the question the project was built
around, so it is worth running even when the expression steps are not.

Step 7 is by far the slowest: it refits every model `glm_de.n_perm` times to
build the permutation null. Lower that value in the config for a quick pass.

## Configuration

Everything tunable lives in `config/config.yml`: input paths, sample IDs and
their hashtag-to-age maps, per-sample PCA dimensions, QC thresholds, clustering
resolution, the cluster-to-stage mapping, and the parameters for each
statistical method. No script contains a threshold or a path of its own.

One setting changes results rather than performance:

- **`analysis.neutrophil_cluster_labels`**. A manual mapping from cluster id to
  maturation stage. Cluster numbering is not stable across Seurat versions, so
  re-check it against `results/figures/gmp_neutrophil_modules.pdf` whenever
  step 3 is re-run.

## Interactive steps

`monocle3::order_cells()` opens a Shiny window when no root is supplied, which
stalls a batch run. `learn_trajectory()` therefore takes a `root_group`
argument and picks the root programmatically as the principal-graph node
closest to the GMP cells — the standard monocle3 recipe. Pass
`root_group = NULL` in an interactive session to choose the root by hand
instead.

## Troubleshooting

**`could not find function "assay<-"`**, raised from inside `counts<-` or
another Bioconductor call. The named package is not the problem;
SummarizedExperiment is loaded but not *attached*. S4 method dispatch resolves
generics through the search path, so calling everything as `Pkg::fun()` is not
sufficient — a method running inside another package can still fail to find a
generic it needs. `require_packages()` attaches rather than merely checking,
and each step declares the packages it needs at the top of the script. If a
new step hits this, add the package to that step's `require_packages()` call.

**`FindClusters` fails on `algorithm = 4`.** That is Leiden, which needs the
Python `leidenalg` package through reticulate. Set `clustering.algorithm: 1`
in the config for Louvain, which needs nothing extra.

**A step stalls with no output.** `monocle3::order_cells()` opens a Shiny
window when it has no root. Every pipeline call passes `root_group`, so this
should not happen in a batch run; if it does, the root selection is falling
through — check that `fine_neu_labels` carries a `GMPs` level.

## Reading the results

`docs/analysis-notes.md` covers what each statistic means, which comparisons
are and are not calibrated for pseudoreplication, and the confound between sex
and batch. Read it before quoting a number.

`docs/original/bone_marrowaging.qmd` is the notebook this project was built
from, kept for provenance. `docs/changes-from-original.md` lists what moved,
what was deduplicated, and the handful of genuine bugs that were fixed.
