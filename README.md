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
| - | `diagnose.R` | inspects the saved objects and reports what any later step would fail on |

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

**`The total size of the N globals exported for future expression is ... MiB.
This exceeds the maximum allowed size 500.00 MiB`**, usually from
`run_sct_reduction`. Seurat routes SCTransform through the `future` package,
whose default ceiling on data captured by a closure is 500 MiB; SCTransform's
`conserve.memory` path captures the count matrix and exceeds it on a full
sample. Raise `compute.globals_max_size_gb` in the config (default 8). It is a
transfer ceiling, not an allocation — raising it reserves no memory.

If the machine is short on RAM, note that a parallel `compute.future_plan`
makes this *worse*, not better: every worker gets its own copy of those same
large globals. `sequential`, the default, is usually both faster and lighter
here.

**`'X' is not an exported object from 'namespace:Seurat'`**, e.g. for
`JoinLayers`. Those functions are defined in **SeuratObject**, and which of
them Seurat re-exports varies between v5 releases. `seurat_fn()` in
`R/dimred.R` looks in both namespaces; use it (or `join_layers()` /
`assay_names()`, which are built on it) rather than hard-coding
`Seurat::thing` for anything that originates in SeuratObject.

**`No cells found` from `subset()`.** Run `Rscript scripts/diagnose.R`: it
reports every saved object's cell count and, for each metadata column the
later steps select on, how many cells carry a value. A column marked absent or
`non_na = 0` is the cause, and the step that should have written it is the one
to re-run.

The usual culprit is `CytoTRACE2_Potency`, which steps 4 and 7 use to pick
differentiated cells. It is written in step 2 and copied onto the merged
object in step 3; both now verify it rather than passing the problem on.
Filters report a per-condition breakdown when they empty out, so the error
names the condition that did it.

**`Tibble columns must have compatible sizes`** from inside `tradeSeq::fitGAM`.
tradeSeq drops cells it cannot model — no finite pseudotime, or zero counts
across the genes supplied — but passes `conditions` through untouched, so a
conditions vector built from the unfiltered object arrives longer than the
model matrix and the fit dies in a `tibble()` call several frames down. Step 5
now filters counts, pseudotime, weights and conditions to one common set of
cells and asserts they agree before calling `fitGAM`. Cells at a timepoint
outside `analysis.age_levels` are the usual source: they become NA conditions.

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
