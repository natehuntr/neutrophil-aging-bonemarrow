# Methods

Single-cell analysis of bone marrow granulopoiesis across age and sex.
Written to describe the pipeline in `scripts/01`–`scripts/10`; every threshold
quoted here is read from `config/config.yml` at run time, and the values given
are those used for the reported analysis.

---

## Experimental design and its constraints

Two CITE-seq libraries were prepared, one per sex. Within each library, four
ages (3, 9, 12 and 18 months) were multiplexed by hashtag antibody, with six
mice pooled per hashtag (48 animals in total).

This design has three consequences that determine what the analysis can and
cannot claim, and they are stated here because several methodological choices
below exist only to accommodate them.

**Sex is completely confounded with library.** It is equally confounded with
sequencing run, probe barcode and staining batch. No analysis can separate a
sex difference from a batch difference, and no statistical treatment recovers
the distinction. Comparisons *between* sexes are therefore reported as ranked
candidates for orthogonal validation rather than as tests (see *Sex-difference
candidates*). Comparisons *across ages within a sex* share a library and are
design-clean; these carry the primary findings.

**There is no biological replication at the analysis level.** One hashtag per
age per sex means cells are pseudoreplicates of a pooled group. Pooling six
animals improves the *precision* of each group estimate — roughly σ/√6 for a
group mean — but does nothing for *accuracy*: the depth asymmetry, the
sex–library confound and any differential cell recovery are unaffected by how
many animals were pooled. Because cells cannot be traced to animals,
between-animal variance is unmeasurable, so precision cannot be quantified
either. No inferential statistic in this pipeline describes mice.

**Individual-mouse resolution is unrecoverable.** The animals are inbred, so
genotype demultiplexing has no SNPs to work with, and mitochondrial
demultiplexing fails for the same reason compounded by the poor mitochondrial
coverage of Flex chemistry. Pseudobulk aggregation by animal is impossible
rather than merely inconvenient.

The two libraries also differ in sequencing depth by a factor of 1.82
(median UMI per cell, female deeper). This is addressed explicitly in *Depth
matching*.

---

## Software

Analysis was performed in R 4.3.2 with Bioconductor 3.18. Core packages:
Seurat v5 and SeuratObject 5.4.0; SoupX 1.6.2; scDblFinder; scuttle;
DropletUtils; SingleR with celldex 1.12.0 (ImmGen reference); CytoTRACE2;
UCell 2.6.2; monocle3 1.4.27; tradeSeq; glmGamPoi; fgsea 1.28.0 with msigdbr.
A fixed seed (42) was set at the start of every step, and all randomised
procedures (permutations, downsampling, random splits) derive their seeds from
it.

---

## Preprocessing and quality control

**Ambient RNA.** Contamination was estimated per library with SoupX
`autoEstCont()` on quick clusters (30 PCs, Louvain resolution 0.5) computed
from the raw and filtered Cell Ranger matrices, and removed with
`adjustCounts(roundToInt = TRUE)`. Estimated contamination was 1.4% (male) and
1.0% (female).

**Hashtag demultiplexing.** HTO counts were CLR-normalised across cells
(`margin = 2`) and demultiplexed with `MULTIseqDemux(autoThresh = TRUE)`. Cells
called Negative or Doublet were removed, and the remaining hashtag calls were
translated into age labels.

**Doublets.** Called with scDblFinder on the RNA assay; 5.3% of cells in both
libraries.

**Cell filtering.** Demultiplexing and doublet calls were applied first, then
quality thresholds, so that retention can be attributed to the correct cause.
Count- and feature-based outliers were identified with
`scuttle::perCellQCFilters()` at 3 MADs, computed **per library**
(`batch = sex`) rather than globally. This matters: the libraries differ
~1.8-fold in depth, and a single global threshold applied to both would
over-filter the shallower library and manufacture exactly the composition
difference the study is trying to measure. The per-library thresholds selected
were 93 UMI / 65 genes (male) and 161 UMI / 128 genes (female). Cells were
additionally required to have mitochondrial fraction < 0.2 and haemoglobin
fraction < 10%.

Retention was 91.0% (male, 12231 → 11129) and 88.3% (female, 11913 → 10519), a
1.30-fold difference in discard rate. This asymmetry is reported alongside the
composition results, since differential recovery is the leading artefactual
explanation for a composition shift; the direction of the selected thresholds
(stricter in the *deeper* library) is opposite to what would be needed to
produce the observed result.

---

## Normalisation, integration and clustering

Each library was normalised twice: log-normalisation with 2000 variable
features for general use, and SCTransform (`vst.flavor = "v2"`) for dimension
reduction. ADT counts were CLR-normalised. PCA and UMAP were computed per assay
using per-sample dimensionalities set in the config (RNA 15 PCs both libraries;
SCT 30 male / 20 female; ADT and WNN dimensionalities likewise per-sample).
Weighted-nearest-neighbour integration of RNA and ADT was performed with
`FindMultiModalNeighbors()`. Clustering used the Leiden algorithm
(`algorithm = 4`); resolution was chosen at 0.4 after inspecting a
clustree across 0.2–1.6.

Libraries were merged for joint analysis. Cell identity was assigned with
SingleR against the ImmGen reference at both main and fine resolution, with
low-confidence calls pruned to NA and retained as an explicit "NA" category
rather than silently dropped. GMPs and neutrophils were subset for all
downstream analysis (4512 cells: 1921 male, 2591 female).

---

## Maturation-stage assignment

Cells were assigned to five granulopoiesis stages — GMP, proNeu, preNeu,
immature and mature — by rank-based scoring of RNA signatures with UCell
(`AddModuleScore_UCell`), each cell taking the highest-scoring signature.
UCell ranks genes within a cell before scoring, which makes it substantially
less sensitive to sequencing depth than a mean-expression module score.

An ADT panel-based assignment was implemented and run in parallel as a
comparison, but was **not** used for the primary stratification. The antibody
panel lacks CD117, CD101 and CD177, which are the markers that separate the
early stages; once unresolvable markers are dropped, the proNeu and proNeu2
panels reduce to the same two antibodies and the GMP panel retains no positive
marker at all. `score_adt_panels()` reports each panel as it actually resolves
and warns when two collapse together. A separate proNeu2 stage was merged into
proNeu for the same reason: neither modality could separate them (3 of 4512
cells were assigned to proNeu2 by RNA), and a stage that exists in the
configuration but not in the data fails every subsequent gate.

The two assignments were crossed and agreed on 63.6% of the 3885 cells labelled
by both. Because a confusion matrix alone cannot say whether a disagreement
matters, assignment margins (the gap between the best and second-best score)
were compared between agreeing and disagreeing cells, per stage. Disagreements
concentrated at low margins indicate a fuzzy boundary neither method can
adjudicate; disagreements at margins comparable to the agreements indicate a
genuine conflict. This is reported per stage rather than pooled, since the
stages behave in opposite directions.

---

## Depth matching

Because the libraries differ 1.82-fold in median depth, and because sequencing
depth influences gene detection, transcriptional complexity and therefore both
pseudotime position and potency scores, a depth-matched count assay
(`RNAmatched`) was constructed once and used by every analysis that compares
across strata.

Counts were thinned by binomial downsampling
(`DropletUtils::downsampleMatrix`) toward the lowest group median, with cells
already at or below the target left unmodified. Matching was performed on the
`age × sex` grouping — the finest grouping any downstream contrast uses — so
that the sex contrast and the within-sex age contrasts read the same assay.
The target was 2225 UMI; 3295 of 4512 cells were thinned, and the post-matching
depth ratio between sexes was 1.00.

CytoTRACE2 potency was re-scored on the matched counts, since it reads
transcriptional complexity directly and is the measure in the pipeline most
exposed to depth. Trajectories (below) were rebuilt from the matched assay with
PCA and UMAP recomputed from those counts rather than inherited from the
log-normalised embedding — the principal graph is learned on the UMAP, so
substituting the count matrix alone would leave pseudotime unmatched.

A per-stratum depth ratio is recomputed inside each stage before scoring.
Global matching thins cells *down* and leaves cells below the target alone, so
a stratum in which one group sits far below the global target is not matched at
all; candidates arising from such strata are flagged and demoted.

---

## Composition analysis

Stage proportions were tabulated by age within each sex and tested by
chi-square, with centred log-ratio transformation applied to remove the
sum-to-one constraint so that one stage rising is not confused with the others
being pushed down. Cross-sex composition comparisons are reported with an
explicit library-confound annotation attached to every row.

---

## Trajectory and developmental position

Trajectories were built with monocle3 on the depth-matched assay, one per
age × sex group and one combined across all cells, rooted programmatically at
the GMP stage (the principal-graph node closest to the majority of GMPs) so
that no interactive root selection was required. Genes associated with each
trajectory were identified by Moran's I. Per-sex generalised additive models
conditioned on age were fitted with tradeSeq `fitGAM`, restricted to genes
graph-associated in any per-age trajectory.

Distributions of pseudotime and of CytoTRACE2 potency were compared across ages
*within* each sex, as a shift relative to the youngest age, summarised by median
shift, Wasserstein distance and a Kolmogorov–Smirnov test with BH adjustment.

Both comparisons were additionally run **within each maturation stage**, gated
on a minimum cell count at every age rather than on the stage total. This
distinction is essential to interpretation: pooled across stages, both
pseudotime and potency move whenever the stage mix moves, and the stage mix does
move with age. A pooled shift is therefore consistent with no individual cell
changing. Only the within-stage comparison speaks to the cells themselves, and
the two are reported together.

A depth diagnostic accompanies these comparisons, reporting median UMI and
median genes detected per age within each sex, both pooled across stages and
per stage. Only the per-stage rows are treated as a verdict: transcriptional
complexity is a property of cell type, so a shift in the stage mix moves the
pooled figures on its own.

---

## Age-associated expression

### Monotonic trends across four ages

Within each stage × sex stratum, per-gene Spearman correlation with age was
computed and calibrated against a null generated by permuting age labels within
the stratum (100 permutations); the threshold was the 95th percentile of the
maximum |rho| across genes. Reporting excess over a stratum-specific null rather
than a fixed correlation cutoff avoids a cutoff whose yield mostly reflects
stratum size.

Permuting age destroys any association between age and library complexity, so
this null cannot calibrate against a complexity trend that is genuinely
correlated with age. Each stratum therefore also reports the Spearman
correlation between age and genes detected, compared against that stratum's own
threshold; a large gene count beside a strong detection trend is one
observation reported twice, not independent findings.

Strata were required to hold at least 50 cells at every age. Under this
requirement only 3 of 10 stage × sex strata were analysable.

### Endpoint contrast with shape check

Because the four-age requirement retires most strata — male 12m holds 24 mature
cells — a second analysis separates discovery from trend confirmation. Genes
were discovered by contrasting the extremes, 3m versus 18m, where the
difference is largest and the cells most numerous; 8 of 10 strata clear the
cell gate on those two ages.

The per-gene statistic is the difference in mean normalised expression between
the two ages. Two thresholds are computed from 100 label permutations: a
family-wise threshold (95th percentile of the maximum absolute difference across
genes) and a false-discovery threshold, obtained by finding the smallest cutoff
at which the expected number of permuted genes exceeding it, divided by the
observed number, falls at or below the target (0.10). Selection uses the
false-discovery threshold, with the family-wise value reported alongside and
genes clearing it flagged separately. Where no cutoff reaches the target, the
family-wise threshold is used and the output records that this occurred.

For the discovered genes, the two intermediate ages were then read off as a
shape check: each intermediate mean is expressed as a fraction of the 3m → 18m
difference, and genes classified as monotonic, between endpoints but out of
order, reversing before the endpoint, or overshooting it. **This is a ranking
signal, not validation** — it reuses the cells that found the genes, and several
intermediate strata are small. Per-age cell counts accompany every table.

---

## Gene set enrichment

Preranked GSEA (fgsea `fgseaMultilevel`, `eps = 0`) was run over GO:BP gene
sets from MSigDB (collection C5, subcollection GO:BP), restricted to sets of
15–500 genes present in the stratum's expressed universe. Two rankings were
used: the sex × age interaction statistic (within the mature and immature
stages, three ages), and the 3m → 18m endpoint effect (all analysable strata).
Redundant sets were collapsed with `collapsePathways()` and the surviving
independent sets flagged, since GO:BP is heavily nested and one signal otherwise
reports as many findings.

---

## Sex-difference candidates

This pipeline does not test sex differences. With sex confounded against
library, run, probe barcode and staining batch, no p-value from these data
describes mice. Step 4 instead produces a ranked, robustness-annotated candidate
list for orthogonal validation, and inferential columns are removed from its
output rather than caveated.

Gene sets were scored per cell with UCell on the depth-matched assay, within
each maturation stage. Only sets for which an orthogonal bench assay could be
proposed were scored — a candidate that cannot be tested is not a candidate.
For each set, the standardised mean difference between sexes was computed with a
bootstrap confidence interval (1000 resamples).

Each candidate carries four robustness annotations:

- **Direction relative to the confound.** The deeper library is female, so a
  female-high result points in the direction the confound already predicts.
  Candidates running *against* the bias are held to a lower evidential bar than
  those running with it, and this is recorded per row.
- **Within-library empirical null.** The same module was scored across 200
  random splits of the deeper library alone, matched on stage and cell number.
  Anything the module finds there is the floor for what it can find between
  sexes. Candidates below the median of this null are demoted outright: a
  random split of one library separating the sexes better than the sexes do is
  evidence against the candidate, not weak evidence for it.
- **Survival under depth matching**, including whether the effect changes sign.
- **Stratum depth ratio**, recomputed after matching within that stratum.

Three depth controls accompany the primary comparison: whether the effect
survives matching the deeper library down; whether downsampling alone
manufactures the effect within a single library, where no sex difference can
exist; and whether the effect is a cell-number artefact rather than a depth one.

---

## Control gates

Analyses are preceded by gates that block output when a known failure mode is
present. Each corresponds to a way these data have been observed to mislead.

- **Sex-chromosome positive control.** Xist, Ddx3y, Uty, Kdm5d and Eif2s3y must
  rank near the top of a male-versus-female contrast. The gate distinguishes
  "absent because the probe set does not target it" from "absent despite being
  targeted", which look identical in a results table and have opposite
  implications. All five are present in this probe set.
- **Depth ratio**, evaluated on the matched counts the analysis reads, not on
  raw counts (limit 1.3-fold).
- **Minimum cells per stratum** (50). Scoped to the strata it names: undersized
  strata are excluded and the analysis continues, and the gate fails only when
  no stratum is usable.
- **ADT isotype controls** must come out non-significant, as a negative control
  on the protein layer.

Gates are fatal by default (`gates.fatal: true`), because the failure mode they
guard against produces results that look correct.

---

## Reproducibility

All parameters are held in a single `config/config.yml`; no threshold is
hard-coded in an analysis script. Each step reads named objects from disk and
writes named objects and tables, so any step can be re-run independently.
Package versions are pinned to a CRAN snapshot and Bioconductor 3.18 via
`scripts/00_install_dependencies.R`. Jobs are submitted as one SLURM job per
step with explicit dependencies (`slurm/submit_all.sh`).
