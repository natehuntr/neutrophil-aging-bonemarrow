# Does this pipeline answer the question?

> **How does neutrophil development change with age, and how does this differ
> between sexes?**

Short version: the code is a competent instrument pointed at the right target,
and it will answer the *age* half descriptively. The *sex* half is not a code
problem — the experimental design cannot separate sex from batch, and no
analysis in this repository or any other can fix that after the fact.

Three separate things determine the answer, and they are worth keeping apart:
what the design can support, what the pipeline actually computes, and what it
does not compute but could.

---

## 1. What the design can support

| | |
|---|---|
| Sexes | 2, **one 10x run each** |
| Ages | 4 hashed *within* each run (3m, 9m, 12m, 18m); 3 used |
| Replicates per age within a sex | **1 hashtag = 1 group, no biological replication** |

**Age comparisons are design-clean.** All four ages share a library, a
dissociation, a reverse transcription and a sequencing run. Anything technical
that varies between runs applies equally to every age, so an age difference
cannot be a batch difference. This is the strong part of the design and it is
what makes the age half of the question answerable at all.

**Sex is perfectly confounded with run.** Every male cell came from one
library and every female cell from another. A male-vs-female difference and a
library-to-library difference are the same number; there is no contrast in the
data that separates them. Reporting "gene X is higher in males" from this
design means "gene X is higher in library A".

The **age x sex interaction is the one cross-sex statistic worth reporting**,
because a constant per-library offset cancels out of a difference of slopes.
That is why the pipeline's headline cross-sex analysis (step 8, and the
interaction model in step 7) is built as an interaction rather than as a
per-age sex contrast. It is not immune: a library whose capture efficiency
compresses or stretches dynamic range produces a difference in *slopes*, not
just in means, and that does not cancel. It is a weaker confound than the main
effect, not an absent one.

**There is no biological replication at all.** One hashtag per age means each
age/sex cell population comes from one animal or one pool. Cells within it are
pseudoreplicates. A gene that differs because that particular mouse had a
subclinical infection is indistinguishable from a gene that differs because of
age. Every p-value in the pipeline — permutation-calibrated ones included —
is computed over cells, so all of them answer "are these cells different?"
rather than "are mice of these ages different?".

### What this means for a claim

| Claim | Supportable? |
|---|---|
| "Granulopoiesis in *this animal set* looks different at 18m than 3m" | Yes |
| "Neutrophil maturation shifts with age in mice" | Hypothesis, needs replicate animals |
| "Males and females differ at 18m" | No — that is a library contrast |
| "The age trajectory differs by sex" | Weakly; interaction cancels offsets but not slope artefacts |

---

## 2. What the pipeline actually computes

It is genuinely well aimed at "development", not just at "neutrophils". Cells
are staged (GMP → proNeu → proNeu2 → preNeu → immature → mature) and most
analyses are run *within* stage, which is the right shape for the question.

| Analysis | Step | What it contributes |
|---|---|---|
| Stage assignment from maturation modules | 3 | Makes "development" an axis, not a label |
| Age trends within each stage (Spearman rho, shape clusters) | 6 | Cell-intrinsic change per stage — the core age result |
| Permutation-calibrated pairwise age DE, RNA and ADT | 6 | Same question, effect-size based, less p-value dependent |
| NB-GLM omnibus `~ age` per sex, and `age x sex` interaction | 7A | Proper count model; the interaction is the defensible cross-sex test |
| Trajectory model `~ ns(pt) * age * sex` | 7B | **The single most on-question analysis**: does a gene's *developmental profile* change with age, with sex, or in shape |
| tradeSeq GAMs per sex, conditioned on age | 5 | Smoothers to visualise what 7B tests |
| Sex x age interaction GSEA | 8 | Pathway-level version of the interaction |

Step 7B is the analysis that most directly asks the question as posed. If
only one result from this repository is to be reported, it is that one.

Two methodological choices deserve credit: selection thresholds come from
refitting on shuffled labels rather than from BH-adjusted per-cell p-values,
and the sex x age GSEA refuses to label a pathway "male-skewed ageing",
instead reporting which per-gene pattern its leading edge is actually made of.
Both push against over-claiming.

---

## 3. What is missing, and could be added

These are gaps in the analysis, not in the design — the data supports all
three, and each is closer to the question than some of what is already there.

### 3.1 Stage composition is never tested — the largest gap
**Addressed** in `scripts/09_development_shifts.R`, part A.

The most direct way neutrophil development changes with age is that the
*proportions* of the stages shift: a shift toward immature cells is the
classic emergency-granulopoiesis / myeloid-bias readout of an aged marrow.

The pipeline assigns every cell a stage and then never tests whether the stage
distribution depends on age or sex. It is not a hard analysis — a
multinomial or Dirichlet-multinomial model of stage counts by age and sex, or
`scCODA`/`propeller` — and with one library per sex it carries the same
confound caveat as everything else cross-sex, while the *age* comparison
within a sex remains clean.

Without it, a strong result reported from the current pipeline could be
composition in disguise. Which leads directly to:

### 3.2 Composition and cell-intrinsic change are not separated
**Addressed**: step 8 now loops over `gsea.stages` and never pools them.

Step 8 computes each sex's age trend on **mature and immature cells pooled**.
If the immature fraction rises with age, a gene that is simply higher in
immature cells will show an age trend in the pooled data while changing in
neither stage. The GSEA ranking, and therefore the headline cross-sex result,
is open to that.

Two fixes, either sufficient: run the rho contrast within a single stage, or
carry stage as a covariate. The per-stage analyses in step 6 do not have this
problem — they are stratified — which makes the inconsistency between step 6
and step 8 worth resolving.

### 3.3 Nothing tests whether cells *move* differently along the trajectory
**Addressed** in `scripts/09_development_shifts.R`, part B.

Pseudotime is fitted and used to test gene expression, but the distribution of
cells along it is never compared between groups. "Development changes with
age" could mean cells pile up before a maturation step, or traverse a
compressed range — both visible as a shift in the pseudotime density, testable
with a two-sample Kolmogorov–Smirnov or Wasserstein distance per group, and
neither currently computed. This is arguably the most literal reading of the
question and it is a handful of lines.

Related: CytoTRACE2 potency is computed and used only as a filter. Comparing
potency distributions across ages is free and on-topic.

### 3.4 Cross-age gene-list overlap is confounded by cell number

Step 5 writes "genes unique to 3m/12m/18m" from Moran's I significance per
age, and draws a Venn. Group sizes differ several-fold, so a gene "unique to
3m" is largely a statement about which age had the most cells. Set overlaps of
per-group significance are a power comparison dressed as a biological one.
Compare effect sizes across ages, or downsample to a common n, or drop the
comparison.

### 3.5 Smaller points

- **Sex-chromosome genes** (`Xist`, `Ddx3y`, `Uty`, …) are flagged in GSEA
  leading edges but left in the ranking, where they are guaranteed top hits of
  any male-vs-female contrast and carry no information about granulopoiesis.
  Exclude them from the ranking and report them separately as a sanity check
  that the hashtag-to-sex assignment is right.
- **Stage labels are a hand-made cluster-id mapping** in the config. Every
  stage-stratified result depends on it, and cluster numbering is not stable
  across re-runs. A marker-score-based assignment would make the staging
  reproducible rather than a step someone has to re-check.
- **The ADT panel is under-used for staging.** CD101, Ly6G and CXCR2 are the
  canonical surface markers of neutrophil maturation, and the panel is used
  for the WNN graph and for age-changing ADT features, but not to define or
  validate the stages themselves — where it would be more robust than
  transcript modules.
- **9m is dropped everywhere**, leaving three timepoints. That is defensible
  (it is absent from some strata), but it makes every "trend" a statement
  about three numbers.

---

## 4. Verdict

**Age half:** the pipeline answers it as well as the data allows. The
within-run design is clean, the analyses are stage-resolved and
permutation-calibrated, and step 7B tests the developmental profile directly.
The results are descriptive of these animals; calling them a property of
ageing requires replicate mice.

**Sex half:** not answerable. One library per sex means every sex difference
is also a batch difference. The interaction analyses are the right way to make
the least-bad cross-sex statement, and they should be reported as
hypothesis-generating with the confound stated in the same sentence.

**Done since this was written:** §3.1, §3.2 and §3.3 are implemented — stage
composition and trajectory position in step 9, and step 8 stratified by stage.

**Still open, in order of value per hour:**

1. Drop or replace the per-age gene-list Venn (§3.4), which currently compares
   power rather than biology.
2. Exclude sex-chromosome genes from the per-age ranking (§3.5), where they
   are guaranteed top hits that say nothing about granulopoiesis.
3. Replace the hand-made cluster-to-stage mapping with a marker-score
   assignment (§3.5), so staging survives a re-run.
4. For a publishable sex claim, there is no analysis: it needs libraries
   containing both sexes, hashed the way the ages already are. Two more runs
   with sexes multiplexed within each would make the sex contrast as clean as
   the age contrast already is.

**What to read first, now that step 9 exists.** `stage_age_trends.csv` and
`pseudotime_shifts.csv` are the most direct answers to the question, and they
are also the two results least dependent on modelling choices: a proportion
and a distributional shift, both computed within a sex, where the design has
no batch confound. If the composition shift and the pseudotime shift point the
same way, and the potency shift agrees with them, that is a real finding about
these animals. If they disagree, the expression-level results should be
treated as unexplained until they are reconciled.
