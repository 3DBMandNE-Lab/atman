# Supplementary — null calibration for the module-DE pipeline

## Motivation

The module-DE claim in the main text — that data-driven module-level testing
recovers biology that per-protein DE disperses into many individually weak
hits — is only publishable if the FDR control at the module level is valid
under the null. A module-DE pipeline that re-learns modules on the same
permuted data used for the test could in principle exhibit "clustering on
noise" inflation: modules constructed to look maximally coherent on a random
realisation of the data will have non-random module-score differences and
apparently significant contrasts. This supplementary tests that
directly.

## Permutation procedure

For each of $B = 1{,}000$ sign-flipping permutations, for each subject
independently with probability 0.5, the paired condition labels are swapped
(A ↔ B). This is the standard paired-design null: it breaks the systematic
A − B effect while preserving within-subject correlation structure.

Under each permutation the **full pipeline is re-run** (not held fixed):

1. Sign-flipped (A, B) measurements are recomputed per subject.
2. Per-subject gene deltas are recomputed under the permuted labels.
3. Module-learning runs on the permuted deltas (complete linkage on
   $(1 - r_{ij})/2$ distance, cut at the same K as observed).
4. Module-level paired-t with BH-FDR across K modules.
5. Q-values of all K modules are recorded.

Implementation: `scripts/calibrate_module_de.py`.

## Results on three contrasts

| dataset | contrast | K | B | $\alpha$ = 0.01 | 0.05 | 0.10 | expected ($\alpha \cdot K$) |
|---|---|---|---|---|---|---|---|
| Dube   | PT1-PR1    | 20 | 1000 | 0.000 | 0.000 | 0.000 | 0.2 / 1.0 / 2.0 |
| Dube   | PT2-PR2    | 20 | 1000 | 0.000 | 0.012 | 0.030 | 0.2 / 1.0 / 2.0 |
| Gisby  | late-early | 15 | 1000 | 0.005 | 0.026 | 0.053 | 0.15 / 0.75 / 1.5 |

Columns under $\alpha$ are the per-permutation mean number of modules
rejected at nominal BH q < $\alpha$. The empirical rejection rate is
**1–2 orders of magnitude below nominal** at every $\alpha$ on every
contrast — BH-FDR is conservative, not anti-conservative, on the null.

### Why so conservative?

The sign-flip null zeros the per-module mean (in expectation) but preserves
the magnitude of within-subject deltas in the denominator of the t-statistic.
Concretely, a module with true signal $\bar d = 2.0$ and sd $= 0.5$ under
the observed data has $t = 12.6$; after sign-flipping each subject with
probability 0.5 the mean drifts to approximately zero while the sd grows
to approximately $|\bar d|$ units ($\sim 2.0$), driving $t$ to near zero
and $p$-values to concentrate near 1. The empirical $p$-value distribution
under the null has median 0.65 (Dube PT2-PR2; computed from
`out/robustness/module_de_null_PT2-PR2.tsv` column `p_value`, n=20,000),
not 0.5, confirming this conservative skew.

The BH procedure is therefore **safe** in this setting: the module-DE
pipeline cannot generate false positives faster than a uniform-null BH
would predict, and on these datasets it generates them much slower. Users
who want an empirically-calibrated q can derive one directly from the
permutation distribution of module p-values; the `--calibrated-q` mode
is straightforward to add if a dataset is shown to be anti-conservative
under null, but that situation has not arisen on the cohorts tested.

## Sensitivity to linkage and K (Dube PT2-PR2, B=100 per cell)

| linkage  | K=5 | K=10 | K=15 | K=20 | K=25 | K=30 | K=40 | K=50 | silhouette (K=20) |
|----------|-----|------|------|------|------|------|------|------|---|
| single   |  1  |  1   |  2   |  2   |  1   |  1   |  2   |  1   | **−0.29** (bad) |
| average  |  0  |  0   |  0   |  0   |  0   |  0   |  0   |  0   |  0.22 |
| complete |  0  |  0   |  0   |  0   |  0   |  0   |  0   |  0   |  0.14 |
| ward     |  0  |  0   |  0   |  0   |  0   |  0   |  0   |  0   |  0.17 |

Integer cells are the observed-data hit count at q<0.05; the null mean
rejection rate is ≤0.12 at every (linkage, K) cell (detailed numbers in
`out/robustness/module_de_sensitivity_PT2-PR2.tsv`).

Observations:

1. **Single linkage is disqualified.** Silhouette scores are strongly
   negative at all K ≥ 10, indicating that single linkage produces
   clusters that are worse than random on these data. Its observed hit
   counts are larger than other linkages, but those "hits" are not
   tracking cohesive regulons — they are an artefact of chaining.

2. **Complete, average, and Ward give equivalent null behaviour** (0 hits
   at q<0.05 under null at every K). Silhouette is comparable across
   these three, with Ward and average slightly higher. Complete was
   retained as the default because it is stable under small-n
   correlation noise: average linkage exhibits severe chaining on
   high-P proteomics correlation matrices (documented in the module
   learning script's rationale), and Ward's squared-Euclidean
   semantics have awkward interpretation on $(1 - r)/2$ distance.

3. **K is elastic in this range.** For complete linkage the observed-data
   hit count is 0 across K=5–50 on this single-contrast module learning
   (contrast with the canonical main-text pipeline using multi-contrast
   features which yields K=20 → 2 hits). The sensitivity analysis uses
   contrast-specific features for clean null interpretation; the main
   text's multi-contrast learning yields slightly different module
   boundaries and correspondingly different observed hit counts. Both
   approaches exhibit well-controlled FDR under the sign-flip null.

## What the calibration does and does not prove

It proves: **the module-DE hits in the main manuscript (e.g., Gisby
`module_14` at q = 3.5 × 10⁻⁴, Dube `module_06` at q = 0.039–0.049) are
not products of FDR inflation under the sign-flip null.** The pipeline's
rejection rate under null is far below nominal at every α tested.

It does not prove: module-level power, biological novelty of the learned
modules, or optimality of the BH-on-modules approach. Those are separate
claims handled in the main Results.

## Reproduce

```bash
# Dube PT1-PR1 and PT2-PR2
python3 scripts/calibrate_module_de.py \
    --measurements out/qc_measurements.tsv --samples out/samples.tsv \
    --contrast PT1-PR1 --k-modules 20 --n-perm 1000 \
    --output out/robustness/module_de_null_PT1-PR1.tsv \
    --summary out/robustness/module_de_null_PT1-PR1_summary.tsv

python3 scripts/calibrate_module_de.py \
    --measurements out/qc_measurements.tsv --samples out/samples.tsv \
    --contrast PT2-PR2 --k-modules 20 --n-perm 1000 \
    --output out/robustness/module_de_null_PT2-PR2.tsv \
    --summary out/robustness/module_de_null_PT2-PR2_summary.tsv

# Gisby late-early
python3 scripts/calibrate_module_de.py \
    --measurements out_gisby/qc_measurements.tsv --samples out_gisby/samples.tsv \
    --contrast late-early --k-modules 15 --n-perm 1000 \
    --output out_gisby/robustness/module_de_null_late-early.tsv \
    --summary out_gisby/robustness/module_de_null_late-early_summary.tsv

# Sensitivity sweep (4 linkages × 8 K values × 100 perms)
python3 scripts/sensitivity_module_de.py \
    --measurements out/qc_measurements.tsv --samples out/samples.tsv \
    --contrast PT2-PR2 --n-perm 100 \
    --output out/robustness/module_de_sensitivity_PT2-PR2.tsv
```

Each permutation runs in ~35 ms; full 1000-perm × 3-contrast calibration
takes ~100 seconds; sensitivity sweep ~2 minutes.
