# Supplementary — simulation study for the S-score

## Design

Synthetic paired proteomics matrices with known ground truth. Each cell in
the grid below generates a per-rep random assignment of true positives
and effect sizes, draws Gaussian per-subject noise, runs the full
per-protein pipeline (paired-t + BH, LOO, S-score) on the synthetic data,
and records top-K precision / recall / F1 where K is the true-positive
count. 50 replicates per cell.

| factor | levels |
|---|---|
| n (subjects) | 8, 12, 16, 20, 24 |
| frac TP      | 0.05, 0.15, 0.30 |
| effect dist  | strong (N(1.5, 0.3)), mixed (N(0.5, 0.5)), weak (N(0.3, 0.1)) |
| variance heterogeneity | homoscedastic σ=0.7 vs per-subject σ ∼ U(0.4, 1.2) |
| fragile-TP fraction   | 0.0 (all robust), 0.25, 0.50 (outlier-driven TPs) |
| P (proteins) | 500 |

"Fragile TPs" are true positives whose effect is driven by a single-subject
spike rather than a cohort-wide signal. They look large on $|\bar d|$ but
their sign match rate under LOO drops because removing the driver subject
zeros the effect.

**Total cells: 5 × 3 × 3 × 2 × 3 = 270.** Implementation:
`scripts/simulate_s_score.py`. Full grid runs in ~10 s.

## Headline result

**$|\bar d|$-ranking dominates discovery F1 in every regime tested.**

| ranking | cells won by F1 (of 270) |
|---|---|
| `\|d\|`     | **258** |
| $S$          |  11 |
| $q$          |   1 |

$S$ never beats $|\bar d|$ by ≥0.02 F1 in any cell. $S$ beats $q$ by ≥0.02
F1 in 229/270 cells. $|\bar d|$ beats $q$ in every cell.

## Where the gap between |d̄| and S widens

| fragile_frac | q mean F1 | \|d\| mean F1 | S mean F1 | S − \|d\| |
|---|---|---|---|---|
| 0.0  | 0.675 | **0.695** | 0.691 | −0.004 |
| 0.25 | 0.538 | **0.644** | 0.577 | −0.067 |
| 0.50 | 0.402 | **0.606** | 0.484 | −0.121 |

At fragile_frac = 0, the three rules are close and the 11 cells where $S$
wins are all in this regime. As fragile-TP proportion grows, $|\bar d|$
pulls ahead substantially because fragile TPs are genuine true positives
in the ground truth and $S$'s SSR factor demotes them. $q$ suffers most
from fragility because single-subject spikes inflate the paired-t SD
denominator and depress all q-values.

## What this means for the S-score

The S-score is **not an F1-optimal discovery ranking**. It is a
stability-prioritised alternative with specific semantics:

- $S$ rewards **cohort-general signals** (consistent across most subjects).
- $|\bar d|$ rewards **signals of any kind**, including single-subject drivers.
- $q$ imposes a correctness filter (BH-FDR) but is fragility-sensitive at small n.

Which is correct depends on the analyst's question. A researcher looking for
**population-level biomarkers** that generalise to new cohorts is served by
$S$ (or by any ranking that penalises single-subject effects): a protein
whose signal depends on one subject will not survive replication. A
researcher looking for **heterogeneity-aware discovery** — including
subgroup-driven or individual-specific response proteins — is served by
$|\bar d|$: these are legitimate real effects in the biology and $S$
throws them away.

The $S$-score's value is therefore **not a bigger true-positive yield** but
a specific inductive prior: give me the proteins whose effect is stable
across subjects, at the cost of some true positives that happen to be
subject-specific. On the real Dube and Gisby datasets in the main text,
every q<0.20 hit has SSR=1.00 and the SSR factor is inactive — so in
practice $S$ and $|\bar d|$ produce near-identical orderings (both
dominated by effect magnitude). The simulation regime where they diverge
(fragile_frac > 0) does not occur in the main-text cohorts; if it did, the
data-consistent rule would be $|\bar d|$, not $S$.

## What the simulation does and does not say

**Says:**
1. $q$-ranking is rarely the best discovery ranking at small n.
2. $|\bar d|$-ranking is empirically the most F1-optimal across a wide grid
   of n, effect size, variance heterogeneity, and TP fragility.
3. $S$-ranking is systematically between $q$ and $|\bar d|$.
4. The gap between $S$ and $|\bar d|$ is driven by TP fragility — when
   all TPs are robust cohort-level signals, the gap is negligible
   (mean F1 difference −0.004).

**Does not say:**
1. That $S$ is useless. LOO rank-stability (Fig 5A in main text)
   remains a property of $S$ that neither $q$ nor $|\bar d|$ has as a
   single column.
2. That a dataset's optimal ranking is always $|\bar d|$. The
   simulation uses Gaussian noise; real proteomics data may have
   structure this does not capture (batch effects, technical outliers,
   non-normal within-subject distributions).
3. That the S-score should be removed from the pipeline. Its role is
   prioritisation with a stability prior, not recovery improvement.

## Recommendation adopted in the main text

- Use $|\bar d|$-ranking (or the q-ordered DE table) for discovery.
- Use $S$-ranking when the analytical goal is explicitly cohort-general
  effects and you want to penalise single-subject drivers.
- Do not present $S$ as a recovery improvement over $|\bar d|$.

The main-text claim for $S$ is reduced to its defensible empirical
property: lower top-K Jaccard variance under LOO resampling than
$q$-ranking on the main-text cohorts, and a guardrail against
direction-unstable proteins dominating top-ranked sets.

## Reproduce

```bash
python3 scripts/simulate_s_score.py --P 500 --reps 50 \
    --output out/simulations/s_score_regime_map.tsv
```

~10 seconds.
