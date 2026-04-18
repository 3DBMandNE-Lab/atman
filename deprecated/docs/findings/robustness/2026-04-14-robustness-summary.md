# Week 2 Robustness Summary (Executed 2026-04-14)

## Scope

Robustness profiling was executed without adding new public CLI commands.

Runs:

1. Leave-one-subject-out (LOO) DE reruns for all 10 subjects.
2. Threshold sweep for `min_pairs = 4, 5, 6`.
3. Top-rank overlap vs baseline for acute contrasts.

Artifacts:

1. `loo_gene_stability.tsv`
2. `threshold_sensitivity.tsv`
3. `rank_stability.tsv`

## Key Results

### 1) Threshold sensitivity (`min_pairs`)

`threshold_sensitivity.tsv` is identical across `min_pairs = 4, 5, 6` on this dataset:

1. `PT1-PR1`: `q<0.05 = 36`, `q<0.10 = 62`
2. `PT2-PR2`: `q<0.05 = 37`, `q<0.10 = 277`
3. `PT2-PT1`: `q<0.05 = 7`, `q<0.10 = 8`
4. `PR2-PR1`: `q<0.05 = 2`, `q<0.10 = 7`

Interpretation: headline counts are not sensitive to this `min_pairs` range.

### 2) LOO sign stability

From `loo_gene_stability.tsv`:

1. `PT1-PR1`: sign match rate `0.716` (across all tested rows), `834` rows with at least one sign flip.
2. `PT2-PR2`: sign match rate `0.731`, `791` rows with at least one sign flip.

Interpretation: many borderline proteins are LOO-volatile, which is expected at `n=10`.

### 3) LOO retention for baseline acute `q<0.05` proteins

1. `PT1-PR1`: baseline `36` proteins at `q<0.05`
   - retained at `q<0.05` in all 10 LOO runs: `11`
   - retained at `q<0.05` in at least 8/10 LOO runs: `11`
   - retained at `q<0.10` in at least 8/10 LOO runs: `30`
2. `PT2-PR2`: baseline `37` proteins at `q<0.05`
   - retained at `q<0.05` in all 10 LOO runs: `5`
   - retained at `q<0.05` in at least 8/10 LOO runs: `7`
   - retained at `q<0.10` in at least 8/10 LOO runs: `33`

Interpretation: strict `q<0.05` calls are fragile; directional and `q<0.10` signals are more stable.

### 4) Top-rank overlap stability (top-20 proteins)

From `rank_stability.tsv`:

1. `PT1-PR1`: mean overlap `16.0 / 20`, mean Jaccard `0.673`.
2. `PT2-PR2`: mean overlap `14.5 / 20`, mean Jaccard `0.589`.

Interpretation: top acute ranking is moderately stable under subject perturbation.

### 5) Canonical heat-stress proteins under LOO

For acute comparisons (`PT1-PR1`, `PT2-PR2`):

1. `HSPB1`: sign-match `10/10` in both acute comparisons.
2. `HSPA1A`: sign-match `10/10` in both acute comparisons.
3. `DNAJB1`: sign-match `9/10` (`PT1-PR1`), `10/10` (`PT2-PR2`).
4. `HSPG2`: sign-match `10/10` in both acute comparisons.

Interpretation: canonical acute heat-stress directionality is robust.

## Claim Guidance

Confirmed:

1. Acute heat response directionality (including canonical heat-shock proteins) is robust.
2. Global headline discovery counts are stable to reasonable `min_pairs` choices.

Exploratory:

1. Individual non-canonical `q<0.05` protein calls are often LOO-fragile at `n=10`.
2. Fine-grained protein ranking near thresholds should be presented with stability context.

