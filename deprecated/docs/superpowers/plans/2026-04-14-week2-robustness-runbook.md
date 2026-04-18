# Week 2 Robustness Runbook (No New CLI Surface)

## Goal

Quantify how stable key acute findings are under subject perturbation and threshold changes.

## Inputs

1. `/tmp/karna_showcase/qc_measurements.tsv`
2. `/tmp/karna_showcase/samples.tsv`
3. `/tmp/karna_showcase/de_results.tsv`

## Outputs (target)

1. `docs/findings/robustness/loo_gene_stability.tsv`
2. `docs/findings/robustness/threshold_sensitivity.tsv`
3. `docs/findings/robustness/rank_stability.tsv`
4. `docs/findings/robustness/2026-04-14-robustness-summary.md`

## Execution Steps

1. Re-run DE ten times, each time excluding one subject from all four conditions.
2. For each rerun, compute:
   - sign consistency for acute contrasts (`PT1-PR1`, `PT2-PR2`),
   - retention of `q<0.05` and `q<0.10` hits.
3. Sweep thresholds:
   - `min_pairs`: `4, 5, 6`,
   - `q`: `0.05, 0.10`.
4. Rank-stability:
   - compare top-20 gene overlap across runs (Jaccard / overlap counts).
5. Write one summary markdown with:
   - stable proteins (high retention),
   - unstable proteins (high volatility),
   - impact on biological headline claims.

## Acceptance Criteria

1. Acute heat headline remains directionally stable in the majority of LOO runs.
2. Volatile findings are explicitly flagged as exploratory.
3. No changes to public command interfaces.

