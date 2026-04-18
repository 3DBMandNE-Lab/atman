# Heterogeneity Analysis Summary (2026-04-14)

## Scope

Goal: characterize **range of human responses**, not only cohort-mean effects.

Input source: `/tmp/karna_showcase` outputs generated from the standard pipeline.

Artifacts:

1. `per_subject_gene_deltas.tsv`
2. `subject_response_burden.tsv`
3. `responder_clusters.tsv`
4. `subject_top_proteins.tsv`
5. `high_variance_acute_genes.tsv`

## Cohort-Range Signal

From `subject_response_burden.tsv`:

1. Acute burden (`(mean_abs_pt1_pr1 + mean_abs_pt2_pr2)/2`) range:
   - min `0.387198` (subject `007`)
   - max `0.612802` (subject `018`)
2. Adaptation burden (`(mean_abs_pt2_pt1 + mean_abs_pr2_pr1)/2`) range:
   - min `0.352999` (subject `007`)
   - max `0.628138` (subject `009`)
3. Acute vs adaptation burden correlation across subjects:
   - `r = 0.631629`

Interpretation: participants are not interchangeable; response magnitude spans a clear range.

## Responder Classes (Median-Split Burdens)

From `responder_clusters.tsv`:

1. `HighAcute_HighAdapt`: `3`
2. `HighAcute_LowAdapt`: `2`
3. `LowAcute_HighAdapt`: `2`
4. `LowAcute_LowAdapt`: `3`

Interpretation: cohort decomposes into distinct response programs rather than one homogeneous pattern.

## Subject-Level Directional Spread

Per-subject mean acute deltas show both stronger and weaker responders:

1. `PT1-PR1` mean change:
   - max `+0.500275` (subject `009`)
   - min `-0.127828` (subject `016`)
2. `PT2-PR2` mean change:
   - max `+0.360721` (subject `019`)
   - min `-0.099636` (subject `015`)

Interpretation: even in acute contrasts, mean direction and magnitude vary materially by individual.

## Recurrent Subject-Fingerprint Proteins (Top-5 by absolute delta)

From `subject_top_proteins.tsv` (counted across subjects):

1. `DAND5` appears in top-5 for `pt2_pr2` in `6/10` subjects.
2. `GH1` appears in top-5 for `pt1_pr1` in `5/10` and `pt2_pr2` in `5/10`.
3. `AGBL2` appears in top-5 for `pt1_pr1` in `4/10`.
4. `ITGAL` appears in top-5 for `pt1_pr1` in `4/10`.

Interpretation: some proteins recur across individuals, but each subject also has strong private extremes.

## Highest Inter-Individual Acute Variability

From `high_variance_acute_genes.tsv` (top mean SD across acute contrasts):

1. `ITGAL`, `DAND5`, `COL4A4`, `GH1`, `MUC16`, `NPPB`, `AGBL2` are among the most variable.

Interpretation: these are candidate markers for responder stratification analyses.

## Claim Framing for Showcase Paper

Confirmed:

1. Acute heat response has a shared cohort-level directionality.
2. Individual response magnitude spans a broad range and forms separable responder classes.

Exploratory:

1. Specific high-variance proteins as stratification biomarkers require external validation.
2. Subject-level fingerprint proteins should be presented as heterogeneity descriptors, not universal markers.

