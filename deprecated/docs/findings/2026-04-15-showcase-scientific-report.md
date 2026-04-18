# Atman Showcase Scientific Report

Date: 2026-04-15  
Project: Human heat-stress proteomics showcase using `Atman`  
Dataset: Dube et al. Olink NPX heat stress / acclimation cohort

## 1. Background

Most proteomics reports on small human cohorts stop at cohort-average significance lists.  
That misses the main biological reality in humans: response heterogeneity.

This report has two goals:

1. verify that `Atman` produces a reproducible and statistically disciplined core analysis from raw NPX inputs,
2. extract biologically meaningful inter-individual response structure (not just group means).

## 2. Research Questions

### RQ1. Can the full raw-to-results pipeline run reproducibly and cleanly?

### RQ2. What is the dominant cohort-level biological signal?

### RQ3. How robust are differential-abundance findings to subject perturbation?

### RQ4. How much inter-individual variation exists, and is it structured?

### RQ5. Do latent response axes map to physiology in a meaningful way?

## 3. Cohort and Samples

From `/tmp/karna_showcase/samples.tsv`:

1. Total samples: `43`
2. Biological samples: `40`
3. Controls: `3`
4. Participants: `10`
5. Balanced biological design:
   - `PR1` = 10
   - `PT1` = 10
   - `PR2` = 10
   - `PT2` = 10

Condition definitions:

1. `PR1`: normothermic, pre-acclimation
2. `PT1`: hyperthermic, pre-acclimation
3. `PR2`: normothermic, post-acclimation
4. `PT2`: hyperthermic, post-acclimation

## 4. Analysis Pipeline

Executed commands (existing public CLI only):

1. `atman ingest`
2. `atman qc`
3. `atman matrix`
4. `atman fold-change`
5. `atman de --test paired-t --paired-by participant --min-pairs 5`

Run workspace: `/tmp/karna_showcase`

Pipeline outcomes:

1. Ingested rows: `122,135`
2. Assays: `2,943`
3. QC-masked rows: `3,725`
4. QC-passed rows: `118,410`
5. DE rows: `11,752` (all computed, none skipped)

## 5. Statistical Methods

### 5.1 Differential abundance

1. Paired Student’s t-test at subject level.
2. BH-FDR applied within each comparison family.
3. Comparisons:
   - `PT1-PR1`
   - `PT2-PR2`
   - `PT2-PT1`
   - `PR2-PR1`

### 5.2 Robustness

1. Leave-one-subject-out reruns (`n=10` reruns).
2. Per-protein direction stability.
3. Strict/relaxed hit retention (`q<0.05`, `q<0.10`).
4. Top-20 ranking overlap and Jaccard stability.
5. Threshold sweep: `min_pairs = 4, 5, 6`.

### 5.3 Heterogeneity

Per-subject per-protein deltas were computed for:

1. `pt1_pr1`
2. `pt2_pr2`
3. `pt2_pt1`
4. `pr2_pr1`

Derived subject metrics:

1. Acute burden = mean absolute delta across `pt1_pr1`, `pt2_pr2`
2. Adaptation burden = mean absolute delta across `pt2_pt1`, `pr2_pr1`
3. Responder classes = median split of acute and adaptation burdens

Extended layers:

1. Recurrent top-5 subject fingerprint proteins
2. Inter-individual variance ranking
3. Latent axes from burden + module trajectory features
4. Physiology mapping by Pearson correlation + BH correction

## 6. Results by Research Question

### 6.1 RQ1: Reproducible pipeline execution

Answer: **Yes**.

Evidence:

1. Full raw-data pipeline executed end-to-end without manual intervention.
2. Core output tables generated consistently (`de_results.tsv`, `de_report.tsv`, matrix/fold-change outputs).
3. No new public CLI surface was required for advanced analysis layers.

### 6.2 RQ2: Dominant cohort-level signal

Answer: **Acute heat contrasts dominate discovery burden**.

Per-comparison discoveries (2,938 proteins/comparison):

1. `PT1-PR1`: `q<0.05 = 36`, `q<0.10 = 62`
2. `PT2-PR2`: `q<0.05 = 37`, `q<0.10 = 277`
3. `PT2-PT1`: `q<0.05 = 7`, `q<0.10 = 8`
4. `PR2-PR1`: `q<0.05 = 2`, `q<0.10 = 7`

Interpretation:

1. Acute hyperthermia transitions (`PR -> PT`) carry most of the differential-abundance signal.
2. Acclimation contrasts are weaker in strict multiple-testing space.

Canonical heat-stress markers (directional check):

1. `HSPB1`, `HSPA1A`, `DNAJB1`, `HSPG2` show positive acute-direction shifts in one or both acute contrasts.
2. Directional coherence is present even where strict q-threshold is not crossed for every marker.

### 6.3 RQ3: Robustness to subject perturbation

Answer: **Core directionality is stable; strict single-protein calls are partly fragile**.

Sign stability (mean over all proteins):

1. `PT1-PR1`: `0.9358`
2. `PT2-PR2`: `0.9398`

LOO retention for baseline acute `q<0.05` proteins:

1. `PT1-PR1` baseline 36:
   - retained `q<0.05` in all 10 reruns: `11`
   - retained `q<0.10` in at least 8/10 reruns: `30`
2. `PT2-PR2` baseline 37:
   - retained `q<0.05` in all 10 reruns: `5`
   - retained `q<0.10` in at least 8/10 reruns: `33`

Top-20 rank stability:

1. `PT1-PR1`: mean overlap `16.0/20`, mean Jaccard `0.673`
2. `PT2-PR2`: mean overlap `14.5/20`, mean Jaccard `0.589`

Threshold sensitivity:

1. `min_pairs` (`4`, `5`, `6`) produced identical discovery counts for all comparisons.

Interpretation:

1. Cohort-level acute biology is robust.
2. Many non-canonical strict-hits near threshold are sensitive at this cohort size.

### 6.4 RQ4: Inter-individual variation and structure

Answer: **Variation is large and structured into responder programs**.

Acute burden:

1. Range: `0.387198` to `0.612802`
2. Mean: `0.488900`
3. SD: `0.075997`

Adaptation burden:

1. Range: `0.352999` to `0.628138`
2. Mean: `0.467690`
3. SD: `0.088434`

Coupling:

1. Acute vs adaptation burden: `r = 0.631629`

Responder class counts:

1. `HighAcute_HighAdapt`: `3`
2. `HighAcute_LowAdapt`: `2`
3. `LowAcute_HighAdapt`: `2`
4. `LowAcute_LowAdapt`: `3`

Recurrent top-5 subject fingerprint proteins:

1. `DAND5` (pt2_pr2): `6/10` subjects
2. `GH1` (pt1_pr1): `5/10`
3. `GH1` (pt2_pr2): `5/10`
4. `ITGAL` (pt1_pr1): `4/10`
5. `AGBL2` (pt1_pr1): `4/10`

Highest inter-individual acute-variance proteins include:

1. `ITGAL`
2. `DAND5`
3. `COL4A4`
4. `GH1`
5. `MUC16`
6. `NPPB`
7. `AGBL2`

Interpretation:

1. Human response is not a single pattern with noise.
2. There is a shared acute core plus distinct individual adaptation programs.

### 6.5 RQ5: Latent axes and physiology mapping

Answer: **Suggestive structure, not definitive after correction**.

Strongest raw mapping signal:

1. `PC1` vs `acc_dLSR`: `r = 0.880`, `p = 0.0039`, `n = 8`

Multiple-testing context:

1. No metric-phenotype pair reached BH `q<0.10` in this cohort.

Interpretation:

1. The latent structure appears biologically aligned with sweat-related adaptation.
2. Evidence is hypothesis-generating rather than confirmatory at current sample size.

## 7. Integrated Biological Interpretation

### 7.1 What is strongly supported

1. Acute heat induces a reproducible proteomic response at cohort level.
2. Canonical heat-stress directionality persists under perturbation.
3. Inter-individual response range is a major biological feature, not a nuisance term.

### 7.2 What is plausible but not yet final

1. Specific high-variance proteins as stratification biomarkers.
2. Specific latent-axis physiology links as adaptation mechanisms.

### 7.3 What should be kept explicitly exploratory

1. Non-canonical strict `q<0.05` single-protein claims near threshold.
2. Any clinical-personalization claim without external validation.

## 8. Limitations

1. `n=10` limits power for high-dimensional correction-heavy inference.
2. Some physiology variables have missing values, reducing effective `n`.
3. Single-cohort setting limits external generalizability.

## 9. Reporting Guidance for Manuscript

### 9.1 Primary narrative

1. Reproducible engine + robust acute biology.
2. Human heterogeneity as a core biological finding.
3. Structured uncertainty treatment (robustness + claim tiers).

### 9.2 Claim register

Strong:

1. Acute response exists and is robust.
2. Heterogeneity is large and structured.

Moderate:

1. Recurrent variable proteins may define responder states.
2. Latent axes may capture adaptation dimensions.

Exploratory:

1. Fine-grained single-protein and physiology links after full correction.

## 10. Artifact Map

### 10.1 Core summaries

1. [2026-04-14-showcase-week1-execution.md](/Users/kevinjoseph/Cursor/atman/docs/findings/2026-04-14-showcase-week1-execution.md)
2. [2026-04-14-robustness-summary.md](/Users/kevinjoseph/Cursor/atman/docs/findings/robustness/2026-04-14-robustness-summary.md)
3. [2026-04-14-heterogeneity-summary.md](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/2026-04-14-heterogeneity-summary.md)

### 10.2 Robustness tables

1. [loo_gene_stability.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/robustness/loo_gene_stability.tsv)
2. [rank_stability.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/robustness/rank_stability.tsv)
3. [threshold_sensitivity.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/robustness/threshold_sensitivity.tsv)

### 10.3 Heterogeneity and latent tables

1. [per_subject_gene_deltas.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/per_subject_gene_deltas.tsv)
2. [subject_response_burden.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/subject_response_burden.tsv)
3. [responder_clusters.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/responder_clusters.tsv)
4. [subject_top_proteins.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/subject_top_proteins.tsv)
5. [high_variance_acute_genes.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/high_variance_acute_genes.tsv)
6. [module_trajectory_scores.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/module_trajectory_scores.tsv)
7. [latent_axes.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/latent_axes.tsv)
8. [subject_archetypes.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/subject_archetypes.tsv)
9. [physiology_mapping.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/physiology_mapping.tsv)

### 10.4 Figures

Figure directory:

[figure_pack](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/figure_pack)

Caption file:

[Figure_Captions.md](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/figure_pack/Figure_Captions.md)

