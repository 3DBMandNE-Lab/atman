# Reproducible and Heterogeneity-Aware Analysis of Human Heat-Stress Proteomics with Atman

## Introduction

Human heat-stress biology combines a shared stress-response core with substantial inter-individual variability in response magnitude and adaptation trajectory. Most small-cohort proteomics analyses emphasize cohort-average significance while underreporting robustness and heterogeneity structure. This creates a recurrent gap between statistical detection and biologically realistic interpretation in human datasets.

We used `Atman` as a local-first analysis engine to process raw Olink NPX data from the Dube heat-stress/acclimation cohort and asked three linked questions: (1) whether the full raw-to-result pipeline is reproducible, (2) what cohort-level response is strongest, and (3) whether inter-individual variation is structured enough to support responder stratification rather than a single pooled narrative.

The dataset contains 10 participants measured at four conditions (`PR1`, `PT1`, `PR2`, `PT2`) plus technical controls, enabling both paired group-level inference and subject-resolved trajectory analysis.

## Methods

### Data and design

The analysis used two raw Olink NPX files from the Dube dataset and corresponding physiological measurements. The generated sample sheet contained 43 total samples: 40 biological samples (10 each at `PR1`, `PT1`, `PR2`, `PT2`) and 3 control samples.

### Core proteomics pipeline

The following `Atman` commands were run from raw inputs:

1. `ingest`
2. `qc`
3. `matrix`
4. `fold-change`
5. `de` (`paired-t`, paired by participant, `min_pairs=5`)

The run produced 122,135 ingested measurement rows, with 3,725 QC-masked rows and 118,410 QC-passed rows.

### Differential abundance statistics

Differential abundance used paired Student’s t-tests at subject level with BH-FDR correction within each comparison family (`PT1-PR1`, `PT2-PR2`, `PT2-PT1`, `PR2-PR1`). Each family included 2,938 proteins.

### Robustness analysis

Robustness was assessed by leave-one-subject-out reruns (10 reruns total), threshold sweep (`min_pairs=4,5,6`), and top-20 ranking overlap metrics (overlap count and Jaccard index).

### Heterogeneity analysis

Per-subject protein deltas were computed for `pt1_pr1`, `pt2_pr2`, `pt2_pt1`, and `pr2_pr1`. Subject-level metrics were defined as:

1. Acute burden: mean absolute delta across `pt1_pr1` and `pt2_pr2`
2. Adaptation burden: mean absolute delta across `pt2_pt1` and `pr2_pr1`

Responder classes were defined by median split of acute and adaptation burdens (high/low x high/low). Additional layers included recurrent top-5 subject fingerprint proteins, high inter-individual variance proteins, latent-axis decomposition from burden/module features, and phenotype correlation mapping.

## Results

### Cohort-level signal is strongest in acute heat contrasts

Differential-abundance discovery counts showed a clear acute dominance:

1. `PT1-PR1`: 36 proteins at `q<0.05` (62 at `q<0.10`)
2. `PT2-PR2`: 37 proteins at `q<0.05` (277 at `q<0.10`)
3. `PT2-PT1`: 7 proteins at `q<0.05` (8 at `q<0.10`)
4. `PR2-PR1`: 2 proteins at `q<0.05` (7 at `q<0.10`)

Canonical heat-stress proteins (`HSPB1`, `HSPA1A`, `DNAJB1`, `HSPG2`) showed directionally positive acute responses in one or both acute contrasts, supporting biological coherence of the acute signature.

### Acute directionality is robust, while strict single-protein significance is partly fragile

LOO analyses showed high per-protein sign stability in acute contrasts:

1. `PT1-PR1`: mean sign-match rate 0.9358
2. `PT2-PR2`: mean sign-match rate 0.9398

However, strict hit retention was lower for many non-canonical proteins. Among baseline acute `q<0.05` proteins:

1. `PT1-PR1`: 11/36 remained `q<0.05` in all 10 LOO reruns
2. `PT2-PR2`: 5/37 remained `q<0.05` in all 10 LOO reruns

At a relaxed threshold (`q<0.10`) and requiring retention in at least 8/10 reruns, stability increased (30/36 for `PT1-PR1`; 33/37 for `PT2-PR2`).

Top-20 ranking stability was moderate:

1. `PT1-PR1`: mean overlap 16.0/20, mean Jaccard 0.673
2. `PT2-PR2`: mean overlap 14.5/20, mean Jaccard 0.589

Threshold sweep results were unchanged across `min_pairs=4,5,6`, indicating that headline counts were not sensitive to this parameter range.

### Inter-individual response heterogeneity is substantial and structured

Subject burden distributions demonstrated clear range:

1. Acute burden: 0.387198 to 0.612802 (mean 0.488900, SD 0.075997)
2. Adaptation burden: 0.352999 to 0.628138 (mean 0.467690, SD 0.088434)
3. Acute-adaptation correlation: `r=0.631629`

Responder classes were balanced across four states:

1. HighAcute_HighAdapt: 3 subjects
2. HighAcute_LowAdapt: 2 subjects
3. LowAcute_HighAdapt: 2 subjects
4. LowAcute_LowAdapt: 3 subjects

This pattern indicates that pooled cohort summaries hide biologically meaningful individual programs.

### Recurrent and high-variance proteins define candidate stratification axes

Across per-subject top-5 acute proteins, recurrence was observed for:

1. `DAND5` (pt2_pr2): 6/10 subjects
2. `GH1` (pt1_pr1): 5/10 subjects
3. `GH1` (pt2_pr2): 5/10 subjects
4. `ITGAL` and `AGBL2` (pt1_pr1): 4/10 subjects each

Highest acute inter-individual variance proteins included `ITGAL`, `DAND5`, `COL4A4`, `GH1`, `MUC16`, `NPPB`, and `AGBL2`, supporting a stratification-oriented interpretation rather than a single canonical responder model.

### Latent axes capture structure; physiology mapping is suggestive

Latent-axis analysis separated subjects into multiple archetypes consistent with burden and module-trajectory features. The strongest raw phenotype association was between `PC1` and acclimation sweat-rate delta (`acc_dLSR`, `r=0.880`, `p=0.0039`, `n=8`). After BH correction across tested phenotype links, no association met `q<0.10`.

Thus, physiology mapping provides biologically plausible structure, but not correction-robust endpoint-specific claims in this cohort size.

## Discussion

This analysis supports a dual-layer interpretation of human heat-stress proteomics.

First, a cohort-level acute response core is reproducible and robust in direction under perturbation. This validates both the biological signal and the computational reliability of the `Atman` pipeline.

Second, and more importantly for human biology, substantial structured heterogeneity is present. Subjects distribute across distinct acute/adaptation response classes, and recurrent high-variance proteins suggest that adaptation is composed of multiple individual programs rather than one shared trajectory.

The robustness analysis clarifies how to report this responsibly: acute directionality is strong; many non-canonical strict-threshold proteins are sensitive to subject removal; relaxed and direction-consistent summaries are more stable in `n=10`.

For showcase framing, the strongest contribution is not merely “detecting significant proteins,” but demonstrating that a reproducible engine can recover both shared biology and person-level architecture from the same raw dataset.

## Limitations

1. Small cohort size limits high-dimensional multiple-testing power.
2. Some physiology variables have missingness, reducing effective sample size for mapping.
3. Results are from a single cohort and require external replication for generalization.

## Data Products and Figures

Core outputs are documented in:

1. [2026-04-15-showcase-results-dossier.md](/Users/kevinjoseph/Cursor/atman/docs/findings/2026-04-15-showcase-results-dossier.md)

Manuscript panels are in:

1. [figure_pack](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/figure_pack)

Captions are in:

1. [Figure_Captions.md](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/figure_pack/Figure_Captions.md)

