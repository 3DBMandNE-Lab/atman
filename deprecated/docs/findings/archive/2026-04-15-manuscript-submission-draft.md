# Reproducible and Heterogeneity-Aware Analysis of Human Heat-Stress Proteomics with Atman

## Introduction

Human heat adaptation has two simultaneous properties: a shared stress-response program and substantial person-to-person variation in magnitude and trajectory. In small cohorts, standard proteomics reporting often emphasizes cohort-average significance while under-reporting robustness and heterogeneity structure. We analyzed the Dube heat-stress Olink dataset to evaluate whether `Atman` can support both reproducible core inference and biologically meaningful subject-level structure from the same raw data.

The study includes 10 participants measured at four conditions (`PR1`, `PT1`, `PR2`, `PT2`), enabling paired acute and acclimation contrasts.

## Methods

### Data and pipeline

Raw NPX data were processed end-to-end with public `Atman` commands (`ingest`, `qc`, `matrix`, `fold-change`, `de`). The resulting sample sheet contained 43 total samples: 40 biological samples (10 per condition) and 3 technical controls.

### Differential abundance

Differential abundance used paired Student’s t-tests at subject level with Benjamini-Hochberg correction within each comparison family. Comparisons were `PT1-PR1`, `PT2-PR2`, `PT2-PT1`, and `PR2-PR1`; `min_pairs=5`.

### Robustness

We ran leave-one-subject-out (LOO) reruns (10 reruns), measured sign consistency, strict/relaxed hit retention (`q<0.05`, `q<0.10`), and top-20 overlap stability. We also tested sensitivity to `min_pairs` (`4`, `5`, `6`).

### Heterogeneity and latent structure

For each subject, we computed per-protein deltas for `pt1_pr1`, `pt2_pr2`, `pt2_pt1`, and `pr2_pr1`. We defined:

1. acute burden = mean absolute delta across `pt1_pr1` and `pt2_pr2`,
2. adaptation burden = mean absolute delta across `pt2_pt1` and `pr2_pr1`.

Responder classes were defined by median splits of acute and adaptation burdens (high/low x high/low). We then computed recurrent subject fingerprint proteins, high-variance proteins, latent axes from burden/module features, and physiology mapping (Pearson correlation + BH correction).

## Results

### 1) Core cohort-level signal

The pipeline processed 122,135 measurements across 2,943 assays; 3,725 rows were QC-masked and 118,410 passed QC. Differential abundance was computed for 2,938 proteins per comparison.

Discovery burden was strongly concentrated in acute contrasts:

1. `PT1-PR1`: 36 proteins at `q<0.05` (62 at `q<0.10`)
2. `PT2-PR2`: 37 proteins at `q<0.05` (277 at `q<0.10`)
3. `PT2-PT1`: 7 proteins at `q<0.05` (8 at `q<0.10`)
4. `PR2-PR1`: 2 proteins at `q<0.05` (7 at `q<0.10`)

Canonical heat-stress proteins (`HSPB1`, `HSPA1A`, `DNAJB1`, `HSPG2`) showed directionally positive acute behavior in one or both acute contrasts, consistent with expected biology.

### 2) Robustness profile

LOO sign stability was high:

1. `PT1-PR1`: mean sign-match rate `0.9358`
2. `PT2-PR2`: mean sign-match rate `0.9398`

Strict-hit retention was lower for many non-canonical proteins. Among baseline acute `q<0.05` proteins:

1. `PT1-PR1`: 11/36 retained `q<0.05` in all 10 LOO reruns
2. `PT2-PR2`: 5/37 retained `q<0.05` in all 10 LOO reruns

Relaxed-threshold retention was stronger:

1. `PT1-PR1`: 30/36 retained `q<0.10` in at least 8/10 reruns
2. `PT2-PR2`: 33/37 retained `q<0.10` in at least 8/10 reruns

Top-20 rank stability was moderate:

1. `PT1-PR1`: mean overlap `16.0/20`, mean Jaccard `0.673`
2. `PT2-PR2`: mean overlap `14.5/20`, mean Jaccard `0.589`

Changing `min_pairs` from 4 to 6 did not change discovery counts.

### 3) Inter-individual heterogeneity

Acute burden ranged from `0.387198` to `0.612802` (mean `0.488900`, SD `0.075997`).  
Adaptation burden ranged from `0.352999` to `0.628138` (mean `0.467690`, SD `0.088434`).  
Acute and adaptation burdens were positively correlated (`r=0.631629`).

Responder classes were distributed across all four states:

1. `HighAcute_HighAdapt`: 3
2. `HighAcute_LowAdapt`: 2
3. `LowAcute_HighAdapt`: 2
4. `LowAcute_LowAdapt`: 3

This indicates structured human variability rather than a single pooled adaptation pattern.

### 4) Recurrent and high-variance proteins

Across per-subject top-5 proteins:

1. `DAND5` appeared in `pt2_pr2` top-5 for 6/10 subjects.
2. `GH1` appeared in top-5 for 5/10 subjects in both `pt1_pr1` and `pt2_pr2`.
3. `ITGAL` and `AGBL2` appeared repeatedly in acute top-5 sets.

Top acute inter-individual variance proteins included `ITGAL`, `DAND5`, `COL4A4`, `GH1`, `MUC16`, `NPPB`, and `AGBL2`.

### 5) Latent axes and physiology mapping

Latent-axis analysis separated subjects into multiple archetypes. The strongest raw phenotype association was `PC1` with acclimation sweat-rate delta (`acc_dLSR`, `r=0.880`, `p=0.0039`, `n=8`). No metric-phenotype association met BH `q<0.10` after correction.

## Discussion

The dataset supports a two-layer interpretation. First, there is a robust cohort-level acute heat response. Second, and biologically central for human studies, response magnitude and adaptation behavior vary substantially across individuals and form structured responder classes. This means the scientifically faithful narrative is not “one average responder,” but “shared acute core plus multiple individual adaptation programs.”

Robustness analyses sharpen this interpretation. Directional biology is stable under participant perturbation, while many non-canonical strict-threshold calls are sensitive at `n=10`. This supports reporting stable directionality and heterogeneity architecture as primary findings, with strict single-protein claims treated more conservatively.

## Limitations

The cohort is small (`n=10`), and some physiology variables have missing values, reducing effective sample size for phenotype mapping. The physiology correlations are therefore informative but not confirmatory after multiple-testing correction. External replication is required before clinical stratification claims.

## Data and Figure Availability

Master index:

1. [2026-04-15-showcase-results-dossier.md](/Users/kevinjoseph/Cursor/atman/docs/findings/2026-04-15-showcase-results-dossier.md)

Figure pack:

1. [figure_pack](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/figure_pack)

Captions:

1. [Figure_Captions.md](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/figure_pack/Figure_Captions.md)

