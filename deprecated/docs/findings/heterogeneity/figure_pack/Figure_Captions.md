# Figure Captions

## Figure 1

### Fig1A. Subject-Level Proteomic Response Map
Scatter of per-subject acute burden versus adaptation burden.  
Acute burden is the mean absolute protein delta across `PT1-PR1` and `PT2-PR2`; adaptation burden is the mean absolute protein delta across `PT2-PT1` and `PR2-PR1`. Points are colored by responder class and labeled by subject ID.

### Fig1B. Responder Class Composition
Bar chart of subject counts per responder class (`HighAcute_HighAdapt`, `HighAcute_LowAdapt`, `LowAcute_HighAdapt`, `LowAcute_LowAdapt`) derived from median splits of acute and adaptation burdens.

## Figure 2

### Fig2A. Individual Burden Profiles
Per-subject paired bars for acute burden and adaptation burden with connector lines to emphasize within-subject profile shape.

### Fig2B. Proteins with Highest Inter-Individual Acute Variability
Heatmap of subject-level acute mean deltas (mean of `PT1-PR1` and `PT2-PR2`) for the top variable proteins ranked by inter-individual standard deviation.

## Figure 3

### Fig3A. Volcano Plot for PT1-PR1
Differential abundance volcano for `PT1-PR1` using mean difference on the x-axis and `-log10(q)` on the y-axis. The dashed line marks `q=0.05`; highlighted points are proteins with `q<0.05`.

### Fig3B. Volcano Plot for PT2-PR2
Differential abundance volcano for `PT2-PR2` with the same visual encoding as Fig3A to enable direct acute-contrast comparison.

### Fig3C. Discovery Counts Across Comparisons
Grouped bars showing number of proteins at `q<0.05` and `q<0.10` for each comparison (`PT1-PR1`, `PT2-PR2`, `PT2-PT1`, `PR2-PR1`).

## Figure 4

### Fig4A. LOO Retention of Baseline Acute q<0.05 Proteins
Per-protein retention counts under leave-one-subject-out reruns for proteins significant at baseline (`q<0.05`) in acute contrasts. Y-axis reports number of LOO reruns (0-10) preserving `q<0.05`.

### Fig4B. Top-Rank Stability Under Leave-One-Subject-Out
Top-20 overlap counts by left-out subject for `PT1-PR1` and `PT2-PR2`, quantifying ranking stability of leading acute signals.

## Figure 5

### Fig5A. Burden Distribution by Responder Class
Class-stratified strip plot of acute and adaptation burdens, showing separation and spread within each responder class.

### Fig5B. Pathway Discovery Density (Restricted Universe)
Heatmap of pathway counts at `FDR<0.10` by comparison and collection from restricted-universe enrichment results.

## Figure 6

### Fig6A. Recurrent Subject-Fingerprint Proteins
Frequency plot of proteins appearing in the per-subject top-5 absolute deltas (separately for `pt1_pr1` and `pt2_pr2`), highlighting recurring versus private response features.

### Fig6B. Distribution of Top Variable Proteins Across Subjects
Horizontal boxplots of acute mean deltas for the highest inter-individual-variance proteins, summarizing spread, median shift, and outlier structure.

## Figure 7

### Fig7A. Subject Archetype Map on Latent Axes
Two-dimensional latent projection (`PC1`, `PC2`) of subjects derived from burden and module-trajectory features.  
Each point is a subject, colored by composite archetype label and split by median axis boundaries.

### Fig7B. Module Trajectory Score Heatmap
Subject-by-feature heatmap of module trajectory scores across contrasts (`pt1_pr1`, `pt2_pr2`, `pt2_pt1`, `pr2_pr1`) for heat-shock, immune, vascular, and ECM-remodeling modules.

## Figure 8

### Fig8A. Dominant PC1 Feature Loadings
Top absolute loadings contributing to latent axis `PC1`, indicating which burden/module features most strongly separate subjects along the primary heterogeneity axis.

### Fig8B. Physiology Mapping Heatmap
Pearson correlation heatmap between latent metrics (`PC1`, `PC2`, acute burden, adaptation burden) and physiological delta endpoints from the Dube phenotype sheet. Values are correlation coefficients.
