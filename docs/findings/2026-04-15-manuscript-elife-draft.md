# Reproducible and heterogeneity-aware analysis of human heat-stress proteomics with karnaProteome

## Abstract

Human heat response in this cohort follows a two-layer architecture: a robust shared acute program plus structured inter-individual adaptation variability. Using karnaProteome to reanalyze the Dube heat-stress/acclimation Olink NPX dataset from raw inputs, we tested whether one deterministic workflow can recover both layers. Differential abundance showed dominant acute signal (PT1-PR1 and PT2-PR2) with weaker acclimation contrasts (PT2-PT1 and PR2-PR1), while leave-one-subject-out reruns preserved high directional stability but revealed threshold fragility for borderline single-protein calls. Subject-level analyses identified broad, organized heterogeneity with all responder classes represented and recurrent high-variation proteins (including DAND5, GH1, ITGAL, and AGBL2). Latent-axis to physiology mapping yielded suggestive raw associations but no BH-significant links after correction. Together, these findings show that heterogeneity-aware analysis changes biological interpretation beyond pooled significance lists and should be treated as a standard output for small-cohort human proteomics.

## Introduction

Human heat acclimation induces coordinated thermoregulatory, cardiovascular, and cellular adaptations, but substantial inter-individual variability is also expected in acclimation trajectories and stress tolerance [1,5,6]. In practice, small human proteomics cohorts are often summarized as pooled significance lists, which can compress biologically meaningful subject-level structure into a single group-average narrative. This risk is amplified in high-dimensional assays where technical and batch structure can distort weak biological effects if not controlled [12]. The Dube dataset provides repeated within-subject plasma proteomics across pre/post-acclimation and normothermic/hyperthermic states and is therefore well suited for separating cohort-average and subject-specific structure [2]. The source publication established the dataset as a reusable resource and reported limited high-confidence physiology-proteomics coupling at this sample size [2], making it an appropriate benchmark for robust heterogeneity-aware reanalysis.

Here we use karnaProteome as a deterministic analysis engine and ask three linked questions: can the full raw-data workflow run reproducibly; what signal is strongest at cohort level; and can subject-level structure be quantified in a way that changes biological interpretation.

## Results

### End-to-end execution from raw files is stable and internally consistent

The complete workflow (`ingest`, `qc`, `matrix`, `fold-change`, `de`) executed without intervention from raw NPX inputs to final differential-abundance tables. The run processed 122,135 measurements across 2,943 assays. QC masking removed 3,725 rows from inferential abundance values while preserving raw values in the long table, leaving 118,410 usable rows after masking logic. The sample table retained the expected experimental structure (40 biological samples, 3 technical controls) and preserved balanced biological design (10 participants x 4 conditions). This confirms that the same engine can be used for both strict reproduction and downstream inferential analysis without pipeline branching, consistent with reproducibility expectations established in high-throughput expression studies [4].

### Cohort-level signal is concentrated in acute transitions, not acclimation contrasts

Differential abundance was evaluated for 2,938 proteins in each comparison family. Signal concentration was asymmetric: acute contrasts dominated discovery burden, whereas acclimation contrasts were sparse. Specifically, PT1-PR1 yielded 36 proteins at q<0.05 (62 at q<0.10), and PT2-PR2 yielded 37 proteins at q<0.05 (277 at q<0.10). In contrast, PT2-PT1 yielded 7 proteins at q<0.05 (8 at q<0.10) and PR2-PR1 yielded 2 proteins at q<0.05 (7 at q<0.10). The primary statistical structure of the dataset is therefore acute heat response rather than broad acclimation-at-rest shift (Fig 3A-C). Notably, acclimation-associated signal was stronger under thermal challenge (PT2-PT1) than at rest (PR2-PR1), consistent with adaptation effects becoming more detectable during provocation than in baseline state.

Canonical heat-stress proteins were directionally concordant with this interpretation. In PT2-PR2, for example, HSPA1A showed mean_diff 0.71389 (q=0.09625), DNAJB1 mean_diff 0.85613 (q=0.15731), and HSPB1 mean_diff 1.28053 (q=0.19879). Although not all reached strict FDR thresholds, directional coherence across canonical markers supports biological validity of the acute-response axis. We treat this as a directional concordance check, not as inferential proof for each marker independently.

### Robustness analysis separates stable directional biology from threshold-sensitive fine structure

LOO reruns (10 total) were used to probe dependence on individual participants. Across all proteins, directionality in acute contrasts was highly stable: mean sign-match rate was 0.9358 for PT1-PR1 and 0.9398 for PT2-PR2. This indicates that the bulk directional geometry of the acute response is not driven by single-subject artifacts.

However, strict-threshold retention was substantially lower for many individual proteins. Of baseline acute q<0.05 hits, 11/36 (PT1-PR1) and 5/37 (PT2-PR2) remained q<0.05 in all LOO reruns. At a less brittle threshold (q<0.10, >=8/10 reruns), retention increased to 30/36 and 33/37, respectively. This pattern shows that the dataset supports strong directional and moderate-threshold reproducibility, while strict single-protein claims are more sensitive near the significance boundary.

Rank robustness showed similar behavior. Top-20 overlap with baseline averaged 16.0 proteins for PT1-PR1 and 14.5 proteins for PT2-PR2 (mean Jaccard 0.673 and 0.589). This is consistent with moderate ranking stability: a stable core plus an unstable margin. Threshold-sweep analysis (min_pairs 4, 5, 6) did not change headline counts in any family, indicating that these conclusions are not artifacts of one min_pairs setting (Fig 4A-B).

To test dependence on model choice, we performed a mixed-model sensitivity analysis for acute proteins with baseline q<0.10 using random-intercept models (abundance ~ condition + (1|subject)). Convergence was achieved for all tested proteins (62/62 in PT1-PR1, 277/277 in PT2-PR2). Effect-direction concordance with paired t-tests was 1.00 in both contrasts, and rank correspondence of absolute effect sizes was also 1.00 (Spearman). Given the paired repeated-measures structure and n=10, this near-complete concordance is expected and should be interpreted as consistency rather than an independent orthogonal validation.

### Inter-individual variation is large and organized, not residual noise

Subject-level burden metrics showed broad range on both axes. Acute burden ranged from 0.387198 to 0.612802 (mean 0.488900), while adaptation burden ranged from 0.352999 to 0.628138 (mean 0.467690). The axes were positively correlated (r=0.631629), indicating partial coupling between immediate stress-response amplitude and longer adaptation signal, but not complete redundancy (Fig 1A, Fig 2A).

Because burden is an aggregation statistic (mean absolute delta) rather than a latent-parameter estimator, we tested alternative burden definitions. Subject rankings were highly preserved under both a top-500-variance burden and an inverse-variance weighted burden (acute Spearman 0.988 and 0.952 vs baseline; adaptation Spearman 0.915 and 0.988 vs baseline), indicating that the heterogeneity pattern is not an artifact of one burden formula.

Subjects populated all four responder bins under median split: HighAcute_HighAdapt (3), HighAcute_LowAdapt (2), LowAcute_HighAdapt (2), and LowAcute_LowAdapt (3). We use these bins as descriptive partitions of structured continuous variation, not as a formal clustering proof with uncertainty intervals (Fig 1B, Fig 5).

To evaluate whether observed heterogeneity exceeded a noise-preserving null, we permuted subject labels within each protein (500 iterations) and recomputed burdens. Observed between-subject dispersion was substantially above null expectation (acute SD: observed 0.075997 vs null mean 0.006365, empirical p=0.001996; adaptation SD: observed 0.088434 vs null mean 0.006893, empirical p=0.001996). Because this equals the permutation-resolution floor (1/501), we report these as p<0.002 at this iteration depth.

### Recurrent fingerprint proteins and high-variance proteins identify candidate stratification axes

Per-subject top-5 acute protein fingerprints were non-random across participants. DAND5 appeared in the PT2-PR2 top-5 for 6/10 subjects. GH1 appeared in top-5 for 5/10 subjects in both acute contrasts. ITGAL and AGBL2 were recurrent in PT1-PR1 top-5 sets. These repeated appearances suggest potential responder-axis features rather than idiosyncratic single-subject outliers (Fig 6A). Biologically, these proteins outline plausible axes of individualized adaptation: DAND5 points to differential engagement of BMP/TGF-beta regulatory tone, GH1 to inter-individual endocrine-metabolic stress handling, ITGAL to variable leukocyte adhesion/immune trafficking, and AGBL2 to cytoskeletal remodeling programs that can shape cellular stress tolerance. The same interpretation extends to high-variance proteins NPPB and COL4A4 (Fig 6B), which are compatible with variability in cardio-vascular load signaling and extracellular-matrix/basement-membrane remodeling under repeated heat exposure.

Protein-level variance ranking reinforced this view. The highest acute inter-individual variance proteins included ITGAL, DAND5, COL4A4, GH1, MUC16, NPPB, and AGBL2 (Fig 2B, Fig 6B). These proteins define a plausible high-heterogeneity feature set for follow-up stratification analyses in larger cohorts.

### Latent axes recover subject archetypes; phenotype coupling remains exploratory after correction

Latent-axis decomposition (burden plus module-trajectory features) separated subjects into multiple archetype regions and preserved the non-uniform responder geometry observed in burden space (Fig 7A-B). PC1 and PC2 explained 29.32% and 20.82% of latent-feature variance, respectively. PC1 loadings indicated that burden/module features materially contributed to this axis geometry (Fig 8A). Physiology mapping provided a biologically suggestive but statistically constrained signal pattern. The strongest raw association was PC1 with acclimation sweat-rate delta (acc_dLSR, r=0.880, p=0.0039, n=8), with additional raw associations in related thermoregulatory/circulatory endpoints (Fig 8B). After BH correction within metric families, no metric-endpoint association reached q<0.10. At n=8 for key acclimation endpoints, these raw correlations are sensitive to leverage by individual observations and should be treated as exploratory signals requiring external validation.

The appropriate interpretation is therefore structural rather than confirmatory: latent axes encode candidate adaptation dimensions, but endpoint-specific mechanistic claims remain underpowered in this cohort.

## Discussion

This study shows that the Dube heat dataset is best interpreted through a two-layer framework: a reproducible cohort-level acute response and a structured person-level adaptation architecture. Limiting interpretation to pooled significance counts would capture the first layer while obscuring the second, and would therefore provide an incomplete account of human response biology.

The robustness profile is central to this interpretation. High proteome-wide sign stability under LOO confirms that acute directional biology is not a single-subject artifact. At the same time, lower all-rerun retention of strict q<0.05 non-canonical hits demonstrates that many near-threshold single-protein claims are fragile at n=10. This is not a contradiction; it is the expected geometry of small-cohort omics data with a stable biological core and a threshold-sensitive margin. In practical terms, this supports confidence in direction-level and pattern-level claims while imposing caution on strict per-protein ranking claims.

The heterogeneity layer changes the biological narrative materially. Burden ranges, full occupancy of responder bins, recurrent high-variation proteins, and permutation-supported excess dispersion indicate that adaptation is not simply weaker or stronger versions of one cohort-average program. Instead, the data are best described as multi-axis heterogeneity over a shared acute core. The recurrent fingerprint proteins provide biologically plausible axis candidates, including BMP/TGF-beta regulation (DAND5), endocrine stress signaling (GH1), immune trafficking (ITGAL), and vascular/ECM load biology (NPPB, COL4A4). This aligns with the broader heat-acclimation literature, where common physiological adaptation coexists with substantial person-level response spread [1,6]. That finding has direct implications for downstream study design: replication should prioritize stratified validation and responder-axis modeling, not only pooled differential-abundance replication.

Latent-axis and physiology mapping results reinforce this point. The presence of strong raw links (for example PC1 with acclimation sweat-rate change) alongside non-significant BH-adjusted outcomes suggests biologically informative structure with insufficient inferential power for endpoint-level confirmation. The source dataset report also emphasized this endpoint-level power limitation [2]. The next rational step is not to discard these axes, but to treat them as candidate a priori hypotheses for external cohorts with larger effective n and lower phenotype missingness.

Methodologically, the key contribution of karnaProteome in this context is not just that it reproduces canonical outputs. The stronger demonstration is that the same deterministic pipeline can support both conservative cohort inference and richer heterogeneity-aware analyses without ad hoc toolchain divergence. For human translational proteomics, that unification of reproducibility and individual-architecture analysis is a practical advantage.

### Limitations

This analysis is constrained by small sample size (n=10 overall; n=8 for several acclimation physiology links), which limits endpoint-level inferential power and increases leverage sensitivity for correlation estimates. All findings derive from one cohort and one platform (Olink NPX), so cross-cohort and cross-platform generalizability remains untested. Module trajectory features are biologically curated and useful for structure discovery, but their module definitions are not yet externally benchmarked as formal prognostic signatures. karnaProteome now includes an optional moderated variance-shrinkage differential mode (limma-style empirical Bayes logic) for sensitivity analysis, but larger follow-up cohorts remain necessary for stable endpoint-level inference and external validation [15]. Phenotype missingness and multiple-testing burden further limit confirmatory physiology coupling in this dataset.

Overall, the evidence supports three conclusions. First, acute heat response is robust and reproducible at cohort level. Second, inter-individual variation is large, structured, and biologically relevant. Third, physiology-coupled latent adaptation axes are plausible but remain exploratory under current multiple-testing constraints.

## Ethics statement

This work is a secondary analysis of an already published, de-identified public dataset [2]. No new participant recruitment, intervention, or re-consent was performed in this study. Original ethics approval and consent procedures are described in the source publication [2]. Data use follows the access and reuse terms provided with the public dataset release [2].

## Author contributions

Author contributions to be finalized at submission using CRediT taxonomy.

## Funding

Funding statement to be completed at submission.

## Competing interests

Competing interests declaration to be completed at submission.

## Materials and methods

### Dataset and sample structure

The analysis used the Dube heat-stress/acclimation Olink Explore NPX dataset [2,7] (mirrored in the analysis workspace under `example_data/dube_heat_2023` for reproducibility; archived repository DOI to be added at submission). Biological sampling was balanced across four within-subject states: PR1, PT1, PR2, PT2 (10 samples each; 40 biological samples total), with three additional technical control samples (CONTROL_SAMPLE identifiers). Sample IDs were parsed using the Dube convention `SSNA-<subject>-<condition>`.

### Core proteomics preprocessing

Raw NPX files were processed with the public `karnaProteome` CLI:

1. `ingest` to produce canonical long-format `measurements.tsv`, `proteins.tsv`, and `samples.tsv`;
2. `qc` with Dube rule (mask rows where sample-level or assay-level QC warning was not PASS);
3. `matrix` to generate panel-wise wide matrices;
4. `fold-change` to compute contrast fold changes;
5. `de` for paired differential abundance.

The public release command surface also includes heterogeneity-oriented operations: `asymmetry` (paired contrast asymmetry metrics), `robustness` (LOO sign/rank stability summaries from DE outputs), and `module-trajectory` (module score extraction from user-supplied `modules.tsv` definitions and per-subject delta tables).

Control samples were retained in intermediate files for faithful reproduction behavior but excluded from differential abundance and subject-level biological burden calculations.

### Differential abundance model

For each protein and each comparison (PT1-PR1, PT2-PR2, PT2-PT1, PR2-PR1), we performed paired Student’s t-tests at subject level:

1. paired differences were computed as `d_i = a_i - b_i` over subjects present in both conditions;
2. tests were run only when the number of valid pairs met `min_pairs` (`min_pairs=5` in main analysis);
3. two-sided p-values were computed from the t-distribution with `df=n_pairs-1`;
4. BH-FDR correction was applied within each comparison family across all proteins tested in that comparison [3,11].

`de` supports two models: `paired-t` (default) and `moderated` (empirical variance-shrinkage on paired-difference variance with configurable prior degrees of freedom).

### Robustness protocol

Three robustness checks were run.

1. Leave-one-subject-out (LOO): 10 reruns, each excluding one participant across all conditions before DE recomputation.
2. Ranking stability: overlap and Jaccard index of top-20 proteins between baseline and each LOO rerun.
3. Threshold sensitivity: DE rerun with `min_pairs` set to `4`, `5`, and `6`.

Stability metrics reported were:

1. per-protein direction sign agreement between baseline and LOO runs,
2. retention count of baseline acute `q<0.05` proteins under LOO at `q<0.05` and `q<0.10`,
3. top-20 overlap distribution by comparison.

Mixed-model sensitivity was additionally run for acute proteins with baseline q<0.10 using random-intercept models (`abundance ~ condition + (1|subject)`), and compared to paired t-test effects by sign concordance and rank concordance [8,10].

### Heterogeneity metrics

Per-subject per-protein deltas were computed for:

1. `pt1_pr1 = PT1 - PR1`,
2. `pt2_pr2 = PT2 - PR2`,
3. `pt2_pt1 = PT2 - PT1`,
4. `pr2_pr1 = PR2 - PR1`.

Subject-level summaries:

1. acute burden = mean absolute value of (`pt1_pr1`, `pt2_pr2`) across proteins;
2. adaptation burden = mean absolute value of (`pt2_pt1`, `pr2_pr1`) across proteins.

Responder classes were assigned by median split on both axes:

1. HighAcute_HighAdapt,
2. HighAcute_LowAdapt,
3. LowAcute_HighAdapt,
4. LowAcute_LowAdapt.

This median-split scheme was used as a descriptive visualization heuristic in this study, not as a confirmatory classifier.

Additional heterogeneity readouts:

1. recurrent subject fingerprint proteins (frequency in per-subject top-5 absolute acute deltas),
2. inter-individual acute variance ranking (protein SD across subjects),
3. module trajectory scores for heat-shock, immune, vascular, and ECM-remodeling feature sets.

Module trajectory scores were computed as mean per-subject delta across pre-specified small protein sets represented on the Olink panel (heat-shock n=4, immune n=10, vascular n=7, ECM-remodeling n=8), with one score per module per contrast. This follows a gene-set-level summarization logic commonly used to stabilize interpretation in high-dimensional expression/proteomics analyses [14]. Module curation was knowledge-driven from canonical heat-response biology (stress-chaperone, immune signaling/adhesion, vascular regulation, and matrix-remodeling proteins present on the panel) rather than data-driven optimization. The heat-shock module used HSPA1A, DNAJB1, HSPB1, and HSPG2, where HSPG2 was included as a heat-responsive ECM proteoglycan feature rather than a classical intracellular chaperone. Full per-module protein membership should be provided in the final supplementary deposition alongside the current score tables.

Burden sensitivity analyses included:

1. top-500 acute-variance restricted burden,
2. inverse-variance weighted burden.

Agreement with baseline burden was quantified with Spearman correlation.

### Latent-axis construction

A feature matrix was constructed from burden variables plus module trajectory contrasts. Features were z-scored and decomposed by singular value decomposition (SVD), and the first two components were retained as latent axes (PC1, PC2) [9]. Subject archetype labels combined latent-axis median splits with burden-class labels.

### Physiology mapping

Physiological endpoints were extracted from `Physiological_data.xlsx` as day-1 heat deltas, day-7 heat deltas, and acclimation deltas. Pearson correlations were computed between each latent/burden metric (PC1, PC2, acute burden, adaptation burden) and each physiological endpoint using complete pairs (`n>=5`). BH correction was applied within each metric family of phenotype tests.

To test whether heterogeneity magnitudes exceeded a label-randomized null, subject labels were permuted within each protein (500 iterations) and burden dispersion statistics were recomputed. Empirical p-values were computed as tail probabilities of permutation statistics relative to observed values, and interpreted with finite-permutation resolution constraints [13].

### Figure generation

Figures were generated from the final robustness and heterogeneity tables and exported to `docs/findings/heterogeneity/figure_pack`. Panel captions are documented in `Figure_Captions.md`.

## Data availability

All processed analysis artifacts, robustness tables, heterogeneity tables, and figure panels will be deposited to a public Zenodo archive at submission. Accession/DOI: pending deposition.

## Code availability

karnaProteome source code and exact analysis scripts/commands used in this study will be archived at submission. Repository DOI/release tag: pending deposition.

## Figure availability

All figure panels and legends are included in the manuscript figure package and will be publicly available with the deposition record at submission.

## Figure legends

### Figure 1

**Fig1A. Subject-level proteomic response map.** Scatter of per-subject acute burden versus adaptation burden. Acute burden is mean absolute protein delta across PT1-PR1 and PT2-PR2; adaptation burden is mean absolute protein delta across PT2-PT1 and PR2-PR1. Points are colored by responder class and labeled by subject ID.

**Fig1B. Responder class composition.** Bar chart of subject counts per responder class (HighAcute_HighAdapt, HighAcute_LowAdapt, LowAcute_HighAdapt, LowAcute_LowAdapt) derived from median splits of acute and adaptation burdens.

### Figure 2

**Fig2A. Individual burden profiles.** Per-subject paired bars for acute burden and adaptation burden with connector lines to emphasize within-subject profile shape.

**Fig2B. Proteins with highest inter-individual acute variability.** Heatmap of subject-level acute mean deltas (mean of PT1-PR1 and PT2-PR2) for top variable proteins ranked by inter-individual standard deviation.

### Figure 3

**Fig3A. Volcano plot for PT1-PR1.** Differential abundance volcano for PT1-PR1 using mean difference on the x-axis and -log10(q) on the y-axis. Dashed line marks q=0.05; highlighted points are proteins with q<0.05.

**Fig3B. Volcano plot for PT2-PR2.** Differential abundance volcano for PT2-PR2 with the same encoding as Fig3A for direct acute-contrast comparison.

**Fig3C. Discovery counts across comparisons.** Grouped bars showing number of proteins at q<0.05 and q<0.10 for each comparison (PT1-PR1, PT2-PR2, PT2-PT1, PR2-PR1).

### Figure 4

**Fig4A. LOO retention of baseline acute q<0.05 proteins.** Per-protein retention counts under leave-one-subject-out reruns for proteins significant at baseline (q<0.05) in acute contrasts. Y-axis reports number of LOO reruns (0-10) preserving q<0.05.

**Fig4B. Top-rank stability under leave-one-subject-out.** Top-20 overlap counts by left-out subject for PT1-PR1 and PT2-PR2, quantifying ranking stability of leading acute signals.

### Figure 5

**Fig5. Burden distribution by responder class.** Class-stratified strip plot of acute and adaptation burdens, showing separation and spread within each responder class.

### Figure 6

**Fig6A. Recurrent subject-fingerprint proteins.** Frequency plot of proteins appearing in per-subject top-5 absolute deltas (separately for PT1-PR1 and PT2-PR2), highlighting recurring versus private response features.

**Fig6B. Distribution of top variable proteins across subjects.** Horizontal boxplots of acute mean deltas for highest inter-individual-variance proteins, summarizing spread, median shift, and outlier structure.

### Figure 7

**Fig7A. Subject archetype map on latent axes.** Two-dimensional latent projection (PC1, PC2) of subjects derived from burden and module-trajectory features. Points are colored by composite archetype labels with median axis boundaries.

**Fig7B. Module trajectory score heatmap.** Subject-by-feature heatmap of module trajectory scores across contrasts (PT1-PR1, PT2-PR2, PT2-PT1, PR2-PR1) for heat-shock, immune, vascular, and ECM-remodeling modules.

### Figure 8

**Fig8A. Dominant PC1 feature loadings.** Top absolute loadings contributing to latent axis PC1, indicating which burden/module features most strongly separate subjects along the primary heterogeneity axis.

**Fig8B. Physiology mapping heatmap.** Pearson correlation heatmap between latent metrics (PC1, PC2, acute burden, adaptation burden) and physiological delta endpoints from the Dube phenotype sheet. Values are correlation coefficients.

## References

1. Périard JD, Racinais S, Sawka MN. Adaptations and mechanisms of human heat acclimation: applications for competitive athletes and sports. Scand J Med Sci Sports. 2015;25 Suppl 1:20-38. doi:10.1111/sms.12408.
2. Gagnon D, Barry H, Barhdadi A, Oussaid E, Mongrain I, Lemieux Perreault LP, et al. A dataset of proteomic changes during human heat stress and heat acclimation. Sci Data. 2023;10:877. doi:10.1038/s41597-023-02809-5.
3. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc B. 1995;57(1):289-300.
4. Shi L, Reid LH, Jones WD, Shippy R, Warrington JA, Baker SC, et al. The MicroArray Quality Control (MAQC) project shows inter- and intraplatform reproducibility of gene expression measurements. Nat Biotechnol. 2006;24(9):1151-61. doi:10.1038/nbt1239.
5. Sawka MN, Leon LR, Montain SJ, Sonna LA. Integrated physiological mechanisms of exercise performance, adaptation, and maladaptation to heat stress. Compr Physiol. 2011;1(4):1883-928. doi:10.1002/cphy.c100082.
6. Tyler CJ, Reeve T, Hodges GJ, Cheung SS. The effects of heat adaptation on physiology, perception and exercise performance in the heat: a meta-analysis. Sports Med. 2016;46(11):1699-724. doi:10.1007/s40279-016-0538-5.
7. Assarsson E, Lundberg M, Holmquist G, Bjorkesten J, Bucht Thorsen S, Ekman D, et al. Homogenous 96-plex PEA immunoassay exhibiting high sensitivity, specificity, and excellent scalability. PLoS One. 2014;9(4):e95192. doi:10.1371/journal.pone.0095192.
8. Bates D, Machler M, Bolker B, Walker S. Fitting linear mixed-effects models using lme4. J Stat Softw. 2015;67(1):1-48. doi:10.18637/jss.v067.i01.
9. Jolliffe IT, Cadima J. Principal component analysis: a review and recent developments. Philos Trans A Math Phys Eng Sci. 2016;374(2065):20150202. doi:10.1098/rsta.2015.0202.
10. Laird NM, Ware JH. Random-effects models for longitudinal data. Biometrics. 1982;38(4):963-74. doi:10.2307/2529876.
11. Storey JD, Tibshirani R. Statistical significance for genomewide studies. Proc Natl Acad Sci U S A. 2003;100(16):9440-5. doi:10.1073/pnas.1530509100.
12. Leek JT, Scharpf RB, Bravo HC, Simcha D, Langmead B, Johnson WE, et al. Tackling the widespread and critical impact of batch effects in high-throughput data. Nat Rev Genet. 2010;11(10):733-9. doi:10.1038/nrg2825.
13. Phipson B, Smyth GK. Permutation p-values should never be zero: calculating exact p-values when permutations are randomly drawn. Stat Appl Genet Mol Biol. 2010;9(1):Article39. doi:10.2202/1544-6115.1585.
14. Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A. 2005;102(43):15545-50. doi:10.1073/pnas.0506580102.
15. Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43(7):e47. doi:10.1093/nar/gkv007.
