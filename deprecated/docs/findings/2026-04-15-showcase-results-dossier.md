# Atman Showcase Results Dossier (2026-04-15)

## What This File Is

This is the main index for all showcase results from the Dube heat-stress dataset.

It links to:

1. the core pipeline run,
2. robustness checks,
3. human-variation results,
4. latent-axis and physiology mapping,
5. the full figure pack.

## Plain-English Summary

1. We ran the full analysis from raw proteomics files to final statistics, and it ran cleanly.
2. The strongest shared signal is acute heat response (`PT1-PR1` and `PT2-PR2`).
3. Known heat-stress proteins stay directionally consistent even when leaving one subject out.
4. People in this cohort do not respond the same way. Response size and pattern vary clearly by person.
5. We can group subjects into responder types (high/low acute response, high/low adaptation response).
6. Some proteins repeatedly show up as high-variation markers, but these need external validation.

## What We Can Claim Safely

### Strong claims

1. There is a real cohort-level acute heat proteomic signal.
2. That acute directionality is robust in sensitivity checks.
3. Inter-individual heterogeneity is a core biological feature in this cohort.

### Medium-confidence claims

1. A small set of proteins may help separate responder types.
2. Latent subject axes may reflect adaptation biology.

### Early/exploratory claims

1. Specific non-canonical single-protein `q<0.05` hits.
2. Specific physiology links after multiple-testing correction in `n=10`.

## Results Documents

1. Week 1 execution:  
   [2026-04-14-showcase-week1-execution.md](/Users/kevinjoseph/Cursor/atman/docs/findings/2026-04-14-showcase-week1-execution.md)
2. Robustness summary:  
   [2026-04-14-robustness-summary.md](/Users/kevinjoseph/Cursor/atman/docs/findings/robustness/2026-04-14-robustness-summary.md)
3. Heterogeneity summary:  
   [2026-04-14-heterogeneity-summary.md](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/2026-04-14-heterogeneity-summary.md)

## Core Data Artifacts

### Robustness

1. [loo_gene_stability.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/robustness/loo_gene_stability.tsv)
2. [rank_stability.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/robustness/rank_stability.tsv)
3. [threshold_sensitivity.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/robustness/threshold_sensitivity.tsv)

### Heterogeneity

1. [per_subject_gene_deltas.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/per_subject_gene_deltas.tsv)
2. [subject_response_burden.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/subject_response_burden.tsv)
3. [responder_clusters.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/responder_clusters.tsv)
4. [subject_top_proteins.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/subject_top_proteins.tsv)
5. [high_variance_acute_genes.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/high_variance_acute_genes.tsv)
6. [module_trajectory_scores.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/module_trajectory_scores.tsv)
7. [latent_axes.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/latent_axes.tsv)
8. [subject_archetypes.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/subject_archetypes.tsv)
9. [physiology_mapping.tsv](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/physiology_mapping.tsv)

## Manuscript Figure Pack

Directory:  
[figure_pack](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/figure_pack)

Panels:

1. `Fig1A_burden_scatter.png`
2. `Fig1B_responder_class_counts.png`
3. `Fig2A_subject_burden_bars.png`
4. `Fig2B_top_variable_heatmap.png`
5. `Fig3A_volcano_PT1_PR1.png`
6. `Fig3B_volcano_PT2_PR2.png`
7. `Fig3C_discovery_counts.png`
8. `Fig4A_loo_q05_retention.png`
9. `Fig4B_rank_overlap_by_subject.png`
10. `Fig5A_burden_by_class.png`
11. `Fig5B_pathway_counts_heatmap.png`
12. `Fig6A_recurrent_fingerprint_proteins.png`
13. `Fig6B_top_variable_boxplots.png`
14. `Fig7A_subject_archetype_map.png`
15. `Fig7B_module_trajectory_heatmap.png`
16. `Fig8A_PC1_loadings.png`
17. `Fig8B_physiology_mapping_heatmap.png`

Captions:

1. [Figure_Captions.md](/Users/kevinjoseph/Cursor/atman/docs/findings/heterogeneity/figure_pack/Figure_Captions.md)

## Notes for Manuscript Writing

1. Lead with reproducibility and robust acute signal.
2. Treat heterogeneity as the main human-biology result, not a side note.
3. Be conservative with strict single-protein significance at this sample size.
4. Present latent-axis and physiology links as informative structure, not final causality.
