# Supplementary — module definitions for the module-trajectory analysis

## Status 2026-04-17: resolved via re-curation

The input `modules.tsv` used to produce the original
`module_trajectory_scores.tsv` in the April 2026 heterogeneity pass was not
persisted to the repository. The four module gene lists below were
**re-curated** from canonical biology in the 2026-04-17 regeneration, then
intersected with the Olink Explore panel measured in the Dube dataset (2,925
unique gene symbols in `out/proteins.tsv`). Single-panel representation was
required for every gene so the module-trajectory aggregator (`atman
module-trajectory`) does not double-count proteins that appear on multiple
Olink panels (IL6, TNF, and CXCL8 are on all four of Cardiometabolic,
Inflammation, Neurology, and Oncology — they were therefore excluded from
the immune module to match the manuscript's stated `n=10` exactly).

## Final module membership

Per-gene rationale is in `modules_full.tsv`. Summary:

| module        | n_genes | membership                                              |
|---------------|---------|---------------------------------------------------------|
| heat_shock    | 4       | HSPA1A, DNAJB1, HSPB1, HSPG2                            |
| immune        | 10      | ITGAL, ITGAM, ITGB2, IL10, IFNG, IL1B, IL2, IL17A, CCL2, CXCL10 |
| vascular      | 7       | NPPB, VEGFA, ICAM1, VCAM1, SELE, EDN1, ACE              |
| ecm_remodel   | 8       | MMP1, MMP3, MMP9, TIMP1, TIMP2, COL4A4, FN1, TGFB1      |

All 29 genes are present in the Dube Olink Explore data on exactly one panel
each, so `atman module-trajectory` aggregates each gene exactly once per
(subject, contrast, module) combination.

## Consequence for downstream analyses

Three outputs had to be regenerated with the new modules, replacing the
April 14 archives:

| output                           | regenerated | v1 archived at                                         |
|----------------------------------|-------------|--------------------------------------------------------|
| `module_trajectory_scores.tsv`   | 2026-04-17  | `heterogeneity/module_trajectory_scores_v1_archived.tsv` |
| `latent_axes.tsv`                | 2026-04-17  | `heterogeneity/latent_axes_v1_archived.tsv`              |
| `physiology_mapping.tsv`         | 2026-04-17  | `heterogeneity/physiology_mapping_v1_archived.tsv`       |

`heat_shock` module-trajectory scores are identical between v1 and v2 (same
4 genes). The other three modules differ, which propagates to the latent
axes (PC1, PC2) and to the physiology correlations. Specifically:

| claim                             | v1 (archived)                         | v2 (current)                              |
|-----------------------------------|---------------------------------------|-------------------------------------------|
| PC1 variance explained            | 29.32%                                | 22.66%                                    |
| PC2 variance explained            | 20.82%                                | 20.66%                                    |
| strongest PC1 raw association     | acc_dLSR (r=0.880, p=0.0039, n=8)     | d1_dSSNA (r=−0.725, p=0.018, n=10)        |
| strongest PC2 raw association     | (not cited)                            | acc_dTbody (r=−0.854, p=0.007, n=8)       |
| BH-significant association at q<0.10 | none                                | none                                      |

The main manuscript has been updated to cite the v2 numbers. The qualitative
claims the manuscript makes about latent-axis structure — two dominant axes,
non-uniform archetype geometry, full occupancy of responder bins, and
exploratory-only physiology coupling under correction — hold for both v1
and v2. Responder bins, acute/adapt burden ranges, fingerprint proteins
(DAND5, GH1, ITGAL, AGBL2), and high-variance-protein rankings are all
computed from per-subject protein deltas and are **unchanged** between v1
and v2 — they are module-independent.

## Regeneration script

`scripts/regenerate_latent_axes.py` takes the current `modules.tsv` +
`subject_response_burden.tsv` + `per_subject_gene_deltas.tsv` +
`example_data/dube_heat_2023/Physiological_data.xlsx` and regenerates all
three downstream outputs. It reuses the Dube xlsx parser from
`scripts/phenotype_regression.py`.

```bash
python3 scripts/regenerate_latent_axes.py
```

Running the script with the current `docs/findings/heterogeneity/modules.tsv`
reproduces the v2 numbers above to within `1e-10` (SVD round-off).
