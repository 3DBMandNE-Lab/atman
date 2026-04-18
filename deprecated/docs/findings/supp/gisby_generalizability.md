# Supplementary — generalizability on Gisby 2021 (Olink Target 96, COVID dialysis)

All numbers in this document are computed directly from primary TSVs in
`out_gisby/` and `out_gisby_loo/`. The exact source path for each value is
given in brackets next to the value.

## Cohort

- Source dataset: Gisby *et al.* 2021, *eLife* 10:e64827, longitudinal
  Olink Target 96 profiling of ESKD dialysis patients with COVID-19.
- Analysis cohort: 23 patients with both an early (≤7 days from first
  swab) and a late (≥14 days) sample. [`out_gisby/samples.tsv` — 46 rows
  of (patient × timepoint)]
- Adapter: `scripts/ingest_gisby.py` (~100 lines of Python). No Atman
  engine changes were needed — the adapter writes the standard long-format
  `measurements.tsv` + `samples.tsv` + `proteins.tsv` contract.
- Contrast: within-subject `late − early`.

## Differential abundance: counts

Counts computed from `out_gisby/de_results.tsv`, filter
`comparison="late-early"` and `skip_reason` empty.

| quantity                           | value           |
|------------------------------------|-----------------|
| proteins tested                    | 436             |
| hits at q<0.05                     | 120             |
| hits at q<0.10                     | 157             |

## Leading hits

Computed from `out_gisby/de_results.tsv`.

**Top 10 down-regulated hits (late < early)**, sorted by `mean_diff`
ascending among rows with q<0.05:

| gene_symbol | mean_diff | bh_q                |
|-------------|-----------|---------------------|
| IFN-gamma   | −2.541    | 8.74 × 10⁻⁴         |
| CXCL10      | −1.880    | 5.48 × 10⁻⁶         |
| DDX58       | −1.711    | 4.43 × 10⁻⁵         |
| MCP-2       | −1.508    | 1.99 × 10⁻⁵         |
| KRT19       | −1.036    | 2.22 × 10⁻³         |
| TRIM21      | −1.031    | 1.51 × 10⁻⁴         |
| MCP-3       | −0.928    | 1.63 × 10⁻²         |
| IL10        | −0.890    | 3.24 × 10⁻²         |
| IGFBP-1     | −0.887    | 6.73 × 10⁻³         |
| TNFSF13B    | −0.885    | 5.29 × 10⁻⁴         |

Leading down-regulation is the canonical interferon-response /
acute-inflammation axis, consistent with resolution of acute viral
inflammation between the early and late timepoints.

**Top 10 up-regulated hits (late > early)**, sorted by `mean_diff`
descending among rows with q<0.05:

| gene_symbol | mean_diff | bh_q                |
|-------------|-----------|---------------------|
| CCL24       | +1.596    | 2.24 × 10⁻⁵         |
| CCL16       | +1.262    | 1.65 × 10⁻⁴         |
| PSP-D       | +0.820    | 9.45 × 10⁻³         |
| CCL17       | +0.728    | 2.31 × 10⁻²         |
| SERPINA12   | +0.688    | 6.19 × 10⁻⁴         |
| SERPINA5    | +0.667    | 5.21 × 10⁻³         |
| MCP-4       | +0.617    | 2.02 × 10⁻²         |
| TNFRSF9     | +0.583    | 3.97 × 10⁻⁴         |
| TNFRSF10C   | +0.569    | 5.48 × 10⁻⁶         |
| CCL18       | +0.555    | 1.48 × 10⁻²         |

Leading up-regulation: convalescent-phase chemokines (CCL24, CCL16, CCL17,
CCL18) and matrix-regulation factors (SERPINA family, TNFRSF co-stimulators).

## LOO robustness across 23 patients

Computed from `out_gisby/robustness/loo_sign_stability.tsv` and
`out_gisby/robustness/rank_stability.tsv`.

| quantity                                                    | value  |
|-------------------------------------------------------------|--------|
| mean per-protein sign-match rate (across all 436 assays)    | 0.9718 |
| mean top-20 rank overlap (/20 across 23 LOO runs)           | 19.04  |
| mean top-20 Jaccard (across 23 LOO runs)                    | 0.9112 |
| baseline q<0.05 hits retained at q<0.05 in *all* 23 LOO     | 98 / 120 |
| baseline q<0.05 hits retained at q<0.10 in *all* 23 LOO     | 120 / 120 |

## Comparison to Dube (n=10) robustness geometry

All Dube numbers computed from `out/robustness/loo_sign_stability.tsv` and
`out/robustness/rank_stability.tsv` by the same rules.

| metric                                                   | Dube PT1-PR1 (n=10) | Dube PT2-PR2 (n=10) | Gisby late-early (n=23) |
|----------------------------------------------------------|---------------------|---------------------|-------------------------|
| mean sign-match rate                                     | 0.9358              | 0.9398              | **0.9718**              |
| mean top-20 overlap (/20)                                | 16.70               | 18.30               | **19.04**               |
| mean top-20 Jaccard                                      | 0.7238              | 0.8473              | **0.9112**              |
| baseline q<0.05 retained q<0.05 in all LOO reruns        | 11 / 36 (31%)       | 5 / 37 (14%)        | **98 / 120 (82%)**      |
| baseline q<0.05 retained q<0.10 in all LOO reruns        | 21 / 36 (58%)       | 26 / 37 (70%)       | **120 / 120 (100%)**    |

Stability margins tighten across every metric at n=23 vs n=10, as
expected. The improvement is visible in the standard output schema — the
stability readouts themselves become informative about sample-size adequacy.

## Data-driven multiscale inference

Computed from `out_gisby/module_de_results.tsv` and
`out_gisby/modules_learned.tsv`.

- Modules learned: 15, by complete-linkage clustering of per-subject
  delta correlations on the filtered assay set (cut height tuned to K=15).
  [`scripts/learn_modules.py`; modules written to
  `out_gisby/modules_learned.tsv`, 436 rows for 15 modules]
- Module-DE hits: 6 / 15 at q<0.05; 6 / 15 at q<0.10.
- Top module hit: `module_14` (14 genes), `mean_diff = −0.462`,
  `bh_q = 3.46 × 10⁻⁴`. [`out_gisby/module_de_results.tsv` row
  `(module_14, late-early)`]

`module_14` membership
[`out_gisby/modules_learned.tsv` rows where `module="module_14"`; 14 rows]:

```
CCL19, CCL23, CXCL10, GRN, IGFBP-1, IL-20RA, IL-2RB, IL18, IL33, LAG3,
MCP-3, RARRES2, TIMP1, TNC
```

This is the canonical type-II interferon / chemokine axis. The module
emerges unsupervised from correlation structure of per-patient late−early
deltas; no gene lists or pathway databases were used.

## Module-scale uncertainty via subject bootstrap

Computed from `out_gisby/robustness/module_bootstrap.tsv` (B=1,000 subject
resamples, contrast `late-early`).

| module    | n_genes | boot_mean | CI₉₅ lower | CI₉₅ upper | P(sign stable) |
|-----------|---------|-----------|------------|------------|----------------|
| module_03 | 40      | −0.165    | −0.266     | −0.070     | 1.000          |
| module_04 | 37      | −0.273    | −0.393     | −0.158     | 1.000          |
| module_09 | 24      | +0.252    | +0.107     | +0.413     | 1.000          |
| module_14 | 14      | −0.460    | −0.642     | −0.299     | 1.000          |

Four modules reach P(sign stable)=1.000. Two more reach ≥0.95 (see
`module_uncertainty.md` for the full 15-module table).

## Protein-scale uncertainty via subject bootstrap

Computed from `out_gisby/robustness/protein_bootstrap.tsv` (B=1,000,
contrast `late-early`, 435 proteins).

| threshold               | count |
|-------------------------|-------|
| total proteins          | 435   |
| P(sign stable) ≥ 0.95   | 220   |
| P(sign stable) ≥ 0.999  | 114   |
| P(sign stable) = 1.000  | 101   |

Canonical interferon-response hits match the classical q-ranking with
machine-precision directional certainty:

| gene         | boot_mean | CI₉₅               | P(sign stable) |
|--------------|-----------|--------------------|----------------|
| IFN-gamma    | −2.51     | [−3.50, −1.52]     | 1.000          |
| CXCL10       | −1.87     | [−2.29, −1.47]     | 1.000          |
| TNFSF13B     | −0.89     | [−1.20, −0.59]     | 1.000          |
| IRF9         | −0.84     | [−1.12, −0.56]     | 1.000          |
| GRN          | −0.61     | [−0.76, −0.46]     | 1.000          |

## Deposited artefacts

```
example_data/gisby_covid_2021/   # raw CSVs mirrored from jackgisby/longitudinal_olink_proteomics
scripts/ingest_gisby.py          # adapter, ~100 lines
out_gisby/                       # pipeline outputs
  measurements.tsv, qc_measurements.tsv, proteins.tsv, samples.tsv
  <panel>_log2_fc.csv × 5
  de_results.tsv, de_report.tsv, module_de_results.tsv, modules_learned.tsv
  robustness/loo_sign_stability.tsv, rank_stability.tsv
  robustness/module_bootstrap.tsv, protein_bootstrap.tsv
  robustness/module_de_null_late-early_summary.tsv
out_gisby_loo/                   # 23 LOO reruns (one subdirectory per dropped patient)
```

## Reproducibility

```
bash scripts/reproduce_gisby.sh   # end-to-end rerun of this demonstration
```
