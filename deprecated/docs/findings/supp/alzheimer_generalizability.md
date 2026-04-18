# Supplementary — unpaired generalizability on the Wei et al. 2025 Alzheimer's plasma Olink Target 96 cohort

## Dataset

Wei J, Chen G, Ni H. *NPX values from Olink proteomic analysis of Alzheimer's disease plasma samples*. Figshare dataset, 2025. [10.6084/m9.figshare.28829720.v1](https://doi.org/10.6084/m9.figshare.28829720.v1).

- Platform: **Olink Target 96 Inflammation (v.3024)**
- Design: **Unpaired case-control, three groups**
  - A = Alzheimer's disease (AD), n=10
  - B = Mild cognitive impairment (MCI), n=10
  - C = Healthy control (HC), n=10
- Matrix: plasma
- Assays: 92 inflammation proteins
- QC: 4 samples carry "Warning" flags in the source export (groupA9, groupB9, groupC4, groupC9); these are dropped at the `effective_abundance()` layer of atman and excluded from DE, leaving n=9 AD, 9 MCI, 8 HC for analysis.
- Data source: Figshare (mirrored locally at `example_data/alzheimer_olink_2024/NPX_values_Alzheimers_study_Samples.xlsx`).

This is the **first unpaired** dataset exercised through atman and validates the newly added `--test welch-t` mode.

## Ingest adapter

`scripts/ingest_alzheimer.py` (~90 lines of Python) parses the Olink NPX Manager xlsx export ("Npx Data" sheet) into atman's internal long-format TSVs. Mapping:

| xlsx row / col                          | atman TSV field                         |
|-----------------------------------------|-----------------------------------------|
| row 3 "Assay", cols 1..92               | `gene_symbol` in `proteins.tsv`         |
| row 4 "Uniprot ID", cols 1..92          | `uniprot`                               |
| row 5 "OlinkID", cols 1..92             | `assay_id`                              |
| rows 7..36, col 0 ("groupA1" etc.)      | `sample_id` in `samples.tsv`            |
| derived from group letter (A/B/C)       | `condition` ∈ {AD, MCI, HC}             |
| rows 7..36, col 93 "Plate ID"           | `plate_id` in `measurements.tsv`        |
| rows 7..36, col 94 "QC Warning"         | `qc_sample` ("Pass" → PASS, else WARN)  |

No changes to atman's Rust engine are required; the adapter produces the same file contract used by the Dube and Gisby pipelines.

## Welch's two-sample t-test in atman

Added as `atman_core::de::welch_t` (~60 lines Rust, 5 unit tests). Used when the CLI flag `--test welch-t` is supplied; `--paired-by` is ignored in this mode. Two groups $A$ and $B$ with $n_A$ and $n_B$ independent samples:

$$
t = \frac{\bar x_A - \bar x_B}{\sqrt{s_A^2/n_A + s_B^2/n_B}}, \qquad
\nu = \frac{(s_A^2/n_A + s_B^2/n_B)^2}{(s_A^2/n_A)^2/(n_A-1) + (s_B^2/n_B)^2/(n_B-1)}.
$$

Two-sided $p$-values under Student's $t$-distribution with $\nu$ degrees of freedom, BH-FDR within each contrast family. Minimum samples per group enforced by `--min-pairs` (overloaded for unpaired mode to mean `min(|A|, |B|)`).

### Numerical verification

Atman's Welch-t output was cross-validated against `scipy.stats.ttest_ind(..., equal_var=False)` on all 276 protein × contrast cells of this dataset (AD-HC, MCI-HC, AD-MCI × 92 proteins):

| metric                           | max |Δ| |
|----------------------------------|----------|
| `|t_atman − t_scipy|`            | 3.9 × 10⁻¹⁴ |
| `|df_atman − df_scipy|`          | 8.9 × 10⁻¹⁵ |
| `|p_atman − p_scipy|`            | 8.4 × 10⁻¹⁴ |

i.e., agreement to IEEE 754 machine precision.

## Pipeline run

```bash
python3 scripts/ingest_alzheimer.py example_data/alzheimer_olink_2024 out_ad
atman de --input-dir out_ad --output-dir out_ad \
    --test welch-t --groups "AD-HC,MCI-HC,AD-MCI" --min-pairs 5
```

or one-shot: `bash scripts/reproduce_alzheimer.sh`.

## Differential abundance — Welch-t, three contrasts

Per-contrast hit counts at `q < 0.05` and `q < 0.10`:

| contrast | n_tested | q<0.05 | q<0.10 |
|----------|----------|--------|--------|
| AD-HC    | 92       | **10** | 13     |
| MCI-HC   | 92       | **10** | 10     |
| AD-MCI   | 92       |  0     |  0     |

Zero hits in AD vs MCI is biologically expected: MCI is an intermediate state on the AD trajectory and plasma inflammation profiles tend to be similar at this cohort size.

### AD-HC — top hits (q < 0.05)

| gene | mean_diff (AD − HC) | t | df | q-value |
|------|--------|------|------|---------|
| IL18     | −3.61 | −12.66 | 11.1 | 5.8×10⁻⁶ |
| CASP-8   | −3.55 | −6.84  | 9.1  | 2.3×10⁻³ |
| EN-RAGE (S100A12) | −3.55 | −11.62 | 10.2 | 1.5×10⁻⁵ |
| IL8      | −3.13 | −4.55  | 7.9  | 2.0×10⁻² |
| ST1A1 (SULT1A1) | −1.52 | −4.27 | 12.8 | 1.7×10⁻² |
| 4E-BP1   | −1.28 | −4.15  | 13.3 | 1.7×10⁻² |
| TNFSF14 (LIGHT) | −1.24 | −3.88 | 14.9 | 1.7×10⁻² |
| OSM      | −1.15 | −4.08  | 12.9 | 1.7×10⁻² |

Direction: the source dataset reports higher NPX values in the HC group for all hits. Whether this reflects cohort-level selection effects (HC group composition, pre-analytic differences) or genuine inverse regulation is a question for the dataset authors and follow-up work; atman records the direction the data support and does not correct for it.

### MCI-HC — top hits (q < 0.05)

Largely overlapping with AD-HC: CASP-8, IL18, EN-RAGE, IL8, OSM, TNFSF14, ST1A1, TGF-alpha. Consistent with MCI being a prodromal state sharing a subset of the inflammatory profile differences.

Full per-protein TSVs: `out_ad/de_results.tsv`, `out_ad/de_report.tsv`.

## Interpretive discipline at n=10 per group

Welch-t is the standard default for case-control designs when group variances may differ. The small cohort (8–10 per group after QC) limits per-protein power, but the Welch denominator correctly absorbs the variance asymmetry visible in these data (IL18: σ_HC=0.82 vs σ_AD=0.42). The three-group design allows a dose-response-style sanity check: proteins that differ in AD-HC should also differ in MCI-HC if the biology is a gradient rather than an AD-specific step change, and that is what we observe. The null AD-MCI contrast is consistent with this reading.

No mechanistic claims are made from this reanalysis. The demonstration's purpose is methodological: atman's Welch-t implementation produces scipy-identical output on real data, runs end-to-end on a case-control Olink dataset, and emits the same TSV-in / TSV-out contract as the paired pipelines.

## Deposited artefacts

```
example_data/alzheimer_olink_2024/   ← raw xlsx (mirrored from figshare)
scripts/ingest_alzheimer.py          ← ~90-line adapter
scripts/reproduce_alzheimer.sh       ← one-command end-to-end rerun
out_ad/                              ← pipeline outputs
  measurements.tsv, qc_measurements.tsv, proteins.tsv, samples.tsv
  Inflammation_log2_fc.csv
  de_results.tsv, de_report.tsv
```
