# Atman Tutorial

This walkthrough runs Atman's Dube et al. 2023 Olink Explore NGS reproduction
path against the fixture data shipped under `example_data/dube_heat_2023/`.

Atman's built-in direct ingest currently targets Olink Explore long CSV files.
Other proteomics formats can be used by writing Atman's canonical TSV files
(`samples.tsv`, `proteins.tsv`, `measurements.tsv`, and optionally
`qc_measurements.tsv`) and then running the downstream analysis commands on
that directory.

## 0. Install and verify

```bash
git clone https://github.com/kevinjoseph/atman.git
cd atman
cargo install --path crates/atman
atman --help
```

Expected: help text listing `ingest`, `validate`, `qc`, `report`, `matrix`,
`fold-change`, `de`, `asymmetry`, `robustness`, `module-trajectory`, and
`module-de`.

Run the test suite to confirm the reproduction base is intact:

```bash
cargo test --workspace --release
```

Expected: all test binaries report `ok` with zero failures.

## 1. Ingest

```bash
atman ingest \
    --platform olink-explore-ngs --parser dube \
    --output-dir out \
    example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv \
    example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv
```

Expected artefacts in `out/`:
- `measurements.tsv` — long-format rows (≈122,135)
- `proteins.tsv` — protein/assay catalogue
- `samples.tsv` — 43 samples (40 biological + 3 technical controls)

Sanity check:
```bash
wc -l out/measurements.tsv   # ≈122,136 including header
```

## 2. QC (Dube rule)

```bash
atman qc --input-dir out --output-dir out --rule dube
```

Masks rows where `QC_Warning` or `Assay_Warning` is not `PASS`. Raw values are
preserved; only the inferential abundance is nulled.

Expected: 3,725 rows masked, 118,410 usable rows remain.

For non-Olink sources, an adapter can write `qc_measurements.tsv` directly if
the source matrix is already filtered or if QC has been handled upstream.

## 3. Validate

```bash
atman validate --input-dir out \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" \
    --min-pairs 5 \
    --report out/validate_report.tsv
```

Checks the canonical TSV files for required columns, duplicate keys, orphan
sample or assay references, parseable platforms, parseable QC/boolean fields,
known abundance units, and enough effective samples for requested contrasts.

Expected: `0 error(s), 0 warning(s)` and a report file containing only the
header row.

For adapter-generated canonical TSVs, run `validate` before downstream
analysis. Unknown abundance units are warnings by default; add `--strict` to
fail CI on warnings.

## 4. QC report

```bash
atman report qc --input-dir out --output-dir out/report
```

Outputs:
- `out/report/qc_summary.tsv` — dataset-level counts and fractions
- `out/report/sample_qc.tsv` — per-sample missingness and QC masking
- `out/report/protein_qc.tsv` — per-protein missingness and support
- `out/report/condition_counts.tsv` — samples, subjects, and effective
  measurements by condition

The command emits warnings for small effective condition groups and sparse
proteins. Tune those thresholds with `--min-subjects` and `--sparse-threshold`.

## 5. Matrix (Dube-wide pivot)

```bash
atman matrix --input-dir out --output-dir out --format dube-wide --split-by panel
```

Emits per-panel wide NPX CSVs matching Dube's published filtered NPX files
cell-exactly (8 panels).

Verify against the published reference:
```bash
diff out/cardiometabolic_npx.csv \
     example_data/dube_heat_2023/filtered_npx/npx/cardiometabolic_npx.csv
```
Expected: empty output (byte-exact).

## 6. Fold change

```bash
atman fold-change --input-dir out --output-dir out \
    --groups "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1"
```

Emits 8 per-panel log2 fold-change CSVs. Maximum delta against Dube's
published FC tables is 1.05 × 10⁻¹⁵ (IEEE 754 last-bit drift).

## 7. Differential abundance

Paired Student's t-test at subject level:

```bash
atman de --input-dir out --output-dir out \
    --test paired-t --paired-by participant \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" --min-pairs 5
```

Output: `out/de_results.tsv` with columns
`panel, assay_id, gene_symbol, uniprot, comparison, n_pairs, mean_a, mean_b, mean_diff, t, df, p_value, bh_q, skip_reason`.

Expected hit counts at `q < 0.05`:
| contrast | proteins at q<0.05 |
|---|---|
| PT1-PR1 | 36 |
| PT2-PR2 | 37 |
| PT2-PT1 | 7 |
| PR2-PR1 | 2 |

Optional moderated variance-shrinkage alternative:
```bash
atman de --input-dir out --output-dir out_mod \
    --test moderated --moderation-prior-df 4 --paired-by participant \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" --min-pairs 5
```

## 8. Asymmetry

```bash
atman asymmetry \
    --de-results out/de_results.tsv \
    --pairs "PT2-PT1:PR2-PR1,PT2-PR2:PT1-PR1" \
    --output out/asymmetry.tsv
```

Emits per-pair asymmetry metrics (challenge-vs-rest) at the contrast-family level.

## 9. Robustness

Robustness consumes the baseline DE table plus one or more rerun DE tables.
The example below shows the command shape using placeholder LOO outputs:

```bash
atman robustness \
    --baseline out/de_results.tsv \
    --loo "out_loo_001/de_results.tsv,out_loo_002/de_results.tsv" \
    --top-k 20 \
    --output-dir out/robustness
```

Output in `out/robustness/`:
- `loo_gene_stability.tsv` — per-protein sign/q retention
- `rank_stability.tsv` — top-20 overlap and Jaccard by excluded subject
- `threshold_sensitivity.tsv` — `min_pairs` sweep

For a compact end-to-end check, run `cargo test --workspace --release`; the
integration tests execute this ingest/QC/matrix/fold-change path and compare
against the bundled reference tables.
