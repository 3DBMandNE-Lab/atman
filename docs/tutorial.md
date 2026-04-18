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

Expected: help text listing `ingest`, `ingest-matrix`, `validate`, `qc`,
`report`, `matrix`, `fold-change`, `de`, `bootstrap`, `asymmetry`,
`robustness`, `module-trajectory`, and `module-de`.

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

Matrix ingest alternative for non-Olink data:

```bash
atman ingest-matrix \
    --matrix matrix.tsv \
    --samples sample_metadata.tsv \
    --orientation proteins-rows \
    --platform spectronaut_report \
    --abundance-unit log2_intensity \
    --sample-id-col sample_id \
    --condition-col condition \
    --assay-id-col Protein.Group \
    --gene-col Genes \
    --output-dir out_ms
```

Use `--orientation samples-rows --proteins protein_metadata.tsv` when rows are
samples and protein IDs are matrix columns. Add `--log2-transform` for positive
linear intensities; non-positive values are written as QC-masked measurements.

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

## 8. Protein bootstrap

Subject-level bootstrap uncertainty for protein effects:

```bash
atman bootstrap protein \
    --input-dir out \
    --groups "PT2-PT1" \
    --test paired-t \
    --n 2000 \
    --seed 1 \
    --output out/protein_bootstrap.tsv
```

For unpaired designs, use `--test welch-t`. The command resamples subjects,
not measurement rows. In paired mode, only matched subjects are resampled.

Output columns include point effect, bootstrap mean effect, confidence
interval, sign stability, effective sample counts, and skip reason.

## 9. Module bootstrap

Subject-level bootstrap uncertainty for user-defined module scores:

```bash
atman bootstrap module \
    --input-dir out \
    --modules-tsv modules.tsv \
    --groups "PT2-PT1" \
    --test paired-t \
    --n 2000 \
    --seed 1 \
    --output out/module_bootstrap.tsv
```

The modules file is tab-separated with `module` and `gene_symbol` columns.
Atman scores each module as the mean effective abundance across observed member
genes in each sample, then resamples subjects. The output reports module effect
intervals, sign stability, sample counts, declared gene count, observed gene
count, and skip reason.

For unpaired designs, use `--test welch-t`.

## 10. Null calibration

Permutation and sign-flip calibration estimate how many discoveries the current
design produces under a no-effect null:

```bash
atman null \
    --input-dir out \
    --output-dir out/null \
    --groups "PT2-PT1" \
    --test paired-t \
    --n 1000 \
    --seed 1
```

Use `--test welch-t` for unpaired label permutation. Use `--test paired-t` for
paired sign flips over matched subjects.

Output in `out/null/`:
- `null_summary.tsv` — observed hit counts, null hit counts, and empirical FDR
  by contrast and q threshold
- `empirical_p.tsv` — per-protein observed p/q values plus empirical p-values

## 11. Asymmetry

```bash
atman asymmetry \
    --de-results out/de_results.tsv \
    --pairs "PT2-PT1:PR2-PR1,PT2-PR2:PT1-PR1" \
    --output out/asymmetry.tsv
```

Emits per-pair asymmetry metrics (challenge-vs-rest) at the contrast-family level.

## 12. Robustness

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
