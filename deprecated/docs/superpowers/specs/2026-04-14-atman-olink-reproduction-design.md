# Atman v0.1 — Olink Explore Reproduction Base

**Status:** Draft design, pending review
**Author:** Kevin Joseph (with Claude)
**Date:** 2026-04-14
**Milestone target:** v0.1
**Reference dataset:** Gagnon, Barry, Barhdadi, Oussaid, Mongrain, Lemieux Perreault, Dubé. *A dataset of proteomic changes during human heat stress and heat acclimation.* Scientific Data (2023). DOI [10.1038/s41597-023-02809-5](https://doi.org/10.1038/s41597-023-02809-5). PMID 38062080. PMC10703874. Figshare: [22811915](https://doi.org/10.6084/m9.figshare.22811915.v1) (raw), [22811891](https://figshare.com/articles/dataset/The_filtered_NPX_values/22811891) (filtered), [22811888](https://figshare.com/articles/dataset/Fold_Change/22811888) (fold change). All CC BY 4.0.

---

## 1. Summary

Atman is a new Rust engine that extends the local-first bioinformatics platform (karna, karnadata, 34 Rust CLIs) into the proteomics modality. It is designed to support multiple proteomics platforms (Olink PEA, mass spectrometry, affinity aptamer assays) via a pluggable ingest layer, but v0.1 implements exactly one adapter: **Olink Explore NGS** in the long-CSV NPX format.

v0.1 scope is deliberately narrow: **reproduce the Dube heat-acclimation published outputs byte-for-byte (filtered NPX) and within 1e-4 (fold change) from the raw NPX files.** No statistical testing, no normalization beyond what Olink already applied, no integrity gate (deferred to v0.2), no non-Olink platforms.

Passing the reproduction integration test is the definition of "v0.1 done."

## 2. Non-goals

Explicitly out of scope for v0.1:

- Karnadata extension (`proteome_sample` table, `validate-proteome` command). Ingest reads raw CSVs directly.
- Differential abundance testing, multiple testing correction, censored/Tobit regression, mixed models.
- Bridge normalization across plates/batches.
- Protein set enrichment.
- Platforms other than Olink Explore NGS long-CSV: Olink Target qPCR, SomaScan, MaxQuant, DIA-NN, Spectronaut, Luminex, Simoa. Their `Platform` enum variants exist but are not implemented.
- Claim emission into knowledge.db or any integration with the evidence pipeline (Crucible, Oracle, Atlas).
- Plotting, reporting, HTML outputs.
- Python or R bindings.
- `peek` / `inspect` subcommands.

Each of these is listed in §13 with its target milestone.

## 3. Architectural pattern

Hybrid: a dedicated engine binary (like karna for spatial) over a generic, platform-agnostic core library. The core is **proteomics-level**, not Olink-level — Olink is one adapter in the ingest layer. Downstream tools (tsvkit, karnaplot, genesets) stay as decoupled consumers via TSV/CSV files.

Workspace layout:

```
atman/
├── Cargo.toml                        # workspace
├── crates/
│   ├── atman/                # engine binary (thin CLI + I/O)
│   │   ├── src/
│   │   │   ├── main.rs               # clap dispatch
│   │   │   ├── commands/
│   │   │   │   ├── ingest.rs
│   │   │   │   ├── qc.rs
│   │   │   │   ├── matrix.rs
│   │   │   │   └── fold_change.rs
│   │   │   └── io.rs                 # atomic-rename writers, TSV/CSV helpers
│   │   └── tests/
│   │       └── dube_reproduction.rs  # end-to-end reproduction integration test
│   └── atman-core/                # platform-agnostic library
│       ├── src/
│       │   ├── lib.rs
│       │   ├── types.rs              # Sample, AssayId, ProteinIdentity,
│       │   │                         #  Abundance, DetectionLimit, QcFlag,
│       │   │                         #  MeasurementRecord, Platform, Batch
│       │   ├── ingest/
│       │   │   ├── mod.rs            # trait ProteomeIngest
│       │   │   └── olink_explore.rs  # long-CSV NPX parser (v0.1 adapter)
│       │   ├── qc.rs                 # declarative QC rule engine
│       │   ├── matrix.rs             # long ↔ Dube-wide pivoting
│       │   ├── fold_change.rs        # linear-mean / log2-ratio FC computation
│       │   └── sample_id.rs          # SampleIdParser trait + Dube parser
│       └── tests/                    # unit tests on synthetic fixtures
├── example_data/
│   └── dube_heat_2023/
│       ├── 20212016_Dube_NPX_2021-11-30.csv             (raw file 1)
│       ├── 20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv  (raw file 2)
│       ├── filtered_npx/npx/<panel>_npx.csv × 8          (reference)
│       └── fold_changes/fold_changes/<panel>_log2_fc.csv × 8   (reference)
├── docs/
│   └── superpowers/specs/
│       └── 2026-04-14-atman-olink-reproduction-design.md
└── README.md
```

**Two crates by design:**
- `atman-core` is pure library. No `main`, no stdout, no direct filesystem access beyond what its ingest adapter needs. All algorithms operate on types, making them testable in isolation and reusable from notebooks, future CLIs, or other crates.
- `atman` is a thin CLI wrapper: clap dispatch, file I/O, TSV/CSV emission, atomic tempfile-rename writers. Target line count per command file: <200 LOC.

## 4. Core data model

### 4.1 Types (`atman-core::types`)

```rust
/// A proteomics platform identifier. New variants added as adapters ship.
/// v0.1 ships only OlinkExploreNgs as implemented; others are named stubs.
#[non_exhaustive]
pub enum Platform {
    OlinkExploreNgs,     // v0.1 implemented
    OlinkTargetQpcr,     // stub
    SomaScan,            // stub
    MaxQuantLfq,         // stub
    DiannReport,         // stub
    SpectronautReport,   // stub
}

/// Newtype wrapping a platform-specific assay primary key.
/// For Olink: "OID20838". For MS: peptide or protein-group ID.
pub struct AssayId(pub String);

/// Protein identity, catalog-level. Keyed by (Platform, AssayId).
/// Gene symbol is optional and non-unique across panels.
pub struct ProteinIdentity {
    pub platform: Platform,
    pub assay_id: AssayId,
    pub uniprot: Vec<String>,        // may be combined entries like ["Q9UJY5"]
    pub gene_symbol: Option<String>, // may be combined like "MICB_MICA"
    pub panel: Option<String>,
    pub panel_lot: Option<String>,
}

/// Unit-tagged abundance. Type safety against mixing log2-NPX with raw intensity.
pub enum Abundance {
    Log2Npx(f64),
    Log2Intensity(f64),
    Ibaq(f64),
    Raw(f64),
}

/// LOD is per-assay-per-plate for Olink, None for many MS pipelines.
pub struct DetectionLimit(pub Option<f64>);

/// Three-state QC flag. Warning / Fail carry a reason string.
pub enum QcFlag {
    Pass,
    Warn(String),
    Fail(String),
}

/// Batch / run provenance. All fields optional so MS and Olink share the type.
pub struct Batch {
    pub plate: Option<String>,
    pub lot: Option<String>,
    pub run: Option<String>,
}

/// One measurement: a single sample × single assay observation.
/// Primary key: (platform, assay_id, sample_id). Enforced at ingest.
pub struct MeasurementRecord {
    pub platform: Platform,
    pub assay_id: AssayId,
    pub sample_id: String,
    pub abundance: Abundance,
    pub abundance_raw: Abundance,    // original, before QC masking
    pub qc_sample: QcFlag,
    pub qc_assay: QcFlag,
    pub detection_limit: DetectionLimit,
    pub below_lod: bool,
    pub batch: Batch,
    pub dropped_by_qc: bool,
}

/// A biological or control sample.
pub struct Sample {
    pub sample_id: String,
    pub subject_id: Option<String>,  // populated by SampleIdParser
    pub condition: Option<String>,   // "exposure" in Dube nomenclature
    pub is_control: bool,
    pub sample_type: Option<String>, // plasma, serum, csf, tissue_lysate, ...
}
```

### 4.2 Primary key and uniqueness invariants

- **Measurement primary key:** `(Platform, AssayId, SampleId)`. Enforced at the boundary of every ingest adapter. Violation = `IngestError::DuplicatePrimaryKey`, hard fail, non-zero exit.
- **Gene symbol is NOT unique.** `TNF` appears in Inflammation, Inflammation_II, and potentially others with distinct OlinkIDs and distinct measurements. All downstream code keys on `AssayId`, never `gene_symbol`.
- **UniProt is NOT unique.** Same reason as gene symbol.
- **Protein catalog is normalized out of measurements.** `ProteinIdentity` lives in a separate sidecar (`proteins.tsv`) keyed by `(Platform, AssayId)`. Measurements carry `AssayId` only, not duplicated metadata — the raw Olink CSV repeats UniProt/gene/panel/lot 40 times per assay; ingest normalizes.

### 4.3 Canonical formats on disk

**Long TSV (`measurements.tsv`):** the canonical shape for all analytical steps. One row per measurement. Tab-delimited. First line is header.

| column | type | notes |
|---|---|---|
| `sample_id` | string | NPX `SampleID` |
| `assay_id` | string | `OID20838` |
| `abundance` | float \| "" | post-QC value; empty if masked by QC filter |
| `abundance_raw` | float | original NPX value, never masked |
| `abundance_unit` | string | `log2_npx` |
| `qc_sample` | string | `PASS` / `WARN` |
| `qc_assay` | string | `PASS` / `WARN` |
| `detection_limit` | float | LOD for this assay on this plate |
| `below_lod` | int | 0 / 1 |
| `dropped_by_qc` | int | 0 / 1 |
| `plate_id` | string | |
| `panel` | string | |
| `panel_lot` | string | |

Missing cells are represented as empty strings (`""`), never `NA` or `NaN`. This matches the tsvkit convention.

**Dube-wide per-panel CSV (reference reproduction format):** comma-delimited, one file per panel. Header = `Participant,Exposure,SampleID,<assay1>,<assay2>,...` where assays are sorted alphabetically. Rows ordered by `(participant, exposure)`. Empty cells = missing. Used for ingest/QC reproduction check.

**Protein catalog sidecar (`proteins.tsv`):** tab-delimited. One row per unique `(Platform, AssayId)`. Columns: `platform, assay_id, uniprot, gene_symbol, panel, panel_lot`.

**Sample sheet sidecar (`samples.tsv`):** tab-delimited. One row per unique sample_id. Columns: `sample_id, subject_id, condition, is_control, sample_type, plate_id`.

## 5. Ingest trait and adapter

```rust
pub trait ProteomeIngest {
    fn platform(&self) -> Platform;
    fn read(
        &self,
        inputs: &[PathBuf],
        sample_id_parser: &dyn SampleIdParser,
    ) -> Result<IngestOutput, IngestError>;
}

pub struct IngestOutput {
    pub measurements: Vec<MeasurementRecord>,
    pub proteins: Vec<ProteinIdentity>,
    pub samples: Vec<Sample>,
}

pub struct OlinkExploreLongCsv;
impl ProteomeIngest for OlinkExploreLongCsv { /* v0.1 impl */ }

pub trait SampleIdParser {
    fn parse(&self, sample_id: &str) -> Result<ParsedSampleId, ParseError>;
}

pub struct DubeSampleIdParser;  // regex: ^SSNA-(?P<subject>[^-]+)-(?P<condition>PR1|PR2|PT1|PT2)$
```

### 5.1 OlinkExploreLongCsv parse rules
- Delimiter: semicolon.
- Header row is the 14 columns observed in Dube raw CSVs: `SampleID; Index; OlinkID; UniProt; Assay; MissingFreq; Panel; Panel_Lot_Nr; PlateID; QC_Warning; LOD; NPX; Normalization; Assay_Warning`. Missing or reordered columns = hard fail (`IngestError::SchemaMismatch`).
- Each data row yields one `MeasurementRecord` and one `(Platform, AssayId) → ProteinIdentity` entry (deduplicated across rows).
- `Index` is ignored.
- `MissingFreq` is ignored for v0.1 (reserved for future QC rule).
- `Normalization` column is verified to be `"Plate control"` for every row; any other value = `IngestError::UnexpectedNormalization`.
- `NPX` parses as f64. No coercion. Missing = error.
- `LOD` parses as f64. `below_lod = (npx < lod)`.
- Multiple input files are processed sequentially; their measurements merged. Duplicate `(Platform, AssayId, SampleId)` across files = error.

### 5.2 SampleIdParser (Dube mode)
- Regex: `^SSNA-(?P<subject>[^-]+)-(?P<condition>PR1|PR2|PT1|PT2)$`
- Non-matching IDs are flagged as control candidates. If the ID also matches `^CONTROL_SAMPLE_`, it's marked `is_control=true`. Otherwise `IngestError::UnparseableSampleId`.
- `--parser dube` is the v0.1 default. Other parsers slot in as additional trait impls.

### 5.3 IngestError variants
```
SchemaMismatch { file: PathBuf, missing: Vec<String>, extra: Vec<String> }
DuplicatePrimaryKey { platform, assay_id, sample_id, files: [PathBuf; 2] }
InvalidAbundance { file, line, raw: String }
MissingRequiredColumn { file, column }
UnexpectedNormalization { file, line, value: String }
UnparseableSampleId { sample_id }
InconsistentProteinMetadata { assay_id, field, first: String, second: String }
IoError { inner: std::io::Error }
```
No silent coercion. Any error produces no output files and a non-zero exit.

## 6. Subcommand surface (v0.1)

Four subcommands. Linear pipeline. Each reads a file, writes a file, exits non-zero on any error. No hidden state.

```
raw NPX CSV ─▶ ingest ─▶ measurements.tsv (long) ─▶ qc ─▶ qc_measurements.tsv ─┐
                          proteins.tsv                      qc_report.tsv       │
                          samples.tsv                                            │
                                                                                 ▼
                                                                              matrix ─▶ <panel>_npx.csv × 8
                                                                                 │
                                                                                 ▼
                                                                          fold-change ─▶ <panel>_log2_fc.csv × 8
```

### 6.1 `atman ingest`
```
atman ingest \
    --platform olink-explore-ngs \
    --parser dube \
    --output-dir <dir> \
    <input1.csv> [<input2.csv> ...]
```
Outputs: `measurements.tsv`, `proteins.tsv`, `samples.tsv`, `ingest_report.tsv`.

`ingest_report.tsv` columns: `input_file, rows_read, rows_ingested, rows_rejected, reject_reason, control_samples_detected, biological_samples_detected, panels_detected`.

### 6.2 `atman qc`
```
atman qc \
    --input-dir <dir> \
    --rule dube \
    --output-dir <dir>
```
Applies one or more named QC rules. v0.1 implements `dube`:
- Set `abundance = ""` and `dropped_by_qc = 1` for any row where `qc_sample != "PASS"` or `qc_assay != "PASS"`.
- Original `abundance_raw` is preserved untouched.
- Control samples are **NOT dropped**. They flow through to the output with their measurements intact. This matches Dube's published filtered NPX behavior (verified empirically — controls appear at the tail of every filtered_npx panel file). An opt-in `--drop-controls` flag is deferred to v0.2 when DE analysis needs it.

Outputs: `qc_measurements.tsv`, `qc_report.tsv`.

`qc_report.tsv` columns: `rule, dropped_count, masked_count, affected_samples, affected_assays, pass_count`.

### 6.3 `atman matrix`
```
atman matrix \
    --input-dir <dir> \
    --format dube-wide \
    --split-by panel \
    --output-dir <dir>
```
Produces 8 per-panel CSVs in Dube-wide format:
- One file per unique value of the raw NPX `Panel` column (values observed: `Cardiometabolic`, `Cardiometabolic_II`, `Inflammation`, `Inflammation_II`, `Neurology`, `Neurology_II`, `Oncology`, `Oncology_II`).
- Filename: ASCII lowercase of the panel value + `_npx.csv`. No space→underscore translation (no spaces observed). `Cardiometabolic_II` → `cardiometabolic_ii_npx.csv`.
- Column order: `Participant, Exposure, SampleID, <assays sorted alphabetically by Unicode codepoint>`. Verified empirically on all 8 Dube reference files.
- Row order: **two-pass**. (1) Rows whose `SampleID` is parseable by the active SampleIdParser (i.e., biological samples) come first, sorted by `(Participant, Exposure)` lexicographically. (2) Rows whose `SampleID` is NOT parseable (controls) come second, appended in ingestion order (the order they appeared in the input NPX file). Controls have empty `Participant` and `Exposure` cells.
- Cell value: `abundance` from `qc_measurements.tsv`, or empty string if missing / masked.

The only `--format` in v0.1 is `dube-wide`. A future `wide` format (tab-delimited, `assay_id` rows × `sample_id` columns, karnaplot-friendly) is planned for v0.2.

### 6.4 `atman fold-change`
```
atman fold-change \
    --input-dir <dir> \
    --groups "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1" \
    --group-column exposure \
    --participant-column participant \
    --split-by panel \
    --output-dir <dir>
```
Produces 8 per-panel log2 FC CSVs:
- Filename: `<lowercased_panel_name>_log2_fc.csv`.
- Column order: `Assay, <comparison1>, <comparison2>, ...` in the order given by `--groups`.
- One row per assay in the panel (all assays, not only non-missing).
- Cell value: log2 fold change per the formula in §7.4.

## 7. Transformation algorithms

### 7.1 QC filter (Dube rule)
```
for each row in measurements:
    if row.qc_sample != Pass or row.qc_assay != Pass:
        row.abundance     = MISSING
        row.dropped_by_qc = true
    # abundance_raw is never mutated
```

### 7.2 Control sample handling
A sample is classified as `is_control = true` iff the active `SampleIdParser` fails to parse its `SampleID`. For the Dube parser this is equivalent to the sample not matching the `^SSNA-(?P<subject>[^-]+)-(?P<condition>PR1|PR2|PT1|PT2)$` regex. The `^CONTROL_SAMPLE_` prefix is a strong hint but not authoritative — unparseable-is-control is the single rule.

**v0.1 keeps controls in the output.** The `qc` command does not drop them, and the `matrix` pivot appends them after biological samples per §7.3. This matches Dube's published filtered NPX, which includes 1-2 controls per panel at the tail of the file.

### 7.3 Dube-wide pivot
```
input:  qc_measurements.tsv joined with samples.tsv (to get participant, exposure, is_control, ingest_order)
group:  by panel
for each panel p:
    bio_rows  = rows in p where is_control=false, sorted by (participant, exposure) lex
    ctl_rows  = rows in p where is_control=true,  sorted by ingest_order
    rows      = bio_rows ++ ctl_rows
    columns   = unique assays in p, sorted by Unicode codepoint order
    cell(r,c) = abundance where (sample_id,assay) matches; "" if missing or masked
    header    = Participant, Exposure, SampleID, <assays...>
    for control rows: Participant = "", Exposure = ""
    emit:     <ascii_lowercase(panel)>_npx.csv, comma-delimited
```

Assay column sort is **Unicode codepoint**, not locale-aware. Verified empirically on all 8 Dube reference files.

Row order is **biological-then-control**. This matters because empty-string `Participant` would sort to the top under naive lex ordering, but Dube places controls at the tail. Verified empirically on all 8 reference files.

Ingest must preserve an `ingest_order` integer per row so that control samples' tail ordering is reproducible. For Dube's data, controls from Explore 1536 files (2 each) go into the non-`_II` panel outputs; controls from the Explore Expansion file (1 each) go into the `_II` panel outputs. The per-file control count differs; the pivot does not hardcode 1-or-2 — it just preserves whatever controls the ingested panel rows contain.

### 7.4 Log2 fold change
```
input: qc_measurements.tsv joined with samples.tsv
group: by panel
for each (panel, assay, comparison) where comparison = "A-B":
    values_a = [abundance for each (participant, exposure=A, assay) where abundance != MISSING]
    values_b = [abundance for each (participant, exposure=B, assay) where abundance != MISSING]
    if values_a.is_empty() or values_b.is_empty():
        fc = MISSING
    else:
        mean_a_linear = mean(2^v for v in values_a)
        mean_b_linear = mean(2^v for v in values_b)
        fc = log2(mean_a_linear / mean_b_linear)
output: per-panel CSV, header = "Assay,<comp1>,<comp2>,..."
```
All arithmetic in IEEE 754 f64. No intermediate rounding. Output formatted as full f64 precision (17 digits is fine; tolerance for reproduction diff is 1e-4).

**Known ambiguity to resolve during implementation:** whether Dube's mean is over the 10 declared participants (with missing → error / 0) or over available participants only. The paper states "Missing data were not imputed," which we interpret as *available participants only*. We verify this empirically by diffing against Dube's published fold change files. If the interpretation is wrong, the diff fails on the first row with any missingness, we adjust to the correct rule, and re-test.

## 8. Reproduction test (the acceptance criterion)

The `atman/tests/dube_reproduction.rs` integration test does:

```
1. Run `atman ingest --platform olink-explore-ngs --parser dube \
       --output-dir tmp/ \
       example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv \
       example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv`
2. Run `atman qc --input-dir tmp/ --rule dube --output-dir tmp/`
3. Run `atman matrix --input-dir tmp/ --format dube-wide --split-by panel --output-dir tmp/`
4. Run `atman fold-change --input-dir tmp/ --groups "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1" \
       --group-column exposure --participant-column participant --split-by panel --output-dir tmp/`
5. For each panel p in [cardiometabolic, inflammation, neurology, oncology, cardiometabolic_ii,
                        inflammation_ii, neurology_ii, oncology_ii]:
   a. Diff tmp/<p>_npx.csv against example_data/.../filtered_npx/npx/<p>_npx.csv
      — exact string match per cell (empty string vs empty string ok).
   b. Diff tmp/<p>_log2_fc.csv against example_data/.../fold_changes/fold_changes/<p>_log2_fc.csv
      — Assay column: exact match.
      — Numeric columns: |delta| <= 1e-4, or both missing.
6. Assert zero mismatches across all 16 reference files.
```

The diff harness is a helper in `atman/tests/diff.rs` (or an internal test-only module). It aligns by row key (first column) and reports up to 10 mismatches before failing with a useful diagnostic.

**Tolerance rationale:**
- **Filtered NPX:** exact string match. Values are NPX values we copy through unchanged when QC passes, and empty strings when QC fails. No arithmetic, no drift. A mismatch means we got the filter rule wrong, the pivot wrong, or the column order wrong.
- **Fold change:** 1e-4 tolerance. The computation is `log2(mean(2^x) / mean(2^y))` in f64; if Dube used R's double precision the result should match to ~15 decimal places. 1e-4 is 11 orders of magnitude looser than that, reserved for any R-vs-Rust implementation detail we can't foresee. A real bug would show deltas orders of magnitude larger than 1e-4.

## 9. Error handling

Error surface is narrow and explicit:

- **Parse errors** (`SchemaMismatch`, `InvalidAbundance`, `UnexpectedNormalization`, `MissingRequiredColumn`, `UnparseableSampleId`): produced at ingest. No output files written. Non-zero exit. Error message prints the file path, line number, and the raw offending value.
- **Invariant violations** (`DuplicatePrimaryKey`, `InconsistentProteinMetadata`): produced at ingest. Same behavior.
- **IO errors** (`IoError`): produced anywhere. Propagated with context; no silent failure.
- **QC anomalies:** NOT errors. A sample with `WARN` QC is not an error — it's data. QC rules mask values declaratively and log what they did.

Every subcommand writes output files via `tempfile::NamedTempFile::persist`, so a mid-run failure never leaves a half-written file next to a previous run's output.

Every subcommand's stderr on success contains a one-line summary (`ingest: 122135 rows, 2943 assays, 40 samples, 0 rejected`). Verbose details go to the `*_report.tsv` sidecar.

## 10. Testing strategy

### 10.1 Unit tests (`atman-core/tests/`)
Synthetic fixtures, one algorithm each:

| test file | what it tests |
|---|---|
| `ingest_olink.rs` | 1-row / 2-row / 10-row synthetic NPX CSV → parsed `IngestOutput` matches expectations |
| `ingest_errors.rs` | each `IngestError` variant: schema mismatch, duplicate PK, invalid NPX, unparseable SampleID, unexpected normalization |
| `qc_dube.rs` | dube rule on synthetic: PASS/PASS → kept, WARN/PASS → masked, PASS/WARN → masked, WARN/WARN → masked |
| `pivot_dube.rs` | 2-panel, 3-sample, 4-assay synthetic → correct pivot, alphabetical assay order, missing as empty |
| `fold_change.rs` | hand-computed FC on 3-assay × 2-participant × 2-exposure fixture → matches expected |
| `sample_id_dube.rs` | valid IDs parse; control IDs classify; malformed IDs error |

Target: ~20-25 unit tests, each <50 LOC, runtime <1s total.

### 10.2 Property tests (`atman-core/tests/properties.rs`)
Invariants that hold for any valid input:
1. No measurement appears twice in long output (PK uniqueness).
2. Every `(participant, exposure)` in input appears as a row in each panel's Dube-wide output.
3. For any `exposure`, `fold_change(exposure, exposure) == 0.0` across all assays (self-vs-self is zero).
4. If all values are missing in either group, FC is missing.
5. Ingest → wide → FC is commutative under assay reordering in the input CSV (output must be identical).
6. QC masking is idempotent: running the `dube` rule twice gives the same result as running it once.

Tool: `proptest` for randomized generation, or hand-rolled for small cases.

### 10.3 Reproduction integration test (`atman/tests/dube_reproduction.rs`)
Described in §8. The single test that defines "v0.1 done."

### 10.4 Runtime budget
Full `cargo test --workspace`: <30 seconds on a Mac Studio M-series. If exceeded, investigate — the reproduction test is the dominant cost and should run in <5 seconds given the dataset size (120K NPX rows).

## 11. Dependencies

Minimize. Every dep is a liability.

| crate | purpose | rationale |
|---|---|---|
| `clap` (derive) | CLI parsing | standard in karna, 34 CLIs |
| `csv` | CSV/TSV reading/writing | standard, handles quoting, delimiter config |
| `serde` + `serde_derive` | struct (de)serialization | standard |
| `thiserror` | error type derivation | standard for library error enums |
| `anyhow` | binary error wrapping | atman binary only, not atman-core |
| `regex` | SampleID parsing | stable, no alternatives |
| `tempfile` | atomic writes | standard |
| `proptest` (dev) | property tests | dev-only |

**Explicitly NOT taking:**
- No async runtime. No tokio. No rayon in v0.1 (the dataset is small enough).
- No HDF5, no parquet, no arrow. TSV/CSV is the wire format.
- No SQL in atman-core. atman binary may gain rusqlite in v0.2 for karnadata integration.
- No statistics crate. `statrs`, `ndarray`, `nalgebra`, `polars` all deferred until v0.2 needs them.

## 12. What v0.1 does NOT implement — cross-references

For each of these, the v0.2+ section below names the milestone and the rough approach.

- Karnadata extension (§13.1)
- Differential abundance testing (§13.2)
- Bridge normalization (§13.3)
- Protein set enrichment (§13.4)
- Non-Olink ingest adapters (§13.5)
- karnaplot integration (§13.6)
- Claim emission into knowledge.db (§13.7)

## 13. Roadmap beyond v0.1

### 13.1 v0.2 — Karnadata extension
Add `proteome_sample` and `proteome_assay` tables to karnadata. Add sibling command `karnadata validate-proteome --strict`. Wire into `atman ingest` as a pre-flight hook that hard-fails if the cohort is not registered. Mandatory metadata fields: `subject_id`, `sample_type`, `platform`, `panels`, `plate_id`, `condition`, `is_control`, `verification_status`. Subject-map TSV authored by a human enforces the biological-replicate-as-statistical-unit discipline at registration time.

### 13.2 v0.2 — Differential abundance
`atman de` subcommand. LOD-aware, paired-sample Wilcoxon by default; unpaired when samples have distinct subject_ids. Multiple testing correction: BH-FDR. Tests operate on the sample-level biological replicate, not individual measurements — enforced by requiring `--paired-by subject_id` or `--independent` explicitly.

### 13.3 v0.2 — Bridge normalization
`atman normalize --method bridge --bridge-samples <file>`. Requires a file listing sample IDs shared across plates/batches; computes per-assay median shift between plates and applies it. The v0.1 `normalize` command is deferred entirely rather than shipped as a type-checked no-op.

### 13.4 v0.2 — Protein set enrichment
Handoff to the existing `genesets` Rust CLI via TSV. No new logic in atman beyond emitting a ranked assay list and a mapping to gene symbols. MSigDB + Reactome + GO-BP as starting sets.

### 13.5 v0.3 — Additional ingest adapters
`OlinkTargetQpcr` (different header schema), `SomaScan` (ADAT format), `MaxQuantLfq` (proteinGroups.txt), `DiannReport` (report.tsv). Each is a new `impl ProteomeIngest`. `atman-core` types should not need to change.

### 13.6 v0.3 — karnaplot integration
Olink-specific plot types via the existing `karnaplot` tool: panel-grouped volcano, per-subject trajectory line plots, panel heatmap with hierarchical clustering, per-assay scatter against a clinical covariate. All via TSV handoff.

### 13.7 v0.4 — Evidence pipeline integration
Claim emission into `knowledge.db`. Protein-level claims flow through Crucible (adjudication), Oracle (external validation against DepMap, CRISPRbrain, LINCS, plus new protein-relevant sources), Atlas (consensus networks), and into the vigil pipeline. Requires deciding how proteomic claims map to the existing claim schema; likely a new claim type.

## 14. Open questions

1. **Fold change missingness rule:** paper says "missing not imputed" but doesn't state the mean's denominator behavior with missing participants. Assumed: mean over available participants only. Verify empirically against Dube's FC files during §8 test run.
2. **f64 vs R numeric drift:** Dube's published FC files have ~15-digit precision. If we're off by more than 1e-12, that's a bug, not drift. 1e-4 tolerance is the coarse outer bound.
3. **Assay column sort in the `_II` panels:** we're assuming Unicode codepoint sort produces the same order Dube used. Empirically verified for `inflammation_npx.csv`; to confirm for all 8 panels via the reproduction test.
4. **Multiple raw files / panel count:** Dube ships 2 raw CSVs containing 8 panels total. Our ingest merges them and `--split-by panel` produces 8 outputs. If a future study splits differently (more files, fewer panels, or vice versa), the ingest is unaffected — only the number of output files changes.

## 15. Acceptance criteria

v0.1 is done when:

1. `cargo test --workspace` passes on a clean checkout (Rust stable, Mac Studio M-series).
2. The reproduction integration test passes: all 16 Dube reference files diff-clean.
3. `cargo clippy --workspace --all-targets -- -D warnings` clean.
4. `cargo fmt --check` clean.
5. No crate has a `main.rs` that panics on any path reachable from the CLI surface; errors surface as `Result` through the binary's `main`.
6. Total runtime of the reproduction test < 30 seconds on a Mac Studio.
7. README contains one-paragraph summary + one-command reproduction recipe.

Passing these seven criteria = v0.1 ships. Not passing any one = v0.1 is not done.

---

**End of design.**
