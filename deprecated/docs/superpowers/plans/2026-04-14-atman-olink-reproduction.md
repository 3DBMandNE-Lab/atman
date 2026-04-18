# Atman v0.1 — Olink Explore Reproduction Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the v0.1 atman engine: a Rust workspace with two crates (`atman-core` library, `atman` binary) that ingests Olink Explore NGS long-CSV NPX files and reproduces the Dube et al. published filtered NPX + log2 fold change outputs byte-for-byte (strings) and within 1e-4 (fold-change numerics).

**Architecture:** Two Rust crates. `atman-core` is a pure library holding platform-agnostic types (`Platform`, `AssayId`, `MeasurementRecord`, …) and all transformation algorithms (ingest, QC, pivot, fold change). `atman` is a thin binary providing four subcommands (`ingest`, `qc`, `matrix`, `fold-change`) that read TSV/CSV from disk, invoke the library, and write TSV/CSV to disk with atomic tempfile renames. The reproduction integration test is the single acceptance gate.

**Tech stack:** Rust stable · `clap` (derive) · `csv` · `serde` + `serde_derive` · `thiserror` · `anyhow` · `regex` · `tempfile` · `proptest` (dev).

**Spec:** `docs/superpowers/specs/2026-04-14-atman-olink-reproduction-design.md` (v2, post-correction commit `62b14fe`). Read §4 (types), §5 (ingest trait), §6 (subcommand surface), §7 (algorithms), §8 (reproduction test) before starting any task.

**Reference data:** `example_data/dube_heat_2023/` contains the 2 raw NPX CSVs and 16 published reference files (8 filtered NPX + 8 fold change). Already downloaded. Do not re-download.

---

## File Structure

The plan builds these files, in the order listed. Each file has one responsibility.

### Workspace root
```
Cargo.toml                              # workspace manifest
.gitignore                              # /target, *.tmp
README.md                               # one-paragraph summary + reproduction recipe
```

### `crates/atman-core/` (library)
```
Cargo.toml
src/
├── lib.rs                              # public re-exports + crate-level docs
├── types.rs                            # Platform, AssayId, ProteinIdentity, Abundance,
│                                       #  DetectionLimit, QcFlag, Batch,
│                                       #  MeasurementRecord, Sample
├── errors.rs                           # IngestError enum (thiserror)
├── sample_id.rs                        # SampleIdParser trait + DubeSampleIdParser + ParsedSampleId
├── ingest/
│   ├── mod.rs                          # ProteomeIngest trait, IngestOutput struct
│   └── olink_explore.rs                # OlinkExploreLongCsv adapter
├── qc.rs                               # dube rule (mask abundance when QC != PASS)
├── matrix.rs                           # long → Dube-wide pivot
└── fold_change.rs                      # log2 FC computation
tests/
├── types_roundtrip.rs                  # basic type construction / serde
├── sample_id_dube.rs                   # parser cases incl. control detection
├── ingest_olink.rs                     # synthetic NPX → IngestOutput happy path
├── ingest_errors.rs                    # each IngestError variant
├── qc_dube.rs                          # PASS/WARN combinatorics
├── pivot_dube.rs                       # 2-panel × 3-sample × 4-assay fixture
├── fold_change.rs                      # hand-computed FC check
└── properties.rs                       # proptest invariants
```

### `crates/atman/` (binary)
```
Cargo.toml
src/
├── main.rs                             # clap dispatch, exit codes
├── io.rs                               # atomic tempfile writer, long TSV read/write, samples.tsv / proteins.tsv
└── commands/
    ├── mod.rs                          # re-export
    ├── ingest.rs                       # atman ingest
    ├── qc.rs                           # atman qc
    ├── matrix.rs                       # atman matrix
    └── fold_change.rs                  # atman fold-change
tests/
├── common/
│   └── mod.rs                          # CSV diff harness
└── dube_reproduction.rs                # the acceptance test (runs the 4-stage pipeline + diffs all 16 files)
```

---

## Task Dependency Overview

Tasks build bottom-up:
- **1–2:** Workspace scaffold and atman-core skeleton.
- **3–6:** Core types → errors → sample parser → ingest trait.
- **7:** Olink ingest adapter (the biggest single piece).
- **8–10:** QC → pivot → fold change algorithms.
- **11:** Property tests (invariant check on the full library).
- **12–13:** atman binary scaffold + IO helpers.
- **14–17:** Four CLI commands.
- **18:** CSV diff harness.
- **19:** Reproduction integration test (the acceptance gate).
- **20:** README + final lint/fmt/test pass.

A later task cannot be started until its dependencies compile and pass their tests. Commit after every task.

---

## Task 1: Workspace Scaffold and `.gitignore`

**Files:**
- Create: `Cargo.toml`
- Create: `.gitignore`
- Create: `crates/` (directory, will be populated by later tasks)

- [ ] **Step 1: Create workspace `Cargo.toml`**

```toml
[workspace]
members = ["crates/atman-core", "crates/atman"]
resolver = "2"

[workspace.package]
version = "0.1.0"
edition = "2021"
rust-version = "1.75"
license = "MIT OR Apache-2.0"
authors = ["Kevin Joseph"]
repository = "https://github.com/kevinjoseph/atman"

[workspace.dependencies]
anyhow = "1"
clap = { version = "4", features = ["derive"] }
csv = "1"
regex = "1"
serde = { version = "1", features = ["derive"] }
tempfile = "3"
thiserror = "1"
proptest = "1"
```

- [ ] **Step 2: Create `.gitignore`**

```
/target
**/*.rs.bk
*.tmp
.DS_Store
```

- [ ] **Step 3: Create empty `crates/` directory (no content yet)**

```bash
mkdir -p crates
```

- [ ] **Step 4: Commit**

```bash
git add Cargo.toml .gitignore
git commit -m "chore: workspace scaffold"
```

---

## Task 2: `atman-core` Crate Skeleton

**Files:**
- Create: `crates/atman-core/Cargo.toml`
- Create: `crates/atman-core/src/lib.rs`

- [ ] **Step 1: Create `crates/atman-core/Cargo.toml`**

```toml
[package]
name = "atman-core"
version.workspace = true
edition.workspace = true
license.workspace = true

[dependencies]
csv = { workspace = true }
regex = { workspace = true }
serde = { workspace = true }
thiserror = { workspace = true }

[dev-dependencies]
proptest = { workspace = true }
tempfile = { workspace = true }
```

- [ ] **Step 2: Create `crates/atman-core/src/lib.rs`**

```rust
//! atman-core: platform-agnostic primitives and algorithms for proteomics data.
//!
//! v0.1 implements exactly one ingest adapter (Olink Explore NGS long CSV) and
//! the transformations required to reproduce the Dube et al. Scientific Data 2023
//! published filtered NPX + log2 fold-change files. See
//! `docs/superpowers/specs/2026-04-14-atman-olink-reproduction-design.md`
//! for the full design.

pub mod errors;
pub mod fold_change;
pub mod ingest;
pub mod matrix;
pub mod qc;
pub mod sample_id;
pub mod types;

pub use errors::IngestError;
pub use ingest::{IngestOutput, ProteomeIngest};
pub use sample_id::{DubeSampleIdParser, ParsedSampleId, SampleIdParser};
pub use types::{
    Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, ProteinIdentity,
    QcFlag, Sample,
};
```

- [ ] **Step 3: Create empty module files so lib.rs compiles**

Create stub files so the compilation succeeds:

```bash
mkdir -p crates/atman-core/src/ingest crates/atman-core/tests
touch crates/atman-core/src/types.rs \
      crates/atman-core/src/errors.rs \
      crates/atman-core/src/sample_id.rs \
      crates/atman-core/src/ingest/mod.rs \
      crates/atman-core/src/qc.rs \
      crates/atman-core/src/matrix.rs \
      crates/atman-core/src/fold_change.rs
```

At this point each file is empty, so the re-exports in `lib.rs` will fail. We'll populate them in the next tasks. For now, temporarily comment out the re-exports block in `lib.rs` and keep only `pub mod` lines. Uncomment as types are added.

- [ ] **Step 4: Verify the empty crate compiles**

Replace `lib.rs` temporarily with this minimal version so Task 2 ends green:

```rust
//! atman-core: platform-agnostic primitives and algorithms for proteomics data.
pub mod errors;
pub mod fold_change;
pub mod ingest;
pub mod matrix;
pub mod qc;
pub mod sample_id;
pub mod types;
```

Run:

```bash
cargo build -p atman-core
```

Expected: success with warnings about unused empty modules.

- [ ] **Step 5: Commit**

```bash
git add crates/atman-core
git commit -m "feat(atman-core): crate skeleton"
```

---

## Task 3: Core Types

**Files:**
- Modify: `crates/atman-core/src/types.rs`
- Modify: `crates/atman-core/src/lib.rs` (add re-exports)
- Test: `crates/atman-core/tests/types_roundtrip.rs`

Rationale: types are used by every other module. Build them first, then everything compiles against a stable foundation.

- [ ] **Step 1: Write the failing test**

Create `crates/atman-core/tests/types_roundtrip.rs`:

```rust
use atman_core::{
    Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, ProteinIdentity,
    QcFlag, Sample,
};

fn make_record() -> MeasurementRecord {
    MeasurementRecord {
        platform: Platform::OlinkExploreNgs,
        assay_id: AssayId("OID20838".to_string()),
        gene_symbol: Some("GGA1".into()),
        sample_id: "SSNA-001B-PR1".to_string(),
        abundance: Abundance::Log2Npx(0.1234),
        abundance_raw: Abundance::Log2Npx(0.1234),
        npx_source_str: "0.1234".to_string(),
        qc_sample: QcFlag::Pass,
        qc_assay: QcFlag::Pass,
        detection_limit: DetectionLimit(Some(0.2178)),
        below_lod: true,
        batch: Batch {
            plate: Some("plate1".into()),
            lot: Some("B04414".into()),
            run: None,
        },
        dropped_by_qc: false,
        ingest_order: 0,
        panel: Some("Neurology".into()),
    }
}

#[test]
fn construct_measurement_record() {
    let rec = make_record();
    assert_eq!(rec.assay_id.0, "OID20838");
    assert_eq!(rec.gene_symbol.as_deref(), Some("GGA1"));
    assert_eq!(rec.npx_source_str, "0.1234");
    assert!(matches!(rec.abundance, Abundance::Log2Npx(v) if (v - 0.1234).abs() < 1e-12));
    assert!(matches!(rec.qc_sample, QcFlag::Pass));
}

#[test]
fn effective_abundance_accessor_honors_qc() {
    let mut rec = make_record();
    assert_eq!(rec.effective_abundance(), Some(0.1234));
    rec.dropped_by_qc = true;
    assert_eq!(rec.effective_abundance(), None);
}

#[test]
fn qc_flag_warn_carries_reason() {
    let f = QcFlag::Warn("plate deviation".to_string());
    match f {
        QcFlag::Warn(r) => assert_eq!(r, "plate deviation"),
        _ => panic!("expected Warn"),
    }
}

#[test]
fn protein_identity_allows_multi_uniprot() {
    let p = ProteinIdentity {
        platform: Platform::OlinkExploreNgs,
        assay_id: AssayId("OID30000".to_string()),
        uniprot: vec!["P01375".into(), "P01376".into()],
        gene_symbol: Some("MICB_MICA".into()),
        panel: Some("Inflammation".into()),
        panel_lot: Some("B04414".into()),
    };
    assert_eq!(p.uniprot.len(), 2);
    assert_eq!(p.gene_symbol.as_deref(), Some("MICB_MICA"));
}

#[test]
fn sample_struct_supports_control() {
    let s = Sample {
        sample_id: "CONTROL_SAMPLE_US_CS_AS_2-1".into(),
        subject_id: None,
        condition: None,
        is_control: true,
        sample_type: Some("plasma".into()),
        ingest_order: 40,
    };
    assert!(s.is_control);
}
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
cargo test -p atman-core --test types_roundtrip
```

Expected: compilation errors — `Platform`, `AssayId`, etc. are unresolved.

- [ ] **Step 3: Implement `crates/atman-core/src/types.rs`**

```rust
//! Core data model. Platform-agnostic. No I/O, no CSV, no filesystem.

use serde::{Deserialize, Serialize};

/// A proteomics platform identifier. v0.1 implements only `OlinkExploreNgs`.
/// Other variants are named stubs whose ingest adapters will be added in v0.3+.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Platform {
    OlinkExploreNgs,
    OlinkTargetQpcr,
    SomaScan,
    MaxQuantLfq,
    DiannReport,
    SpectronautReport,
}

impl Platform {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::OlinkExploreNgs => "olink_explore_ngs",
            Self::OlinkTargetQpcr => "olink_target_qpcr",
            Self::SomaScan => "somascan",
            Self::MaxQuantLfq => "maxquant_lfq",
            Self::DiannReport => "diann_report",
            Self::SpectronautReport => "spectronaut_report",
        }
    }
}

/// Newtype wrapping a platform-specific assay primary key.
/// Olink: `OID20838`. MS: peptide or protein-group ID.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct AssayId(pub String);

/// Protein catalog entry, keyed by `(Platform, AssayId)`. Gene symbol and UniProt
/// are non-unique across panels — they are metadata, not keys.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ProteinIdentity {
    pub platform: Platform,
    pub assay_id: AssayId,
    pub uniprot: Vec<String>,
    pub gene_symbol: Option<String>,
    pub panel: Option<String>,
    pub panel_lot: Option<String>,
}

/// Unit-tagged abundance. Type-system defense against mixing log2 NPX with raw
/// intensity. v0.1 only emits `Log2Npx` but the enum is ready for MS units.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Abundance {
    Log2Npx(f64),
    Log2Intensity(f64),
    Ibaq(f64),
    Raw(f64),
}

impl Abundance {
    /// Returns the underlying f64 regardless of unit. Useful for the fold-change
    /// path which already knows the unit is log2.
    pub fn as_f64(&self) -> f64 {
        match self {
            Self::Log2Npx(v) | Self::Log2Intensity(v) | Self::Ibaq(v) | Self::Raw(v) => *v,
        }
    }
}

/// Limit of detection. `None` for platforms without a per-assay LOD concept.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DetectionLimit(pub Option<f64>);

/// Three-state QC flag. `Warn` and `Fail` carry a reason string.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum QcFlag {
    Pass,
    Warn(String),
    Fail(String),
}

impl QcFlag {
    pub fn is_pass(&self) -> bool {
        matches!(self, Self::Pass)
    }
    pub fn as_str(&self) -> &str {
        match self {
            Self::Pass => "PASS",
            Self::Warn(_) => "WARN",
            Self::Fail(_) => "FAIL",
        }
    }
}

/// Batch / run provenance. All fields optional for cross-platform use.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Batch {
    pub plate: Option<String>,
    pub lot: Option<String>,
    pub run: Option<String>,
}

/// One sample × one assay measurement. Primary key is `(platform, assay_id, sample_id)`.
///
/// We carry BOTH the canonical f64 value (`abundance`) and the original
/// NPX source string (`npx_source_str`). The f64 drives all numeric paths
/// (QC logic, fold change); the source string is written verbatim to the
/// Dube-wide CSV to guarantee byte-equivalent reproduction — no f64
/// parse+reformat round-trip, no R-vs-Rust decimal drift, no surprises on
/// values like `-0.0000` or ties-to-even rounding.
///
/// `gene_symbol` is carried here (duplicated across every measurement for a
/// given `assay_id`) so the Dube-wide pivot and log2 fold-change row key
/// can use the **gene symbol** (e.g. `ACTN4`, `HLA-DRA`, `MICB_MICA`) —
/// which is what Dube's reference CSVs use as column and row keys — without
/// having to thread the protein catalog through every transformation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MeasurementRecord {
    pub platform: Platform,
    pub assay_id: AssayId,
    /// Gene symbol from the raw Olink `Assay` column. Populated by the
    /// ingest adapter. Used as the wide-format column key and the
    /// fold-change row key, not as a primary key.
    pub gene_symbol: Option<String>,
    pub sample_id: String,
    pub abundance: Abundance,
    pub abundance_raw: Abundance,
    /// Exact source NPX string as it appeared in the raw CSV, e.g. `"0.0980"`
    /// or `"-0.1219"`. Empty for rows that never had a numeric source.
    pub npx_source_str: String,
    pub qc_sample: QcFlag,
    pub qc_assay: QcFlag,
    pub detection_limit: DetectionLimit,
    pub below_lod: bool,
    pub batch: Batch,
    pub dropped_by_qc: bool,
    /// Monotonic ingest-order index assigned at parse time. Used as the stable
    /// tiebreak for control-sample row order in the Dube-wide pivot (§7.3).
    pub ingest_order: u64,
    /// Panel name from the raw Olink `Panel` column (e.g., `Inflammation_II`).
    /// Populated by ingest; consumed by the matrix pivot to group rows.
    pub panel: Option<String>,
}

impl MeasurementRecord {
    /// Returns `None` when the record has been masked by QC, otherwise the f64
    /// abundance. Every consumer that reads `abundance` for computation should
    /// go through this accessor so a future refactor that eliminates the
    /// `dropped_by_qc` flag doesn't silently use stale values.
    pub fn effective_abundance(&self) -> Option<f64> {
        if self.dropped_by_qc {
            None
        } else {
            Some(self.abundance.as_f64())
        }
    }
}

/// A biological or control sample. Fields like `subject_id` and `condition` are
/// populated by a `SampleIdParser`; unparseable IDs produce `is_control=true`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Sample {
    pub sample_id: String,
    pub subject_id: Option<String>,
    pub condition: Option<String>,
    pub is_control: bool,
    pub sample_type: Option<String>,
    /// First-seen ingest order across all inputs.
    pub ingest_order: u64,
}
```

- [ ] **Step 4: Re-run test to verify it passes**

```bash
cargo test -p atman-core --test types_roundtrip
```

Expected: 4 tests passed.

- [ ] **Step 5: Commit**

```bash
git add crates/atman-core/src/types.rs crates/atman-core/tests/types_roundtrip.rs
git commit -m "feat(atman-core): core types (Platform, AssayId, MeasurementRecord, …)"
```

---

## Task 4: `IngestError` Enum

**Files:**
- Modify: `crates/atman-core/src/errors.rs`

- [ ] **Step 1: Write the failing test** (inline in `errors.rs`)

Add to `errors.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn schema_mismatch_display_contains_details() {
        let e = IngestError::SchemaMismatch {
            file: PathBuf::from("/tmp/test.csv"),
            missing: vec!["NPX".to_string()],
            extra: vec![],
        };
        let s = format!("{}", e);
        assert!(s.contains("/tmp/test.csv"));
        assert!(s.contains("NPX"));
    }

    #[test]
    fn duplicate_primary_key_display() {
        let e = IngestError::DuplicatePrimaryKey {
            platform: "olink_explore_ngs".into(),
            assay_id: "OID20838".into(),
            sample_id: "SSNA-001B-PR1".into(),
        };
        let s = format!("{}", e);
        assert!(s.contains("OID20838"));
        assert!(s.contains("SSNA-001B-PR1"));
    }
}
```

- [ ] **Step 2: Run to verify fail**

```bash
cargo test -p atman-core --lib errors
```

Expected: compile error — `IngestError` does not exist.

- [ ] **Step 3: Implement `crates/atman-core/src/errors.rs`**

```rust
//! Error types for atman-core.
use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum IngestError {
    #[error(
        "schema mismatch in {file}: missing columns {missing:?}, unexpected columns {extra:?}"
    )]
    SchemaMismatch {
        file: PathBuf,
        missing: Vec<String>,
        extra: Vec<String>,
    },

    #[error(
        "duplicate primary key: platform={platform} assay_id={assay_id} sample_id={sample_id}"
    )]
    DuplicatePrimaryKey {
        platform: String,
        assay_id: String,
        sample_id: String,
    },

    #[error("invalid abundance value {raw:?} at {file}:{line}")]
    InvalidAbundance {
        file: PathBuf,
        line: usize,
        raw: String,
    },

    #[error("missing required column {column} in {file}")]
    MissingRequiredColumn { file: PathBuf, column: String },

    #[error("unexpected normalization value {value:?} at {file}:{line} (expected 'Plate control')")]
    UnexpectedNormalization {
        file: PathBuf,
        line: usize,
        value: String,
    },

    #[error("unparseable sample id {sample_id:?}")]
    UnparseableSampleId { sample_id: String },

    #[error("inconsistent protein metadata for assay_id {assay_id}: field={field} first={first:?} second={second:?}")]
    InconsistentProteinMetadata {
        assay_id: String,
        field: String,
        first: String,
        second: String,
    },

    #[error("csv parse error in {file}: {source}")]
    CsvError {
        file: PathBuf,
        #[source]
        source: csv::Error,
    },

    #[error("io error in {file}: {source}")]
    IoError {
        file: PathBuf,
        #[source]
        source: std::io::Error,
    },
}
```

- [ ] **Step 4: Re-run test**

```bash
cargo test -p atman-core --lib errors
```

Expected: 2 tests passed.

- [ ] **Step 5: Commit**

```bash
git add crates/atman-core/src/errors.rs
git commit -m "feat(atman-core): IngestError enum"
```

---

## Task 5: `SampleIdParser` + `DubeSampleIdParser`

**Files:**
- Modify: `crates/atman-core/src/sample_id.rs`
- Test: `crates/atman-core/tests/sample_id_dube.rs`

**Spec reference:** §5.2 (regex `^SSNA-(?P<subject>[^-]+)-(?P<condition>PR1|PR2|PT1|PT2)$`).

- [ ] **Step 1: Write the failing test**

Create `crates/atman-core/tests/sample_id_dube.rs`:

```rust
use atman_core::{DubeSampleIdParser, ParsedSampleId, SampleIdParser};

#[test]
fn parses_dube_biological_sample_ids() {
    let p = DubeSampleIdParser;
    let r = p.parse("SSNA-001B-PR1").unwrap();
    assert_eq!(r, ParsedSampleId::Biological { subject: "001B".into(), condition: "PR1".into() });

    let r = p.parse("SSNA-019-PT2").unwrap();
    assert_eq!(r, ParsedSampleId::Biological { subject: "019".into(), condition: "PT2".into() });
}

#[test]
fn classifies_control_sample_as_control() {
    let p = DubeSampleIdParser;
    let r = p.parse("CONTROL_SAMPLE_US_CS_AS_2-1").unwrap();
    assert_eq!(r, ParsedSampleId::Control);

    let r = p.parse("CONTROL_SAMPLE_US_CS_AS_2-4").unwrap();
    assert_eq!(r, ParsedSampleId::Control);
}

#[test]
fn rejects_malformed_non_control() {
    let p = DubeSampleIdParser;
    // Wrong prefix, not a control — must error.
    let r = p.parse("RANDOM-001B-PR1");
    assert!(r.is_err(), "expected parse error for non-SSNA non-control id");
}

#[test]
fn rejects_ssna_with_unknown_timepoint() {
    let p = DubeSampleIdParser;
    let r = p.parse("SSNA-001B-PX9");
    assert!(r.is_err());
}

#[test]
fn rejects_ssna_without_subject() {
    let p = DubeSampleIdParser;
    let r = p.parse("SSNA--PR1");
    assert!(r.is_err());
}
```

- [ ] **Step 2: Run to verify fail**

```bash
cargo test -p atman-core --test sample_id_dube
```

- [ ] **Step 3: Implement `crates/atman-core/src/sample_id.rs`**

```rust
//! Sample ID parsing. The Dube parser matches `^SSNA-<subject>-<PR1|PR2|PT1|PT2>$`;
//! unparseable IDs with a `CONTROL_SAMPLE_` prefix are classified as controls;
//! anything else is an error.

use crate::errors::IngestError;
use regex::Regex;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParsedSampleId {
    Biological { subject: String, condition: String },
    Control,
}

pub trait SampleIdParser {
    fn parse(&self, sample_id: &str) -> Result<ParsedSampleId, IngestError>;
}

pub struct DubeSampleIdParser;

impl SampleIdParser for DubeSampleIdParser {
    fn parse(&self, sample_id: &str) -> Result<ParsedSampleId, IngestError> {
        // Lazy-init regex once per process.
        use std::sync::OnceLock;
        static RE: OnceLock<Regex> = OnceLock::new();
        let re = RE.get_or_init(|| {
            Regex::new(r"^SSNA-(?P<subject>[^-]+)-(?P<condition>PR1|PR2|PT1|PT2)$").unwrap()
        });

        if let Some(caps) = re.captures(sample_id) {
            return Ok(ParsedSampleId::Biological {
                subject: caps["subject"].to_string(),
                condition: caps["condition"].to_string(),
            });
        }

        if sample_id.starts_with("CONTROL_SAMPLE_") {
            return Ok(ParsedSampleId::Control);
        }

        Err(IngestError::UnparseableSampleId {
            sample_id: sample_id.to_string(),
        })
    }
}
```

- [ ] **Step 4: Re-run test**

```bash
cargo test -p atman-core --test sample_id_dube
```

Expected: 5 tests passed.

- [ ] **Step 5: Commit**

```bash
git add crates/atman-core/src/sample_id.rs crates/atman-core/tests/sample_id_dube.rs
git commit -m "feat(atman-core): SampleIdParser + DubeSampleIdParser"
```

---

## Task 6: `ProteomeIngest` Trait + `IngestOutput`

**Files:**
- Modify: `crates/atman-core/src/ingest/mod.rs`

- [ ] **Step 1: Implement `crates/atman-core/src/ingest/mod.rs`**

```rust
//! Ingest boundary. Adapters implement `ProteomeIngest` to produce `IngestOutput`
//! from raw files. v0.1 ships exactly one: `olink_explore::OlinkExploreLongCsv`.

pub mod olink_explore;

use crate::{
    errors::IngestError,
    sample_id::SampleIdParser,
    types::{MeasurementRecord, Platform, ProteinIdentity, Sample},
};
use std::path::PathBuf;

/// Canonical output of any ingest adapter.
#[derive(Debug, Default)]
pub struct IngestOutput {
    pub measurements: Vec<MeasurementRecord>,
    pub proteins: Vec<ProteinIdentity>,
    pub samples: Vec<Sample>,
}

pub trait ProteomeIngest {
    fn platform(&self) -> Platform;

    /// Read one or more input files. `sample_id_parser` is consulted for every
    /// unique `sample_id`; biological samples populate `Sample::subject_id` and
    /// `condition`, controls get `is_control=true` with `subject_id=None`.
    fn read(
        &self,
        inputs: &[PathBuf],
        sample_id_parser: &dyn SampleIdParser,
    ) -> Result<IngestOutput, IngestError>;
}
```

- [ ] **Step 2: Create empty `olink_explore.rs` so the module compiles**

```rust
//! Olink Explore NGS long-CSV adapter. Implemented in Task 7.
// Intentionally empty — Task 7 populates this.
```

- [ ] **Step 3: Verify compilation**

```bash
cargo build -p atman-core
```

Expected: success.

- [ ] **Step 4: Commit**

```bash
git add crates/atman-core/src/ingest
git commit -m "feat(atman-core): ProteomeIngest trait + IngestOutput"
```

---

## Task 7: Olink Explore Long-CSV Ingest Adapter

**Files:**
- Modify: `crates/atman-core/src/ingest/olink_explore.rs`
- Test: `crates/atman-core/tests/ingest_olink.rs`
- Test: `crates/atman-core/tests/ingest_errors.rs`

This is the biggest single task. Break it into two sub-passes: (7a) happy path, (7b) error variants.

### 7a. Happy-path ingest

- [ ] **Step 1: Write the happy-path test with a synthetic 2-row fixture**

Create `crates/atman-core/tests/ingest_olink.rs`:

```rust
use atman_core::ingest::olink_explore::OlinkExploreLongCsv;
use atman_core::{
    Abundance, DubeSampleIdParser, Platform, ProteomeIngest,
};
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

fn write_fixture(dir: &TempDir, name: &str, content: &str) -> PathBuf {
    let p = dir.path().join(name);
    fs::write(&p, content).unwrap();
    p
}

// Synthetic header matches the real Dube format exactly (semicolon-delimited).
const HEADER: &str =
    "SampleID;Index;OlinkID;UniProt;Assay;MissingFreq;Panel;Panel_Lot_Nr;PlateID;QC_Warning;LOD;NPX;Normalization;Assay_Warning";

#[test]
fn ingest_minimal_happy_path() {
    let dir = TempDir::new().unwrap();
    let rows = vec![
        // SSNA-001B-PR1 × OID20838 × GGA1, passes QC
        "SSNA-001B-PR1;1;OID20838;Q9UJY5;GGA1;0.1;Neurology;B04414;PLATE1;PASS;0.20;0.15;Plate control;PASS",
        // SSNA-001B-PR2 × OID20838, WARN sample QC
        "SSNA-001B-PR2;2;OID20838;Q9UJY5;GGA1;0.1;Neurology;B04414;PLATE1;WARN;0.20;0.10;Plate control;PASS",
        // Control sample, PASS
        "CONTROL_SAMPLE_US_CS_AS_2-1;3;OID20838;Q9UJY5;GGA1;0.1;Neurology;B04414;PLATE1;PASS;0.20;0.05;Plate control;PASS",
    ];
    let content = format!("{}\n{}\n", HEADER, rows.join("\n"));
    let path = write_fixture(&dir, "synth.csv", &content);

    let adapter = OlinkExploreLongCsv;
    let parser = DubeSampleIdParser;
    let out = adapter.read(&[path], &parser).expect("ingest ok");

    assert_eq!(out.measurements.len(), 3);
    assert_eq!(out.proteins.len(), 1);
    assert_eq!(out.samples.len(), 3);
    assert_eq!(out.platform_is(Platform::OlinkExploreNgs), true);

    // First row: biological sample, PASS/PASS, abundance preserved.
    let m = &out.measurements[0];
    assert_eq!(m.sample_id, "SSNA-001B-PR1");
    assert_eq!(m.assay_id.0, "OID20838");
    assert_eq!(m.gene_symbol.as_deref(), Some("GGA1"));
    assert_eq!(m.npx_source_str, "0.15");
    assert_eq!(m.panel.as_deref(), Some("Neurology"));
    assert!(matches!(m.abundance, Abundance::Log2Npx(v) if (v - 0.15).abs() < 1e-12));
    assert!(!m.dropped_by_qc);
    assert_eq!(m.below_lod, true); // 0.15 < LOD 0.20

    // Sample sheet: controls classified.
    let ctrl_sample = out.samples.iter().find(|s| s.sample_id.starts_with("CONTROL_SAMPLE_")).unwrap();
    assert!(ctrl_sample.is_control);
    assert!(ctrl_sample.subject_id.is_none());

    let bio_sample = out.samples.iter().find(|s| s.sample_id == "SSNA-001B-PR1").unwrap();
    assert!(!bio_sample.is_control);
    assert_eq!(bio_sample.subject_id.as_deref(), Some("001B"));
    assert_eq!(bio_sample.condition.as_deref(), Some("PR1"));
}

#[test]
fn ingest_two_files_merges_samples() {
    let dir = TempDir::new().unwrap();
    let f1 = format!(
        "{}\nSSNA-001B-PR1;1;OID20838;Q9UJY5;GGA1;0.1;Neurology;B04414;PLATE1;PASS;0.20;0.15;Plate control;PASS\n",
        HEADER
    );
    let f2 = format!(
        "{}\nSSNA-001B-PR1;1;OID30000;P01375;TNF;0.1;Inflammation_II;B04414;PLATE2;PASS;0.10;0.30;Plate control;PASS\n",
        HEADER
    );
    let p1 = write_fixture(&dir, "f1.csv", &f1);
    let p2 = write_fixture(&dir, "f2.csv", &f2);

    let out = OlinkExploreLongCsv
        .read(&[p1, p2], &DubeSampleIdParser)
        .unwrap();

    assert_eq!(out.measurements.len(), 2);
    assert_eq!(out.proteins.len(), 2);
    assert_eq!(out.samples.len(), 1); // Same sample, merged.
}

// Helper extension to IngestOutput used by the first test.
trait Helpers {
    fn platform_is(&self, p: Platform) -> bool;
}
impl Helpers for atman_core::IngestOutput {
    fn platform_is(&self, p: Platform) -> bool {
        self.measurements.iter().all(|m| m.platform == p)
    }
}
```

- [ ] **Step 2: Run to verify fail**

```bash
cargo test -p atman-core --test ingest_olink
```

- [ ] **Step 3: Implement `olink_explore.rs`**

```rust
//! Olink Explore NGS long-CSV ingest adapter. See spec §5.1.

use crate::{
    errors::IngestError,
    ingest::{IngestOutput, ProteomeIngest},
    sample_id::{ParsedSampleId, SampleIdParser},
    types::{
        Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, ProteinIdentity,
        QcFlag, Sample,
    },
};
use serde::Deserialize;
use std::{
    collections::{HashMap, HashSet},
    path::{Path, PathBuf},
};

pub struct OlinkExploreLongCsv;

/// The exact 14-column header expected by this adapter, in canonical order.
const EXPECTED_HEADERS: &[&str] = &[
    "SampleID",
    "Index",
    "OlinkID",
    "UniProt",
    "Assay",
    "MissingFreq",
    "Panel",
    "Panel_Lot_Nr",
    "PlateID",
    "QC_Warning",
    "LOD",
    "NPX",
    "Normalization",
    "Assay_Warning",
];

/// Struct mapped from each row by `csv` + `serde`.
#[derive(Debug, Deserialize)]
#[serde(rename_all = "PascalCase")]
#[allow(non_snake_case)]
struct RawRow {
    #[serde(rename = "SampleID")]
    sample_id: String,
    #[serde(rename = "Index")]
    _index: String,
    #[serde(rename = "OlinkID")]
    olink_id: String,
    #[serde(rename = "UniProt")]
    uniprot: String,
    #[serde(rename = "Assay")]
    assay: String,
    #[serde(rename = "MissingFreq")]
    _missing_freq: String,
    #[serde(rename = "Panel")]
    panel: String,
    #[serde(rename = "Panel_Lot_Nr")]
    panel_lot: String,
    #[serde(rename = "PlateID")]
    plate_id: String,
    #[serde(rename = "QC_Warning")]
    qc_warning: String,
    #[serde(rename = "LOD")]
    lod: String,
    #[serde(rename = "NPX")]
    npx: String,
    #[serde(rename = "Normalization")]
    normalization: String,
    #[serde(rename = "Assay_Warning")]
    assay_warning: String,
}

impl ProteomeIngest for OlinkExploreLongCsv {
    fn platform(&self) -> Platform {
        Platform::OlinkExploreNgs
    }

    fn read(
        &self,
        inputs: &[PathBuf],
        sample_id_parser: &dyn SampleIdParser,
    ) -> Result<IngestOutput, IngestError> {
        let mut out = IngestOutput::default();
        let mut seen_pk: HashSet<(String, String)> = HashSet::new();
        let mut seen_sample: HashMap<String, u64> = HashMap::new();
        let mut protein_catalog: HashMap<String, ProteinIdentity> = HashMap::new();
        let mut ingest_counter: u64 = 0;

        for input in inputs {
            read_one_file(
                input,
                sample_id_parser,
                &mut out,
                &mut seen_pk,
                &mut seen_sample,
                &mut protein_catalog,
                &mut ingest_counter,
            )?;
        }

        // Finalize protein catalog — stable sort by (platform, assay_id).
        let mut proteins: Vec<ProteinIdentity> = protein_catalog.into_values().collect();
        proteins.sort_by(|a, b| a.assay_id.0.cmp(&b.assay_id.0));
        out.proteins = proteins;

        Ok(out)
    }
}

fn read_one_file(
    path: &Path,
    sample_id_parser: &dyn SampleIdParser,
    out: &mut IngestOutput,
    seen_pk: &mut HashSet<(String, String)>,
    seen_sample: &mut HashMap<String, u64>,
    protein_catalog: &mut HashMap<String, ProteinIdentity>,
    ingest_counter: &mut u64,
) -> Result<(), IngestError> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b';')
        .has_headers(true)
        .from_path(path)
        .map_err(|e| IngestError::CsvError {
            file: path.to_path_buf(),
            source: e,
        })?;

    // Validate header exactness.
    let headers: Vec<String> = reader
        .headers()
        .map_err(|e| IngestError::CsvError {
            file: path.to_path_buf(),
            source: e,
        })?
        .iter()
        .map(|s| s.to_string())
        .collect();
    validate_headers(path, &headers)?;

    for (line_idx, result) in reader.deserialize::<RawRow>().enumerate() {
        let raw: RawRow = result.map_err(|e| IngestError::CsvError {
            file: path.to_path_buf(),
            source: e,
        })?;
        // +2 because of 1-based line numbers + header line
        let line_number = line_idx + 2;

        // Normalization must be exactly "Plate control".
        if raw.normalization != "Plate control" {
            return Err(IngestError::UnexpectedNormalization {
                file: path.to_path_buf(),
                line: line_number,
                value: raw.normalization,
            });
        }

        // Parse numeric columns.
        let npx: f64 = raw.npx.parse().map_err(|_| IngestError::InvalidAbundance {
            file: path.to_path_buf(),
            line: line_number,
            raw: raw.npx.clone(),
        })?;
        let lod: f64 = raw.lod.parse().map_err(|_| IngestError::InvalidAbundance {
            file: path.to_path_buf(),
            line: line_number,
            raw: raw.lod.clone(),
        })?;

        // Build protein catalog entry (deduplicated across rows).
        protein_catalog
            .entry(raw.olink_id.clone())
            .or_insert_with(|| ProteinIdentity {
                platform: Platform::OlinkExploreNgs,
                assay_id: AssayId(raw.olink_id.clone()),
                uniprot: vec![raw.uniprot.clone()],
                gene_symbol: Some(raw.assay.clone()),
                panel: Some(raw.panel.clone()),
                panel_lot: Some(raw.panel_lot.clone()),
            });

        // PK uniqueness: (OlinkID, SampleID) — platform is constant.
        let pk = (raw.olink_id.clone(), raw.sample_id.clone());
        if !seen_pk.insert(pk) {
            return Err(IngestError::DuplicatePrimaryKey {
                platform: Platform::OlinkExploreNgs.as_str().to_string(),
                assay_id: raw.olink_id.clone(),
                sample_id: raw.sample_id.clone(),
            });
        }

        // Sample sheet (deduplicated by sample_id).
        if !seen_sample.contains_key(&raw.sample_id) {
            let parsed = sample_id_parser.parse(&raw.sample_id)?;
            let (subject_id, condition, is_control) = match parsed {
                ParsedSampleId::Biological { subject, condition } => {
                    (Some(subject), Some(condition), false)
                }
                ParsedSampleId::Control => (None, None, true),
            };
            let order = *ingest_counter;
            seen_sample.insert(raw.sample_id.clone(), order);
            *ingest_counter += 1;
            out.samples.push(Sample {
                sample_id: raw.sample_id.clone(),
                subject_id,
                condition,
                is_control,
                sample_type: None,
                ingest_order: order,
            });
        }

        let qc_sample = parse_qc_flag(&raw.qc_warning);
        let qc_assay = parse_qc_flag(&raw.assay_warning);

        let npx_source_str = raw.npx.clone();
        let gene_symbol = Some(raw.assay.clone());
        let panel = Some(raw.panel.clone());

        out.measurements.push(MeasurementRecord {
            platform: Platform::OlinkExploreNgs,
            assay_id: AssayId(raw.olink_id),
            gene_symbol,
            sample_id: raw.sample_id,
            abundance: Abundance::Log2Npx(npx),
            abundance_raw: Abundance::Log2Npx(npx),
            npx_source_str,
            qc_sample,
            qc_assay,
            detection_limit: DetectionLimit(Some(lod)),
            below_lod: npx < lod,
            batch: Batch {
                plate: Some(raw.plate_id),
                lot: Some(raw.panel_lot),
                run: None,
            },
            dropped_by_qc: false,
            ingest_order: *ingest_counter,
            panel,
        });
        *ingest_counter += 1;
    }

    Ok(())
}

fn validate_headers(path: &Path, headers: &[String]) -> Result<(), IngestError> {
    let got: HashSet<&str> = headers.iter().map(|s| s.as_str()).collect();
    let want: HashSet<&str> = EXPECTED_HEADERS.iter().copied().collect();
    let missing: Vec<String> = want.difference(&got).map(|s| (*s).to_string()).collect();
    let extra: Vec<String> = got.difference(&want).map(|s| (*s).to_string()).collect();
    if !missing.is_empty() || !extra.is_empty() {
        return Err(IngestError::SchemaMismatch {
            file: path.to_path_buf(),
            missing,
            extra,
        });
    }
    Ok(())
}

fn parse_qc_flag(s: &str) -> QcFlag {
    match s {
        "PASS" => QcFlag::Pass,
        "WARN" => QcFlag::Warn(String::new()),
        "FAIL" => QcFlag::Fail(String::new()),
        other => QcFlag::Warn(other.to_string()),
    }
}
```

- [ ] **Step 4: Run happy-path test**

```bash
cargo test -p atman-core --test ingest_olink
```

Expected: 2 tests passed.

- [ ] **Step 5: Commit the happy path**

```bash
git add crates/atman-core/src/ingest/olink_explore.rs crates/atman-core/tests/ingest_olink.rs
git commit -m "feat(atman-core): Olink Explore long-CSV ingest (happy path)"
```

### 7b. Error variants

- [ ] **Step 6: Write error-variant tests**

Create `crates/atman-core/tests/ingest_errors.rs`:

```rust
use atman_core::ingest::olink_explore::OlinkExploreLongCsv;
use atman_core::{DubeSampleIdParser, IngestError, ProteomeIngest};
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

const HEADER: &str =
    "SampleID;Index;OlinkID;UniProt;Assay;MissingFreq;Panel;Panel_Lot_Nr;PlateID;QC_Warning;LOD;NPX;Normalization;Assay_Warning";

fn write(dir: &TempDir, name: &str, content: &str) -> PathBuf {
    let p = dir.path().join(name);
    fs::write(&p, content).unwrap();
    p
}

#[test]
fn schema_mismatch_missing_column() {
    let dir = TempDir::new().unwrap();
    let bad_header = "SampleID;Index;OlinkID;UniProt;Assay;MissingFreq;Panel;Panel_Lot_Nr;PlateID;QC_Warning;LOD;NPX;Normalization";
    let content = format!("{}\nSSNA-001B-PR1;1;OID1;P1;X;0.1;P;L;PL;PASS;0.1;0.2;Plate control\n", bad_header);
    let p = write(&dir, "f.csv", &content);
    let err = OlinkExploreLongCsv.read(&[p], &DubeSampleIdParser).unwrap_err();
    assert!(matches!(err, IngestError::SchemaMismatch { .. }));
}

#[test]
fn duplicate_primary_key_across_two_rows() {
    let dir = TempDir::new().unwrap();
    let content = format!(
        "{}\nSSNA-001B-PR1;1;OID1;P1;X;0.1;Neurology;L;PL;PASS;0.1;0.2;Plate control;PASS\nSSNA-001B-PR1;2;OID1;P1;X;0.1;Neurology;L;PL;PASS;0.1;0.3;Plate control;PASS\n",
        HEADER
    );
    let p = write(&dir, "f.csv", &content);
    let err = OlinkExploreLongCsv.read(&[p], &DubeSampleIdParser).unwrap_err();
    assert!(matches!(err, IngestError::DuplicatePrimaryKey { .. }));
}

#[test]
fn invalid_abundance() {
    let dir = TempDir::new().unwrap();
    let content = format!(
        "{}\nSSNA-001B-PR1;1;OID1;P1;X;0.1;Neurology;L;PL;PASS;0.1;notanumber;Plate control;PASS\n",
        HEADER
    );
    let p = write(&dir, "f.csv", &content);
    let err = OlinkExploreLongCsv.read(&[p], &DubeSampleIdParser).unwrap_err();
    assert!(matches!(err, IngestError::InvalidAbundance { .. }));
}

#[test]
fn unexpected_normalization() {
    let dir = TempDir::new().unwrap();
    let content = format!(
        "{}\nSSNA-001B-PR1;1;OID1;P1;X;0.1;Neurology;L;PL;PASS;0.1;0.2;Quantile;PASS\n",
        HEADER
    );
    let p = write(&dir, "f.csv", &content);
    let err = OlinkExploreLongCsv.read(&[p], &DubeSampleIdParser).unwrap_err();
    assert!(matches!(err, IngestError::UnexpectedNormalization { .. }));
}

#[test]
fn unparseable_non_control_sample_id() {
    let dir = TempDir::new().unwrap();
    let content = format!(
        "{}\nRANDOM-X-Y;1;OID1;P1;X;0.1;Neurology;L;PL;PASS;0.1;0.2;Plate control;PASS\n",
        HEADER
    );
    let p = write(&dir, "f.csv", &content);
    let err = OlinkExploreLongCsv.read(&[p], &DubeSampleIdParser).unwrap_err();
    assert!(matches!(err, IngestError::UnparseableSampleId { .. }));
}
```

- [ ] **Step 7: Run error tests**

```bash
cargo test -p atman-core --test ingest_errors
```

Expected: 5 tests passed.

- [ ] **Step 8: Commit**

```bash
git add crates/atman-core/tests/ingest_errors.rs
git commit -m "test(atman-core): ingest error variant coverage"
```

---

## Task 8: QC Rule (Dube)

**Files:**
- Modify: `crates/atman-core/src/qc.rs`
- Test: `crates/atman-core/tests/qc_dube.rs`

**Spec reference:** §7.1.

- [ ] **Step 1: Write the failing test**

Create `crates/atman-core/tests/qc_dube.rs`:

```rust
use atman_core::qc::apply_dube_rule;
use atman_core::{
    Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, QcFlag,
};

fn rec(qc_s: QcFlag, qc_a: QcFlag) -> MeasurementRecord {
    MeasurementRecord {
        platform: Platform::OlinkExploreNgs,
        assay_id: AssayId("OID1".into()),
        gene_symbol: Some("X".into()),
        sample_id: "SSNA-001B-PR1".into(),
        abundance: Abundance::Log2Npx(0.5),
        abundance_raw: Abundance::Log2Npx(0.5),
        npx_source_str: "0.5".into(),
        qc_sample: qc_s,
        qc_assay: qc_a,
        detection_limit: DetectionLimit(Some(0.2)),
        below_lod: false,
        batch: Batch { plate: None, lot: None, run: None },
        dropped_by_qc: false,
        ingest_order: 0,
        panel: Some("Neurology".into()),
    }
}

#[test]
fn pass_pass_unchanged() {
    let mut r = rec(QcFlag::Pass, QcFlag::Pass);
    apply_dube_rule(&mut r);
    assert!(!r.dropped_by_qc);
    assert!(matches!(r.abundance, Abundance::Log2Npx(v) if (v - 0.5).abs() < 1e-12));
}

#[test]
fn warn_sample_masks_abundance() {
    let mut r = rec(QcFlag::Warn("".into()), QcFlag::Pass);
    apply_dube_rule(&mut r);
    assert!(r.dropped_by_qc);
    // abundance_raw must remain intact
    assert!(matches!(r.abundance_raw, Abundance::Log2Npx(v) if (v - 0.5).abs() < 1e-12));
}

#[test]
fn warn_assay_masks_abundance() {
    let mut r = rec(QcFlag::Pass, QcFlag::Warn("".into()));
    apply_dube_rule(&mut r);
    assert!(r.dropped_by_qc);
}

#[test]
fn warn_both_masks() {
    let mut r = rec(QcFlag::Warn("".into()), QcFlag::Warn("".into()));
    apply_dube_rule(&mut r);
    assert!(r.dropped_by_qc);
}

#[test]
fn idempotent() {
    let mut r = rec(QcFlag::Warn("".into()), QcFlag::Pass);
    apply_dube_rule(&mut r);
    let snapshot = r.clone();
    apply_dube_rule(&mut r);
    assert_eq!(r, snapshot);
}
```

- [ ] **Step 2: Run to verify fail**

```bash
cargo test -p atman-core --test qc_dube
```

- [ ] **Step 3: Implement `qc.rs`**

```rust
//! QC rules. v0.1 implements only the Dube rule (§7.1): mask `abundance` to a
//! sentinel when either `qc_sample` or `qc_assay` is not `Pass`. `abundance_raw`
//! is never mutated.

use crate::types::MeasurementRecord;

/// Apply the Dube filter rule to a single record in-place. Idempotent.
pub fn apply_dube_rule(record: &mut MeasurementRecord) {
    let masked = !record.qc_sample.is_pass() || !record.qc_assay.is_pass();
    if masked {
        record.dropped_by_qc = true;
        // abundance is masked by setting dropped_by_qc=true. The on-disk long TSV
        // writer treats `dropped_by_qc=1` rows as empty-cell when serializing the
        // `abundance` column to the Dube wide format. In memory, we keep the
        // original value as well, so higher layers can reconstruct if needed.
    }
}

/// Convenience: apply rule to all records in a slice.
pub fn apply_dube_rule_all(records: &mut [MeasurementRecord]) {
    for r in records.iter_mut() {
        apply_dube_rule(r);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, QcFlag,
    };

    #[test]
    fn rule_sets_dropped_by_qc_on_warn() {
        let mut r = MeasurementRecord {
            platform: Platform::OlinkExploreNgs,
            assay_id: AssayId("X".into()),
            gene_symbol: Some("X".into()),
            sample_id: "S".into(),
            abundance: Abundance::Log2Npx(1.0),
            abundance_raw: Abundance::Log2Npx(1.0),
            npx_source_str: "1.0".into(),
            qc_sample: QcFlag::Warn("".into()),
            qc_assay: QcFlag::Pass,
            detection_limit: DetectionLimit(None),
            below_lod: false,
            batch: Batch { plate: None, lot: None, run: None },
            dropped_by_qc: false,
            ingest_order: 0,
            panel: None,
        };
        apply_dube_rule(&mut r);
        assert!(r.dropped_by_qc);
    }
}
```

- [ ] **Step 4: Run tests**

```bash
cargo test -p atman-core --test qc_dube
cargo test -p atman-core --lib qc
```

Expected: both green.

- [ ] **Step 5: Commit**

```bash
git add crates/atman-core/src/qc.rs crates/atman-core/tests/qc_dube.rs
git commit -m "feat(atman-core): dube QC rule"
```

---

## Task 9: Long → Dube-Wide Pivot

**Files:**
- Modify: `crates/atman-core/src/matrix.rs`
- Test: `crates/atman-core/tests/pivot_dube.rs`

**Spec reference:** §7.3 (post-correction version).

**Key design points driven by empirical inspection of Dube's reference files:**
1. Column headers are **gene symbols** (`ACTN4`, `HLA-DRA`, `MICB_MICA`, …), not OlinkIDs. Carry them on `MeasurementRecord.gene_symbol` (already added in Task 3) and use that as the within-panel column key.
2. Cell values are the **original NPX source strings** (`0.098`, `-0.1219`, …), not a reformatted f64. Carry them on `MeasurementRecord.npx_source_str` and write verbatim to avoid any round-trip drift.
3. Rows partition into biological (Dube-parseable sample IDs) first, then control samples in ingest order.

The pivot output therefore carries `Vec<Option<String>>` per row, not `Vec<Option<f64>>` — it's a display pivot, not a numeric one. Fold change consumes the numeric long TSV directly.

- [ ] **Step 1: Write the failing test**

Create `crates/atman-core/tests/pivot_dube.rs`:

```rust
use atman_core::matrix::{dube_wide_pivot, DubeWidePanel};
use atman_core::{
    Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, QcFlag, Sample,
};

fn m(
    sample: &str,
    assay_id: &str,
    gene: &str,
    panel: &str,
    source: &str,
    val: f64,
    dropped: bool,
    order: u64,
) -> MeasurementRecord {
    MeasurementRecord {
        platform: Platform::OlinkExploreNgs,
        assay_id: AssayId(assay_id.into()),
        gene_symbol: Some(gene.into()),
        sample_id: sample.into(),
        abundance: Abundance::Log2Npx(val),
        abundance_raw: Abundance::Log2Npx(val),
        npx_source_str: source.into(),
        qc_sample: QcFlag::Pass,
        qc_assay: QcFlag::Pass,
        detection_limit: DetectionLimit(None),
        below_lod: false,
        batch: Batch { plate: None, lot: None, run: None },
        dropped_by_qc: dropped,
        ingest_order: order,
        panel: Some(panel.into()),
    }
}

fn sample_bio(id: &str, subj: &str, cond: &str, order: u64) -> Sample {
    Sample {
        sample_id: id.into(),
        subject_id: Some(subj.into()),
        condition: Some(cond.into()),
        is_control: false,
        sample_type: None,
        ingest_order: order,
    }
}

fn sample_ctl(id: &str, order: u64) -> Sample {
    Sample {
        sample_id: id.into(),
        subject_id: None,
        condition: None,
        is_control: true,
        sample_type: None,
        ingest_order: order,
    }
}

#[test]
fn pivot_two_panels_bio_then_control_by_gene_symbol() {
    let measurements = vec![
        m("SSNA-001B-PR1", "OID2", "ZZZ", "P1", "1.0", 1.0, false, 0),
        m("SSNA-001B-PR1", "OID1", "AAA", "P1", "2.0", 2.0, false, 1),
        m("SSNA-001B-PR2", "OID1", "AAA", "P1", "3.0", 3.0, false, 2),
        m("SSNA-001B-PR2", "OID2", "ZZZ", "P1", "4.0", 4.0, false, 3),
        m("CONTROL_SAMPLE_X", "OID1", "AAA", "P1", "9.0", 9.0, false, 4),
        m("CONTROL_SAMPLE_X", "OID2", "ZZZ", "P1", "8.0", 8.0, false, 5),
        m("SSNA-001B-PR1", "OID3", "BBB", "P2", "5.0", 5.0, false, 6),
    ];
    let samples = vec![
        sample_bio("SSNA-001B-PR1", "001B", "PR1", 0),
        sample_bio("SSNA-001B-PR2", "001B", "PR2", 1),
        sample_ctl("CONTROL_SAMPLE_X", 2),
    ];

    let panels = dube_wide_pivot(&measurements, &samples);
    assert_eq!(panels.len(), 2);

    let p1: &DubeWidePanel = panels.iter().find(|p| p.panel == "P1").unwrap();
    // Columns by GENE SYMBOL, sorted
    assert_eq!(p1.assays, vec!["AAA".to_string(), "ZZZ".to_string()]);
    // Row order: (001B, PR1), (001B, PR2), control
    assert_eq!(p1.rows.len(), 3);
    assert_eq!(p1.rows[0].participant, "001B");
    assert_eq!(p1.rows[0].exposure, "PR1");
    assert_eq!(p1.rows[0].sample_id, "SSNA-001B-PR1");
    // Values are the SOURCE STRINGS, not reformatted numbers
    assert_eq!(p1.rows[0].values[0], Some("2.0".to_string())); // AAA
    assert_eq!(p1.rows[0].values[1], Some("1.0".to_string())); // ZZZ

    assert_eq!(p1.rows[1].sample_id, "SSNA-001B-PR2");
    assert_eq!(p1.rows[2].participant, "");
    assert_eq!(p1.rows[2].exposure, "");
    assert_eq!(p1.rows[2].sample_id, "CONTROL_SAMPLE_X");
}

#[test]
fn dropped_by_qc_becomes_missing() {
    let m1 = m("SSNA-001B-PR1", "OID1", "AAA", "P1", "1.0", 1.0, true, 0);
    let m2 = m("SSNA-001B-PR1", "OID2", "BBB", "P1", "2.0", 2.0, false, 1);
    let s = vec![sample_bio("SSNA-001B-PR1", "001B", "PR1", 0)];
    let panels = dube_wide_pivot(&[m1, m2], &s);
    assert_eq!(panels.len(), 1);
    let row = &panels[0].rows[0];
    assert_eq!(row.values[0], None); // AAA masked
    assert_eq!(row.values[1], Some("2.0".to_string())); // BBB kept
}

#[test]
fn control_rows_preserve_ingest_order() {
    // Two controls, second ingested later. Must appear in ingest order.
    let measurements = vec![
        m("CONTROL_B", "OID1", "AAA", "P1", "1.0", 1.0, false, 0),
        m("CONTROL_A", "OID1", "AAA", "P1", "2.0", 2.0, false, 1),
    ];
    let samples = vec![
        sample_ctl("CONTROL_B", 0),
        sample_ctl("CONTROL_A", 1),
    ];
    let panels = dube_wide_pivot(&measurements, &samples);
    assert_eq!(panels[0].rows[0].sample_id, "CONTROL_B"); // ingest order, not alpha
    assert_eq!(panels[0].rows[1].sample_id, "CONTROL_A");
}
```

- [ ] **Step 2: Run to verify fail**

```bash
cargo test -p atman-core --test pivot_dube
```

Expected: compilation error — `dube_wide_pivot` / `DubeWidePanel` / `DubeWideRow` don't exist yet.

- [ ] **Step 3: Implement `crates/atman-core/src/matrix.rs`**

```rust
//! Long → Dube-wide pivot. See spec §7.3.
//!
//! Produces one `DubeWidePanel` per raw Olink `Panel` value, keyed by gene
//! symbol (the raw `Assay` column), with string cell values copied verbatim
//! from `MeasurementRecord.npx_source_str`. This is a display pivot: no
//! numeric arithmetic, no reformatting, so the output can reproduce Dube's
//! published filtered NPX byte-for-byte.

use crate::types::{MeasurementRecord, Sample};
use std::collections::{BTreeMap, BTreeSet, HashMap};

#[derive(Debug, Clone, PartialEq)]
pub struct DubeWideRow {
    /// Empty string for control samples; subject_id for biological samples.
    pub participant: String,
    /// Empty string for control samples; condition for biological samples.
    pub exposure: String,
    pub sample_id: String,
    /// Aligned with `DubeWidePanel::assays`. `None` = missing / masked.
    pub values: Vec<Option<String>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DubeWidePanel {
    /// Raw panel name from the Olink CSV, e.g. `Inflammation_II`. The file
    /// writer lowercases this to produce `inflammation_ii_npx.csv`.
    pub panel: String,
    /// Gene symbols, sorted by Unicode codepoint order.
    pub assays: Vec<String>,
    pub rows: Vec<DubeWideRow>,
}

pub fn dube_wide_pivot(
    measurements: &[MeasurementRecord],
    samples: &[Sample],
) -> Vec<DubeWidePanel> {
    let sample_by_id: HashMap<&str, &Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    let mut by_panel: BTreeMap<String, PanelBuilder> = BTreeMap::new();

    for m in measurements {
        let panel = match &m.panel {
            Some(p) => p.clone(),
            None => continue,
        };
        let gene = match &m.gene_symbol {
            Some(g) => g.clone(),
            None => continue, // no gene symbol → cannot place in wide column
        };
        let entry = by_panel.entry(panel).or_insert_with(PanelBuilder::default);
        entry.note_assay(gene.clone());
        entry.set_value(&m.sample_id, &gene, m);
    }

    by_panel
        .into_iter()
        .map(|(panel_name, builder)| builder.finalize(panel_name, &sample_by_id))
        .collect()
}

#[derive(Default)]
struct PanelBuilder {
    assays: BTreeSet<String>, // sorted by insertion-independent codepoint order
    /// sample_id → gene_symbol → Option<source string>
    cells: HashMap<String, HashMap<String, Option<String>>>,
}

impl PanelBuilder {
    fn note_assay(&mut self, gene: String) {
        self.assays.insert(gene);
    }

    fn set_value(&mut self, sample_id: &str, gene: &str, m: &MeasurementRecord) {
        let value = if m.dropped_by_qc {
            None
        } else {
            Some(m.npx_source_str.clone())
        };
        self.cells
            .entry(sample_id.to_string())
            .or_default()
            .insert(gene.to_string(), value);
    }

    fn finalize(
        self,
        panel: String,
        sample_by_id: &HashMap<&str, &Sample>,
    ) -> DubeWidePanel {
        let assays: Vec<String> = self.assays.into_iter().collect(); // sorted

        let mut bio: Vec<DubeWideRow> = Vec::new();
        let mut ctl: Vec<(u64, DubeWideRow)> = Vec::new();

        for (sample_id, cell_map) in self.cells {
            let s = match sample_by_id.get(sample_id.as_str()) {
                Some(s) => *s,
                None => continue,
            };
            let values: Vec<Option<String>> = assays
                .iter()
                .map(|gene| cell_map.get(gene).cloned().unwrap_or(None))
                .collect();

            if s.is_control {
                let row = DubeWideRow {
                    participant: String::new(),
                    exposure: String::new(),
                    sample_id: s.sample_id.clone(),
                    values,
                };
                ctl.push((s.ingest_order, row));
            } else {
                let row = DubeWideRow {
                    participant: s.subject_id.clone().unwrap_or_default(),
                    exposure: s.condition.clone().unwrap_or_default(),
                    sample_id: s.sample_id.clone(),
                    values,
                };
                bio.push(row);
            }
        }

        bio.sort_by(|a, b| {
            a.participant
                .cmp(&b.participant)
                .then_with(|| a.exposure.cmp(&b.exposure))
        });
        ctl.sort_by_key(|(order, _)| *order);

        let mut rows = bio;
        rows.extend(ctl.into_iter().map(|(_, r)| r));

        DubeWidePanel {
            panel,
            assays,
            rows,
        }
    }
}
```

- [ ] **Step 4: Run pivot tests**

```bash
cargo test -p atman-core --test pivot_dube
cargo test -p atman-core
```

Expected: all tests green.

- [ ] **Step 5: Commit**

```bash
git add crates/atman-core/src/matrix.rs crates/atman-core/tests/pivot_dube.rs
git commit -m "feat(atman-core): long → Dube-wide pivot (gene-symbol keyed, string values)"
```

---

## Task 10: Log2 Fold Change

**Files:**
- Modify: `crates/atman-core/src/fold_change.rs`
- Test: `crates/atman-core/tests/fold_change.rs`

**Spec reference:** §7.4.

- [ ] **Step 1: Write the failing test**

Create `crates/atman-core/tests/fold_change.rs`:

```rust
use atman_core::fold_change::{compute_log2_fc, FoldChangeInput, FoldChangeOutput, Comparison};

#[test]
fn hand_computed_three_participants_two_exposures() {
    // Participants: P1, P2, P3
    // Exposures: PR1, PT1
    // Panel P, Assay A
    // PR1 values: 0.0, 1.0, 2.0  -> 2^v = 1, 2, 4 -> mean = 7/3
    // PT1 values: 1.0, 2.0, 3.0  -> 2^v = 2, 4, 8 -> mean = 14/3
    // FC(PT1-PR1) = log2((14/3) / (7/3)) = log2(2) = 1.0
    let input = FoldChangeInput::from_cells(vec![
        // (panel, assay, participant, exposure, value)
        ("P", "A", "P1", "PR1", Some(0.0)),
        ("P", "A", "P2", "PR1", Some(1.0)),
        ("P", "A", "P3", "PR1", Some(2.0)),
        ("P", "A", "P1", "PT1", Some(1.0)),
        ("P", "A", "P2", "PT1", Some(2.0)),
        ("P", "A", "P3", "PT1", Some(3.0)),
    ]);
    let comps = vec![Comparison { a: "PT1".into(), b: "PR1".into() }];
    let out: FoldChangeOutput = compute_log2_fc(&input, &comps);
    assert_eq!(out.panels.len(), 1);
    let panel = &out.panels[0];
    assert_eq!(panel.panel, "P");
    assert_eq!(panel.assays, vec!["A".to_string()]);
    let fc = panel.values[0][0];
    assert!((fc.unwrap() - 1.0).abs() < 1e-12);
}

#[test]
fn missing_in_one_group_drops_participant_not_panel() {
    // Only P1 has PT1 data; P2 has only PR1. Fold change should still compute
    // over AVAILABLE participants per group.
    let input = FoldChangeInput::from_cells(vec![
        ("P", "A", "P1", "PR1", Some(0.0)), // 2^0 = 1
        ("P", "A", "P2", "PR1", Some(2.0)), // 2^2 = 4; mean = 2.5
        ("P", "A", "P1", "PT1", Some(1.0)), // 2^1 = 2; mean = 2
        // ("P","A","P2","PT1", None) intentionally absent
    ]);
    let comps = vec![Comparison { a: "PT1".into(), b: "PR1".into() }];
    let out = compute_log2_fc(&input, &comps);
    let fc = out.panels[0].values[0][0].unwrap();
    let expected = (2.0_f64 / 2.5_f64).log2();
    assert!((fc - expected).abs() < 1e-12, "fc={}, expected={}", fc, expected);
}

#[test]
fn empty_group_produces_missing_fc() {
    let input = FoldChangeInput::from_cells(vec![
        ("P", "A", "P1", "PR1", Some(1.0)),
        // no PT1
    ]);
    let comps = vec![Comparison { a: "PT1".into(), b: "PR1".into() }];
    let out = compute_log2_fc(&input, &comps);
    assert!(out.panels[0].values[0][0].is_none());
}

#[test]
fn self_vs_self_is_zero() {
    let input = FoldChangeInput::from_cells(vec![
        ("P", "A", "P1", "PR1", Some(2.0)),
        ("P", "A", "P2", "PR1", Some(3.0)),
    ]);
    let comps = vec![Comparison { a: "PR1".into(), b: "PR1".into() }];
    let out = compute_log2_fc(&input, &comps);
    assert_eq!(out.panels[0].values[0][0], Some(0.0));
}
```

- [ ] **Step 2: Run to verify fail**

```bash
cargo test -p atman-core --test fold_change
```

- [ ] **Step 3: Implement `fold_change.rs`**

```rust
//! Log2 fold change per spec §7.4. Compute per (panel, assay, comparison):
//! 1. collect values in group A and group B, excluding missing
//! 2. if either group is empty → FC is None
//! 3. else FC = log2(mean(2^a) / mean(2^b))

use std::collections::BTreeMap;

#[derive(Debug, Clone)]
pub struct Comparison {
    pub a: String,
    pub b: String,
}

pub struct FoldChangeInput {
    /// (panel, assay) -> (exposure -> Vec<(participant, value)>)
    panels: BTreeMap<(String, String), BTreeMap<String, Vec<(String, f64)>>>,
}

impl FoldChangeInput {
    /// Construct from a sequence of cells. Accepts any iterator of tuples
    /// where the four string fields deref to `&str` — `&str` works for
    /// tests with string literals, `String` works for CLI code that owns
    /// the strings. Missing values (`None`) are dropped silently.
    pub fn from_cells<I, S1, S2, S3, S4>(cells: I) -> Self
    where
        I: IntoIterator<Item = (S1, S2, S3, S4, Option<f64>)>,
        S1: AsRef<str>,
        S2: AsRef<str>,
        S3: AsRef<str>,
        S4: AsRef<str>,
    {
        let mut panels: BTreeMap<(String, String), BTreeMap<String, Vec<(String, f64)>>> =
            BTreeMap::new();
        for (panel, assay, participant, exposure, value) in cells {
            if let Some(v) = value {
                panels
                    .entry((panel.as_ref().to_string(), assay.as_ref().to_string()))
                    .or_default()
                    .entry(exposure.as_ref().to_string())
                    .or_default()
                    .push((participant.as_ref().to_string(), v));
            }
        }
        Self { panels }
    }
}

pub struct FoldChangePanel {
    pub panel: String,
    pub assays: Vec<String>,
    /// Outer index: assay (aligned with `assays`). Inner: comparisons (aligned
    /// with the `comparisons` argument passed to `compute_log2_fc`).
    pub values: Vec<Vec<Option<f64>>>,
}

pub struct FoldChangeOutput {
    pub panels: Vec<FoldChangePanel>,
}

pub fn compute_log2_fc(
    input: &FoldChangeInput,
    comparisons: &[Comparison],
) -> FoldChangeOutput {
    // Group (panel, assay) -> assay rows within that panel.
    let mut per_panel: BTreeMap<String, Vec<(String, &BTreeMap<String, Vec<(String, f64)>>)>> =
        BTreeMap::new();
    for ((panel, assay), by_exposure) in &input.panels {
        per_panel
            .entry(panel.clone())
            .or_default()
            .push((assay.clone(), by_exposure));
    }

    let mut out_panels = Vec::new();
    for (panel_name, mut assay_rows) in per_panel {
        assay_rows.sort_by(|a, b| a.0.cmp(&b.0));
        let assays: Vec<String> = assay_rows.iter().map(|(a, _)| a.clone()).collect();
        let mut values: Vec<Vec<Option<f64>>> = Vec::with_capacity(assays.len());
        for (_assay, by_exposure) in &assay_rows {
            let mut row: Vec<Option<f64>> = Vec::with_capacity(comparisons.len());
            for c in comparisons {
                row.push(compute_one_fc(by_exposure, &c.a, &c.b));
            }
            values.push(row);
        }
        out_panels.push(FoldChangePanel {
            panel: panel_name,
            assays,
            values,
        });
    }
    FoldChangeOutput { panels: out_panels }
}

fn compute_one_fc(
    by_exposure: &BTreeMap<String, Vec<(String, f64)>>,
    a: &str,
    b: &str,
) -> Option<f64> {
    let va = by_exposure.get(a)?;
    let vb = by_exposure.get(b)?;
    if va.is_empty() || vb.is_empty() {
        return None;
    }
    let mean_a = mean_linear(va);
    let mean_b = mean_linear(vb);
    if mean_b == 0.0 {
        return None;
    }
    Some((mean_a / mean_b).log2())
}

fn mean_linear(values: &[(String, f64)]) -> f64 {
    let sum: f64 = values.iter().map(|(_, v)| 2.0_f64.powf(*v)).sum();
    sum / (values.len() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mean_linear_basic() {
        let v = vec![("a".into(), 0.0), ("b".into(), 1.0), ("c".into(), 2.0)];
        assert!((mean_linear(&v) - 7.0 / 3.0).abs() < 1e-12);
    }
}
```

- [ ] **Step 4: Run tests**

```bash
cargo test -p atman-core --test fold_change
```

Expected: 4 tests passed.

- [ ] **Step 5: Commit**

```bash
git add crates/atman-core/src/fold_change.rs crates/atman-core/tests/fold_change.rs
git commit -m "feat(atman-core): log2 fold change (linear-mean, log2-ratio)"
```

---

## Task 11: Property Tests

**Files:**
- Test: `crates/atman-core/tests/properties.rs`

**Spec reference:** §10.2.

- [ ] **Step 1: Write the property tests**

Create `crates/atman-core/tests/properties.rs`:

```rust
use atman_core::fold_change::{compute_log2_fc, Comparison, FoldChangeInput};
use atman_core::matrix::dube_wide_pivot;
use atman_core::{
    Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, QcFlag, Sample,
};
use proptest::prelude::*;

fn arb_value() -> impl Strategy<Value = f64> {
    (-5.0_f64..5.0_f64).prop_filter("finite", |v| v.is_finite())
}

// Property: fold_change(x, x) == 0 whenever both sides use the same participant set.
proptest! {
    #[test]
    fn self_vs_self_is_zero_when_symmetric(values in prop::collection::vec(arb_value(), 1..6)) {
        let cells: Vec<_> = values
            .iter()
            .enumerate()
            .map(|(i, v)| ("P", "A", format!("P{}", i).leak() as &'static str, "PR1", Some(*v)))
            .collect();
        let input = FoldChangeInput::from_cells(cells.iter().map(|(p,a,pt,ex,v)| (*p,*a,*pt,*ex,*v)));
        let comps = vec![Comparison { a: "PR1".into(), b: "PR1".into() }];
        let out = compute_log2_fc(&input, &comps);
        let fc = out.panels[0].values[0][0].unwrap();
        prop_assert!(fc.abs() < 1e-12, "self vs self was {}", fc);
    }
}

// Property: QC masking is idempotent.
#[test]
fn pivot_is_deterministic_across_input_order() {
    use atman_core::matrix::DubeWidePanel;

    let s = vec![
        Sample { sample_id: "SSNA-1-PR1".into(), subject_id: Some("1".into()),
                 condition: Some("PR1".into()), is_control: false, sample_type: None, ingest_order: 0 },
        Sample { sample_id: "SSNA-2-PR1".into(), subject_id: Some("2".into()),
                 condition: Some("PR1".into()), is_control: false, sample_type: None, ingest_order: 1 },
    ];
    let base = vec![
        mk("SSNA-1-PR1", "Z", "P", 1.0, 0),
        mk("SSNA-1-PR1", "A", "P", 2.0, 1),
        mk("SSNA-2-PR1", "A", "P", 3.0, 2),
        mk("SSNA-2-PR1", "Z", "P", 4.0, 3),
    ];
    let mut reversed = base.clone();
    reversed.reverse();
    let a = dube_wide_pivot(&base, &s);
    let b = dube_wide_pivot(&reversed, &s);
    // Drop ingest-order-sensitive fields (there are none for biological rows) —
    // compare directly.
    assert_eq!(a, b, "pivot output must not depend on input ordering");
}

fn mk(sample: &str, assay: &str, panel: &str, val: f64, order: u64) -> MeasurementRecord {
    MeasurementRecord {
        platform: Platform::OlinkExploreNgs,
        assay_id: AssayId(assay.into()),
        gene_symbol: Some(assay.into()),
        sample_id: sample.into(),
        abundance: Abundance::Log2Npx(val),
        abundance_raw: Abundance::Log2Npx(val),
        npx_source_str: format!("{}", val),
        qc_sample: QcFlag::Pass,
        qc_assay: QcFlag::Pass,
        detection_limit: DetectionLimit(None),
        below_lod: false,
        batch: Batch { plate: None, lot: None, run: None },
        dropped_by_qc: false,
        ingest_order: order,
        panel: Some(panel.into()),
    }
}
```

- [ ] **Step 2: Run**

```bash
cargo test -p atman-core --test properties
```

Expected: all properties pass (proptest will run 256 default cases per property).

- [ ] **Step 3: Commit**

```bash
git add crates/atman-core/tests/properties.rs
git commit -m "test(atman-core): property tests for pivot and fold change"
```

---

## Task 12: `atman` Binary Scaffold

**Files:**
- Create: `crates/atman/Cargo.toml`
- Create: `crates/atman/src/main.rs`
- Create: `crates/atman/src/commands/mod.rs`

- [ ] **Step 1: Create `Cargo.toml`**

```toml
[package]
name = "atman"
version.workspace = true
edition.workspace = true
license.workspace = true

[dependencies]
anyhow = { workspace = true }
clap = { workspace = true }
csv = { workspace = true }
atman-core = { path = "../atman-core" }
serde = { workspace = true }
tempfile = { workspace = true }

[[bin]]
name = "atman"
path = "src/main.rs"
```

- [ ] **Step 2: Create minimal `main.rs` with four subcommand stubs**

```rust
//! atman: engine binary for proteomics data. v0.1 ships four subcommands
//! targeting reproduction of the Dube et al. Olink Explore published outputs.

use anyhow::Result;
use clap::{Parser, Subcommand};

mod commands;
mod io;

#[derive(Parser, Debug)]
#[command(name = "atman", version, about = "Proteomics engine (Olink Explore NGS, v0.1)")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Parse raw NPX files into canonical long TSV + catalog + sample sheet.
    Ingest(commands::ingest::Args),
    /// Apply QC rule (Dube: mask rows where QC_Warning or Assay_Warning ≠ PASS).
    Qc(commands::qc::Args),
    /// Pivot QC'd long TSV into per-panel Dube-wide CSVs.
    Matrix(commands::matrix::Args),
    /// Compute log2 fold change per panel for the given comparisons.
    FoldChange(commands::fold_change::Args),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Ingest(args) => commands::ingest::run(args),
        Command::Qc(args) => commands::qc::run(args),
        Command::Matrix(args) => commands::matrix::run(args),
        Command::FoldChange(args) => commands::fold_change::run(args),
    }
}
```

- [ ] **Step 3: Create empty command stubs so the binary compiles**

`crates/atman/src/commands/mod.rs`:

```rust
pub mod fold_change;
pub mod ingest;
pub mod matrix;
pub mod qc;
```

Create four files with stub `Args` and `run`:

`crates/atman/src/commands/ingest.rs`:
```rust
use anyhow::Result;
use clap::Args as ClapArgs;

#[derive(ClapArgs, Debug)]
pub struct Args {}
pub fn run(_args: Args) -> Result<()> {
    anyhow::bail!("ingest: not yet implemented (Task 14)")
}
```

Repeat analogous stubs for `qc.rs`, `matrix.rs`, `fold_change.rs`.

`crates/atman/src/io.rs`:
```rust
// Populated in Task 13.
```

- [ ] **Step 4: Verify compilation**

```bash
cargo build -p atman
```

Expected: success with unused-code warnings.

- [ ] **Step 5: Commit**

```bash
git add crates/atman
git commit -m "feat(atman): binary scaffold with four subcommand stubs"
```

---

## Task 13: `io.rs` — Atomic Writer + Long TSV Read/Write

**Files:**
- Modify: `crates/atman/src/io.rs`

- [ ] **Step 1: Write the failing test (inline)**

Add to `io.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn atomic_write_creates_file_with_content() {
        let d = TempDir::new().unwrap();
        let p = d.path().join("out.tsv");
        atomic_write(&p, b"hello\n").unwrap();
        assert_eq!(std::fs::read_to_string(&p).unwrap(), "hello\n");
    }

    #[test]
    fn atomic_write_overwrites_existing() {
        let d = TempDir::new().unwrap();
        let p = d.path().join("out.tsv");
        atomic_write(&p, b"old").unwrap();
        atomic_write(&p, b"new").unwrap();
        assert_eq!(std::fs::read(&p).unwrap(), b"new");
    }
}
```

- [ ] **Step 2: Implement `io.rs`**

```rust
//! File I/O helpers: atomic writes, long-TSV read/write, CSV write.
//!
//! The long TSV is header-indexed on read (not positional) so adding columns
//! in the future does not silently mis-read older files. The writer emits a
//! stable column order but the reader tolerates either ordering.

use anyhow::{anyhow, Context, Result};
use atman_core::{
    fold_change::FoldChangePanel,
    matrix::DubeWidePanel,
    Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, Platform, ProteinIdentity,
    QcFlag, Sample,
};
use std::{
    collections::HashMap,
    io::Write,
    path::{Path, PathBuf},
};

/// Atomic file write: write to a sibling tempfile, fsync, then rename over target.
pub fn atomic_write(path: &Path, bytes: &[u8]) -> Result<()> {
    let parent = path.parent().unwrap_or(Path::new("."));
    let mut tmp = tempfile::NamedTempFile::new_in(parent)
        .with_context(|| format!("creating tempfile in {:?}", parent))?;
    tmp.write_all(bytes)
        .with_context(|| format!("writing tempfile for {:?}", path))?;
    tmp.persist(path)
        .with_context(|| format!("atomic rename to {:?}", path))?;
    Ok(())
}

const LONG_TSV_HEADER: &[&str] = &[
    "sample_id",
    "assay_id",
    "gene_symbol",
    "panel",
    "npx_source_str",
    "abundance",
    "abundance_raw",
    "abundance_unit",
    "qc_sample",
    "qc_assay",
    "detection_limit",
    "below_lod",
    "dropped_by_qc",
    "plate_id",
    "panel_lot",
    "ingest_order",
];

/// Write the canonical long-format measurements TSV (spec §4.3).
/// The `abundance` column is empty for rows with `dropped_by_qc=true`;
/// `abundance_raw` and `npx_source_str` are always written (never masked).
pub fn write_measurements_long(path: &Path, records: &[MeasurementRecord]) -> Result<()> {
    let mut buf = String::new();
    buf.push_str(&LONG_TSV_HEADER.join("\t"));
    buf.push('\n');
    for r in records {
        let abundance = if r.dropped_by_qc {
            String::new()
        } else {
            format_f64(r.abundance.as_f64())
        };
        let abundance_raw = format_f64(r.abundance_raw.as_f64());
        let lod = match r.detection_limit.0 {
            Some(v) => format_f64(v),
            None => String::new(),
        };
        let plate = r.batch.plate.clone().unwrap_or_default();
        let panel = r.panel.clone().unwrap_or_default();
        let lot = r.batch.lot.clone().unwrap_or_default();
        let gene = r.gene_symbol.clone().unwrap_or_default();
        // Column order must match LONG_TSV_HEADER exactly.
        buf.push_str(&r.sample_id);             buf.push('\t');
        buf.push_str(&r.assay_id.0);            buf.push('\t');
        buf.push_str(&gene);                    buf.push('\t');
        buf.push_str(&panel);                   buf.push('\t');
        buf.push_str(&r.npx_source_str);        buf.push('\t');
        buf.push_str(&abundance);               buf.push('\t');
        buf.push_str(&abundance_raw);           buf.push('\t');
        buf.push_str("log2_npx");               buf.push('\t');
        buf.push_str(r.qc_sample.as_str());     buf.push('\t');
        buf.push_str(r.qc_assay.as_str());      buf.push('\t');
        buf.push_str(&lod);                     buf.push('\t');
        buf.push_str(&(r.below_lod as u8).to_string());          buf.push('\t');
        buf.push_str(&(r.dropped_by_qc as u8).to_string());      buf.push('\t');
        buf.push_str(&plate);                   buf.push('\t');
        buf.push_str(&lot);                     buf.push('\t');
        buf.push_str(&r.ingest_order.to_string());
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

/// Read the long-format measurements TSV, indexing columns by header name.
pub fn read_measurements_long(path: &Path) -> Result<Vec<MeasurementRecord>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;

    // Map header name -> column index.
    let headers = reader
        .headers()
        .with_context(|| format!("reading headers from {:?}", path))?
        .clone();
    let col: HashMap<String, usize> = headers
        .iter()
        .enumerate()
        .map(|(i, h)| (h.to_string(), i))
        .collect();
    let need = |name: &str| -> Result<usize> {
        col.get(name)
            .copied()
            .ok_or_else(|| anyhow!("missing column {:?} in {:?}", name, path))
    };
    let c_sample = need("sample_id")?;
    let c_assay  = need("assay_id")?;
    let c_gene   = need("gene_symbol")?;
    let c_panel  = need("panel")?;
    let c_src    = need("npx_source_str")?;
    let c_abund  = need("abundance")?;
    let c_abund_raw = need("abundance_raw")?;
    let c_qcs    = need("qc_sample")?;
    let c_qca    = need("qc_assay")?;
    let c_lod    = need("detection_limit")?;
    let c_below  = need("below_lod")?;
    let c_drop   = need("dropped_by_qc")?;
    let c_plate  = need("plate_id")?;
    let c_lot    = need("panel_lot")?;
    let c_order  = need("ingest_order")?;

    let mut out = Vec::new();
    for result in reader.records() {
        let row = result.with_context(|| format!("reading record from {:?}", path))?;
        let abundance_raw: f64 = row[c_abund_raw]
            .parse()
            .context("abundance_raw parse")?;
        let dropped_by_qc: bool = row[c_drop].parse::<u8>().map(|v| v != 0).unwrap_or(false);
        let abund_str = &row[c_abund];
        let abundance = if abund_str.is_empty() {
            Abundance::Log2Npx(abundance_raw)
        } else {
            Abundance::Log2Npx(abund_str.parse().context("abundance parse")?)
        };
        let detection_limit = if row[c_lod].is_empty() {
            DetectionLimit(None)
        } else {
            DetectionLimit(Some(row[c_lod].parse().context("detection_limit parse")?))
        };
        let below_lod: bool = row[c_below].parse::<u8>().map(|v| v != 0).unwrap_or(false);
        out.push(MeasurementRecord {
            platform: Platform::OlinkExploreNgs,
            assay_id: AssayId(row[c_assay].to_string()),
            gene_symbol: empty_to_none(&row[c_gene]),
            sample_id: row[c_sample].to_string(),
            abundance,
            abundance_raw: Abundance::Log2Npx(abundance_raw),
            npx_source_str: row[c_src].to_string(),
            qc_sample: parse_qc(&row[c_qcs]),
            qc_assay: parse_qc(&row[c_qca]),
            detection_limit,
            below_lod,
            batch: Batch {
                plate: empty_to_none(&row[c_plate]),
                lot: empty_to_none(&row[c_lot]),
                run: None,
            },
            dropped_by_qc,
            ingest_order: row[c_order].parse().unwrap_or(0),
            panel: empty_to_none(&row[c_panel]),
        });
    }
    Ok(out)
}

/// samples.tsv writer (spec §4.3).
pub fn write_samples(path: &Path, samples: &[Sample]) -> Result<()> {
    let mut buf = String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for s in samples {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            s.sample_id,
            s.subject_id.clone().unwrap_or_default(),
            s.condition.clone().unwrap_or_default(),
            s.is_control as u8,
            s.sample_type.clone().unwrap_or_default(),
            s.ingest_order,
        ));
    }
    atomic_write(path, buf.as_bytes())
}

pub fn read_samples(path: &Path) -> Result<Vec<Sample>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)?;
    let mut out = Vec::new();
    for row in reader.records() {
        let row = row?;
        out.push(Sample {
            sample_id: row[0].to_string(),
            subject_id: empty_to_none(&row[1]),
            condition: empty_to_none(&row[2]),
            is_control: row[3].parse::<u8>()? != 0,
            sample_type: empty_to_none(&row[4]),
            ingest_order: row[5].parse().unwrap_or(0),
        });
    }
    Ok(out)
}

/// proteins.tsv writer (spec §4.3).
pub fn write_proteins(path: &Path, proteins: &[ProteinIdentity]) -> Result<()> {
    let mut buf = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for p in proteins {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            p.platform.as_str(),
            p.assay_id.0,
            p.uniprot.join(","),
            p.gene_symbol.clone().unwrap_or_default(),
            p.panel.clone().unwrap_or_default(),
            p.panel_lot.clone().unwrap_or_default(),
        ));
    }
    atomic_write(path, buf.as_bytes())
}

/// Dube-wide panel CSV writer. Writes source strings verbatim — no f64
/// round-trip — so reproduction against Dube's published files is byte-exact.
pub fn write_dube_wide_panel(
    out_dir: &Path,
    panel: &DubeWidePanel,
) -> Result<PathBuf> {
    let filename = format!("{}_npx.csv", panel.panel.to_ascii_lowercase());
    let path = out_dir.join(filename);
    let mut buf = String::new();
    buf.push_str("Participant,Exposure,SampleID");
    for a in &panel.assays {
        buf.push(',');
        buf.push_str(a);
    }
    buf.push('\n');
    for row in &panel.rows {
        buf.push_str(&row.participant);
        buf.push(',');
        buf.push_str(&row.exposure);
        buf.push(',');
        buf.push_str(&row.sample_id);
        for v in &row.values {
            buf.push(',');
            if let Some(s) = v {
                buf.push_str(s); // verbatim source string
            }
        }
        buf.push('\n');
    }
    atomic_write(&path, buf.as_bytes())?;
    Ok(path)
}

/// Fold-change panel CSV writer: one file per panel, format `Assay,<comps...>`.
pub fn write_fold_change_panel(
    out_dir: &Path,
    panel: &FoldChangePanel,
    comparisons: &[atman_core::fold_change::Comparison],
) -> Result<PathBuf> {
    let filename = format!("{}_log2_fc.csv", panel.panel.to_ascii_lowercase());
    let path = out_dir.join(filename);
    let mut buf = String::from("Assay");
    for c in comparisons {
        buf.push(',');
        buf.push_str(&format!("{}-{}", c.a, c.b));
    }
    buf.push('\n');
    for (i, assay) in panel.assays.iter().enumerate() {
        buf.push_str(assay);
        for v in &panel.values[i] {
            buf.push(',');
            if let Some(x) = v {
                buf.push_str(&format_f64_fc(*x));
            }
        }
        buf.push('\n');
    }
    atomic_write(&path, buf.as_bytes())?;
    Ok(path)
}

fn parse_qc(s: &str) -> QcFlag {
    match s {
        "PASS" => QcFlag::Pass,
        "WARN" => QcFlag::Warn(String::new()),
        "FAIL" => QcFlag::Fail(String::new()),
        _ => QcFlag::Pass,
    }
}

fn empty_to_none(s: &str) -> Option<String> {
    if s.is_empty() { None } else { Some(s.to_string()) }
}

/// Format an f64 for the canonical long TSV (17 significant digits, no locale).
fn format_f64(v: f64) -> String {
    format!("{}", v)
}

/// Format an f64 for the fold-change CSV. Rust's default `{}` for f64 produces
/// the "shortest roundtrip" representation, which reproduces the input f64
/// exactly on parse. Dube's fold-change file uses R's default full-precision
/// output (~15-17 digits). The reproduction test uses a 1e-4 tolerance via
/// parse+delta, so format drift here does not affect correctness.
fn format_f64_fc(v: f64) -> String {
    format!("{}", v)
}
```

- [ ] **Step 3: Run tests**

```bash
cargo test -p atman --lib io
```

Expected: 2 tests pass.

- [ ] **Step 4: Commit**

```bash
git add crates/atman/src/io.rs
git commit -m "feat(atman): io module (atomic write, long TSV, Dube-wide, FC)"
```

---

## Task 14: `ingest` Command

**Files:**
- Modify: `crates/atman/src/commands/ingest.rs`

- [ ] **Step 1: Implement the command**

```rust
use anyhow::{Context, Result};
use clap::Args as ClapArgs;
use atman_core::{
    ingest::olink_explore::OlinkExploreLongCsv, DubeSampleIdParser, ProteomeIngest,
};
use std::path::PathBuf;

use crate::io::{write_measurements_long, write_proteins, write_samples};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Proteomics platform. v0.1 supports only `olink-explore-ngs`.
    #[arg(long, default_value = "olink-explore-ngs")]
    platform: String,

    /// Sample ID parser. v0.1 supports only `dube`.
    #[arg(long, default_value = "dube")]
    parser: String,

    /// Output directory for measurements.tsv, proteins.tsv, samples.tsv.
    #[arg(long)]
    output_dir: PathBuf,

    /// Input NPX files (one or more).
    inputs: Vec<PathBuf>,
}

pub fn run(args: Args) -> Result<()> {
    if args.platform != "olink-explore-ngs" {
        anyhow::bail!("platform {:?} not supported in v0.1", args.platform);
    }
    if args.parser != "dube" {
        anyhow::bail!("parser {:?} not supported in v0.1", args.parser);
    }
    if args.inputs.is_empty() {
        anyhow::bail!("no input files given");
    }
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let adapter = OlinkExploreLongCsv;
    let parser = DubeSampleIdParser;
    let out = adapter.read(&args.inputs, &parser)?;

    write_measurements_long(&args.output_dir.join("measurements.tsv"), &out.measurements)?;
    write_proteins(&args.output_dir.join("proteins.tsv"), &out.proteins)?;
    write_samples(&args.output_dir.join("samples.tsv"), &out.samples)?;

    eprintln!(
        "ingest: {} rows, {} assays, {} samples, {} inputs",
        out.measurements.len(),
        out.proteins.len(),
        out.samples.len(),
        args.inputs.len()
    );
    Ok(())
}
```

- [ ] **Step 2: Smoke test against the real data**

```bash
cargo run -p atman -- ingest \
    --output-dir /tmp/kpp-ingest \
    example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv \
    example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv
```

Expected: ingest succeeds (exit 0), stderr prints a one-line summary of the form `ingest: <N> rows, <M> assays, <K> samples, 2 inputs`, and `/tmp/kpp-ingest/` contains `measurements.tsv`, `proteins.tsv`, `samples.tsv`.

- [ ] **Step 3: Sanity-check the outputs (no hardcoded numbers)**

```bash
head -3 /tmp/kpp-ingest/measurements.tsv        # verify header contains gene_symbol + npx_source_str columns
wc -l /tmp/kpp-ingest/*.tsv                     # record counts for your reference
head -5 /tmp/kpp-ingest/proteins.tsv            # verify OID / UniProt / gene / panel populated
head -5 /tmp/kpp-ingest/samples.tsv             # verify control samples have is_control=1 and biological samples have subject_id + condition populated
```

Expected patterns (the exact counts will vary — write them down, they become the ground-truth for future regressions):
- `measurements.tsv` row count equals total rows in the two raw CSVs plus one header line.
- `proteins.tsv` row count matches unique `OlinkID` count across both files (around 2,943 per the Dube paper abstract).
- `samples.tsv` contains both biological samples (SSNA-… with subject_id populated) and control samples (CONTROL_SAMPLE_… with is_control=1 and no subject_id).

- [ ] **Step 4: Commit**

```bash
git add crates/atman/src/commands/ingest.rs
git commit -m "feat(atman): ingest command"
```

---

## Task 15: `qc` Command

**Files:**
- Modify: `crates/atman/src/commands/qc.rs`

**Invariant to remember:** `apply_dube_rule` sets `dropped_by_qc = true` on masked rows but does NOT overwrite `abundance` in memory. The long-TSV writer emits an empty `abundance` cell for masked rows, so on re-read the abundance field is repopulated from `abundance_raw`. Every downstream consumer that needs the *effective* value must call `MeasurementRecord::effective_abundance()` (which honors the flag) — never read `.abundance` directly. The pivot and fold-change paths already do this; new code must too.

- [ ] **Step 1: Implement**

```rust
use anyhow::{Context, Result};
use clap::Args as ClapArgs;
use atman_core::qc::apply_dube_rule_all;
use std::path::PathBuf;

use crate::io::{read_measurements_long, write_measurements_long};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing `measurements.tsv` from `ingest`.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output directory for `qc_measurements.tsv`.
    #[arg(long)]
    output_dir: PathBuf,

    /// QC rule name. v0.1 supports only `dube`.
    #[arg(long, default_value = "dube")]
    rule: String,
}

pub fn run(args: Args) -> Result<()> {
    if args.rule != "dube" {
        anyhow::bail!("rule {:?} not supported in v0.1", args.rule);
    }
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let in_path = args.input_dir.join("measurements.tsv");
    let mut records = read_measurements_long(&in_path)?;
    let before = records.iter().filter(|r| !r.dropped_by_qc).count();
    apply_dube_rule_all(&mut records);
    let after = records.iter().filter(|r| !r.dropped_by_qc).count();
    let masked = before - after;

    write_measurements_long(
        &args.output_dir.join("qc_measurements.tsv"),
        &records,
    )?;
    eprintln!(
        "qc: rule={} total={} masked={} passed={}",
        args.rule,
        records.len(),
        masked,
        after,
    );
    Ok(())
}
```

- [ ] **Step 2: Smoke test**

```bash
cargo run -p atman -- qc \
    --input-dir /tmp/kpp-ingest \
    --output-dir /tmp/kpp-qc
```

Expected: stderr prints a count. Output file exists.

- [ ] **Step 3: Commit**

```bash
git add crates/atman/src/commands/qc.rs
git commit -m "feat(atman): qc command (dube rule)"
```

---

## Task 16: `matrix` Command

**Files:**
- Modify: `crates/atman/src/commands/matrix.rs`

- [ ] **Step 1: Implement**

```rust
use anyhow::{Context, Result};
use clap::Args as ClapArgs;
use atman_core::matrix::dube_wide_pivot;
use std::path::PathBuf;

use crate::io::{read_measurements_long, read_samples, write_dube_wide_panel};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing `qc_measurements.tsv` and `samples.tsv`.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output directory for per-panel CSVs.
    #[arg(long)]
    output_dir: PathBuf,

    /// Output format. v0.1 supports only `dube-wide`.
    #[arg(long, default_value = "dube-wide")]
    format: String,

    /// Split dimension. v0.1 supports only `panel`.
    #[arg(long, default_value = "panel")]
    split_by: String,
}

pub fn run(args: Args) -> Result<()> {
    if args.format != "dube-wide" {
        anyhow::bail!("format {:?} not supported in v0.1", args.format);
    }
    if args.split_by != "panel" {
        anyhow::bail!("split-by {:?} not supported in v0.1", args.split_by);
    }
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let measurements = read_measurements_long(&args.input_dir.join("qc_measurements.tsv"))?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;

    let panels = dube_wide_pivot(&measurements, &samples);
    for panel in &panels {
        let p = write_dube_wide_panel(&args.output_dir, panel)?;
        eprintln!("matrix: wrote {:?} ({} assays, {} rows)", p, panel.assays.len(), panel.rows.len());
    }
    Ok(())
}
```

- [ ] **Step 2: Smoke test**

```bash
cargo run -p atman -- matrix \
    --input-dir /tmp/kpp-qc \
    --output-dir /tmp/kpp-matrix
# Note: matrix reads samples.tsv from --input-dir, so copy or symlink it there first:
cp /tmp/kpp-ingest/samples.tsv /tmp/kpp-qc/samples.tsv
cargo run -p atman -- matrix \
    --input-dir /tmp/kpp-qc \
    --output-dir /tmp/kpp-matrix
ls /tmp/kpp-matrix/
```

Expected: 8 `<panel>_npx.csv` files, filenames matching `cardiometabolic_npx.csv`, `cardiometabolic_ii_npx.csv`, etc.

- [ ] **Step 3: Commit**

```bash
git add crates/atman/src/commands/matrix.rs
git commit -m "feat(atman): matrix command (dube-wide pivot)"
```

---

## Task 17: `fold-change` Command

**Files:**
- Modify: `crates/atman/src/commands/fold_change.rs`

- [ ] **Step 1: Implement**

```rust
use anyhow::{Context, Result};
use clap::Args as ClapArgs;
use atman_core::fold_change::{compute_log2_fc, Comparison, FoldChangeInput};
use std::path::PathBuf;

use crate::io::{read_measurements_long, read_samples, write_fold_change_panel};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[arg(long)]
    input_dir: PathBuf,

    #[arg(long)]
    output_dir: PathBuf,

    /// Comma-separated comparisons in `A-B` form. Example: "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1".
    #[arg(long)]
    groups: String,

    /// Column (on the sample sheet) that names the exposure for each sample.
    #[arg(long, default_value = "exposure")]
    group_column: String,

    /// Column on the sample sheet that names the participant.
    #[arg(long, default_value = "participant")]
    participant_column: String,

    /// Split dimension. v0.1 supports only `panel`.
    #[arg(long, default_value = "panel")]
    split_by: String,
}

pub fn run(args: Args) -> Result<()> {
    if args.split_by != "panel" {
        anyhow::bail!("split-by {:?} not supported in v0.1", args.split_by);
    }
    if args.group_column != "exposure" || args.participant_column != "participant" {
        anyhow::bail!(
            "v0.1 expects group_column=exposure participant_column=participant (Dube schema)"
        );
    }
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let comparisons: Vec<Comparison> = args
        .groups
        .split(',')
        .map(|s| {
            let mut it = s.splitn(2, '-');
            let a = it.next().unwrap_or("").trim().to_string();
            let b = it.next().unwrap_or("").trim().to_string();
            Comparison { a, b }
        })
        .collect();
    if comparisons.is_empty() {
        anyhow::bail!("no comparisons given");
    }

    let measurements = read_measurements_long(&args.input_dir.join("qc_measurements.tsv"))?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;

    // Build FoldChangeInput. Only biological samples contribute.
    use std::collections::HashMap;
    let sample_by_id: HashMap<&str, &atman_core::Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    // Fold change uses GENE SYMBOL (not OlinkID) as the assay key so the
    // output's Assay column matches Dube's reference (ACTN4, HLA-DRA, …).
    // Rows without a gene symbol or mapped sample metadata are skipped.
    let cells: Vec<(String, String, String, String, Option<f64>)> = measurements
        .iter()
        .filter_map(|m| {
            let abundance = m.effective_abundance()?; // honors dropped_by_qc
            let s = sample_by_id.get(m.sample_id.as_str())?;
            if s.is_control {
                return None;
            }
            let participant = s.subject_id.clone()?;
            let exposure = s.condition.clone()?;
            let panel = m.panel.clone()?;
            let gene = m.gene_symbol.clone()?;
            Some((panel, gene, participant, exposure, Some(abundance)))
        })
        .collect();
    let input = FoldChangeInput::from_cells(cells);

    let output = compute_log2_fc(&input, &comparisons);
    for panel in &output.panels {
        let p = write_fold_change_panel(&args.output_dir, panel, &comparisons)?;
        eprintln!("fold-change: wrote {:?} ({} assays)", p, panel.assays.len());
    }
    Ok(())
}
```

- [ ] **Step 2: Smoke test**

```bash
cargo run -p atman -- fold-change \
    --input-dir /tmp/kpp-qc \
    --output-dir /tmp/kpp-fc \
    --groups "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1"
ls /tmp/kpp-fc/
```

Expected: 8 `<panel>_log2_fc.csv` files.

- [ ] **Step 3: Commit**

```bash
git add crates/atman/src/commands/fold_change.rs
git commit -m "feat(atman): fold-change command"
```

---

## Task 18: CSV Diff Harness

**Files:**
- Create: `crates/atman/tests/common/mod.rs`

- [ ] **Step 1: Create the harness**

```rust
//! CSV diff helper for reproduction tests. Aligns two CSVs by row key
//! (the first column) and reports the first N mismatches.

use std::collections::HashMap;
use std::path::Path;

pub struct DiffReport {
    pub file_left: String,
    pub file_right: String,
    pub header_mismatch: Option<(Vec<String>, Vec<String>)>,
    pub missing_rows: Vec<String>,
    pub extra_rows: Vec<String>,
    pub cell_mismatches: Vec<CellMismatch>,
}

pub struct CellMismatch {
    pub row_key: String,
    pub column: String,
    pub left: String,
    pub right: String,
}

impl DiffReport {
    pub fn is_clean(&self) -> bool {
        self.header_mismatch.is_none()
            && self.missing_rows.is_empty()
            && self.extra_rows.is_empty()
            && self.cell_mismatches.is_empty()
    }

    pub fn assert_clean(&self) {
        if self.is_clean() {
            return;
        }
        let mut msg = format!("\n=== DIFF FAIL: {} vs {} ===\n", self.file_left, self.file_right);
        if let Some((l, r)) = &self.header_mismatch {
            msg.push_str(&format!(
                "header mismatch:\n  left : {:?}\n  right: {:?}\n",
                l, r
            ));
        }
        if !self.missing_rows.is_empty() {
            msg.push_str(&format!(
                "missing rows (first 5): {:?}\n",
                &self.missing_rows[..self.missing_rows.len().min(5)]
            ));
        }
        if !self.extra_rows.is_empty() {
            msg.push_str(&format!(
                "extra rows (first 5): {:?}\n",
                &self.extra_rows[..self.extra_rows.len().min(5)]
            ));
        }
        if !self.cell_mismatches.is_empty() {
            msg.push_str(&format!(
                "cell mismatches (first 10 of {}):\n",
                self.cell_mismatches.len()
            ));
            for m in self.cell_mismatches.iter().take(10) {
                msg.push_str(&format!(
                    "  row={} col={} left={:?} right={:?}\n",
                    m.row_key, m.column, m.left, m.right
                ));
            }
        }
        panic!("{}", msg);
    }
}

/// Load a CSV as (header, rows_by_key). Key = first column. Remaining columns are values.
pub fn load_csv(path: &Path) -> (Vec<String>, HashMap<String, Vec<String>>) {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)
        .unwrap_or_else(|e| panic!("opening {:?}: {}", path, e));
    let header: Vec<String> = reader.headers().unwrap().iter().map(|s| s.to_string()).collect();
    let mut rows: HashMap<String, Vec<String>> = HashMap::new();
    for r in reader.records() {
        let r = r.unwrap();
        let cells: Vec<String> = r.iter().map(|s| s.to_string()).collect();
        let key = cells[0].clone();
        rows.insert(key, cells);
    }
    (header, rows)
}

/// Exact string-level diff. Use for filtered NPX files.
pub fn diff_strict(left: &Path, right: &Path) -> DiffReport {
    let (lh, lr) = load_csv(left);
    let (rh, rr) = load_csv(right);
    let mut report = DiffReport {
        file_left: left.display().to_string(),
        file_right: right.display().to_string(),
        header_mismatch: None,
        missing_rows: vec![],
        extra_rows: vec![],
        cell_mismatches: vec![],
    };
    if lh != rh {
        report.header_mismatch = Some((lh.clone(), rh.clone()));
        return report;
    }
    for (key, lrow) in &lr {
        match rr.get(key) {
            Some(rrow) => {
                for (i, col) in lh.iter().enumerate() {
                    let lv = lrow.get(i).cloned().unwrap_or_default();
                    let rv = rrow.get(i).cloned().unwrap_or_default();
                    if lv != rv {
                        report.cell_mismatches.push(CellMismatch {
                            row_key: key.clone(),
                            column: col.clone(),
                            left: lv,
                            right: rv,
                        });
                    }
                }
            }
            None => report.extra_rows.push(key.clone()),
        }
    }
    for key in rr.keys() {
        if !lr.contains_key(key) {
            report.missing_rows.push(key.clone());
        }
    }
    report
}

/// Numeric diff: exact match on the Assay column, |delta| <= tol on the rest.
/// Used for fold-change files.
pub fn diff_numeric(left: &Path, right: &Path, tol: f64) -> DiffReport {
    let (lh, lr) = load_csv(left);
    let (rh, rr) = load_csv(right);
    let mut report = DiffReport {
        file_left: left.display().to_string(),
        file_right: right.display().to_string(),
        header_mismatch: None,
        missing_rows: vec![],
        extra_rows: vec![],
        cell_mismatches: vec![],
    };
    if lh != rh {
        report.header_mismatch = Some((lh.clone(), rh.clone()));
        return report;
    }
    for (key, lrow) in &lr {
        match rr.get(key) {
            Some(rrow) => {
                for (i, col) in lh.iter().enumerate() {
                    let lv = lrow.get(i).cloned().unwrap_or_default();
                    let rv = rrow.get(i).cloned().unwrap_or_default();
                    if i == 0 {
                        if lv != rv {
                            report.cell_mismatches.push(CellMismatch {
                                row_key: key.clone(),
                                column: col.clone(),
                                left: lv,
                                right: rv,
                            });
                        }
                    } else {
                        let lf = lv.parse::<f64>().ok();
                        let rf = rv.parse::<f64>().ok();
                        match (lf, rf) {
                            (Some(a), Some(b)) if (a - b).abs() <= tol => {}
                            (None, None) if lv == rv => {}
                            _ => report.cell_mismatches.push(CellMismatch {
                                row_key: key.clone(),
                                column: col.clone(),
                                left: lv,
                                right: rv,
                            }),
                        }
                    }
                }
            }
            None => report.extra_rows.push(key.clone()),
        }
    }
    for key in rr.keys() {
        if !lr.contains_key(key) {
            report.missing_rows.push(key.clone());
        }
    }
    report
}
```

- [ ] **Step 2: Commit**

```bash
git add crates/atman/tests/common/mod.rs
git commit -m "test(atman): CSV diff harness for reproduction test"
```

---

## Task 19: Reproduction Integration Test

**Files:**
- Create: `crates/atman/tests/dube_reproduction.rs`

This is the acceptance test. It runs the full pipeline against the real data and diffs every output against Dube's published files.

- [ ] **Step 1: Write the integration test**

```rust
//! End-to-end reproduction of Dube et al. Scientific Data 2023 Olink Explore
//! published filtered NPX and fold change files. Passing this test is the
//! definition of v0.1 done.

mod common;

use common::{diff_numeric, diff_strict};
use std::path::PathBuf;
use std::process::Command;

fn repo_root() -> PathBuf {
    // CARGO_MANIFEST_DIR = crates/atman
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

fn run_atman(args: &[&str]) {
    let bin = env!("CARGO_BIN_EXE_atman");
    let status = Command::new(bin)
        .args(args)
        .status()
        .expect("run atman");
    assert!(status.success(), "atman {:?} failed", args);
}

const PANELS: &[&str] = &[
    "cardiometabolic",
    "cardiometabolic_ii",
    "inflammation",
    "inflammation_ii",
    "neurology",
    "neurology_ii",
    "oncology",
    "oncology_ii",
];

#[test]
fn dube_reproduction_end_to_end() {
    let root = repo_root();
    let data = root.join("example_data/dube_heat_2023");
    let raw1 = data.join("20212016_Dube_NPX_2021-11-30.csv");
    let raw2 = data.join("20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv");
    assert!(raw1.exists(), "raw file 1 missing: {:?}", raw1);
    assert!(raw2.exists(), "raw file 2 missing: {:?}", raw2);

    let tmp = tempfile::tempdir().unwrap();
    let tmp_path = tmp.path();

    // Stage 1: ingest
    run_atman(&[
        "ingest",
        "--platform", "olink-explore-ngs",
        "--parser", "dube",
        "--output-dir", tmp_path.to_str().unwrap(),
        raw1.to_str().unwrap(),
        raw2.to_str().unwrap(),
    ]);

    // Stage 2: qc
    run_atman(&[
        "qc",
        "--input-dir", tmp_path.to_str().unwrap(),
        "--output-dir", tmp_path.to_str().unwrap(),
        "--rule", "dube",
    ]);

    // Stage 3: matrix
    run_atman(&[
        "matrix",
        "--input-dir", tmp_path.to_str().unwrap(),
        "--output-dir", tmp_path.to_str().unwrap(),
        "--format", "dube-wide",
        "--split-by", "panel",
    ]);

    // Stage 4: fold change
    run_atman(&[
        "fold-change",
        "--input-dir", tmp_path.to_str().unwrap(),
        "--output-dir", tmp_path.to_str().unwrap(),
        "--groups", "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1",
    ]);

    // Stage 5a: strict diff of filtered NPX files
    for panel in PANELS {
        let ours = tmp_path.join(format!("{}_npx.csv", panel));
        let reference = data.join(format!("filtered_npx/npx/{}_npx.csv", panel));
        assert!(ours.exists(), "missing ours: {:?}", ours);
        assert!(reference.exists(), "missing reference: {:?}", reference);
        diff_strict(&ours, &reference).assert_clean();
    }

    // Stage 5b: numeric diff of fold change files
    for panel in PANELS {
        let ours = tmp_path.join(format!("{}_log2_fc.csv", panel));
        let reference = data.join(format!("fold_changes/fold_changes/{}_log2_fc.csv", panel));
        assert!(ours.exists(), "missing ours: {:?}", ours);
        assert!(reference.exists(), "missing reference: {:?}", reference);
        diff_numeric(&ours, &reference, 1e-4).assert_clean();
    }
}
```

- [ ] **Step 2: Run the reproduction test**

```bash
cargo test -p atman --test dube_reproduction -- --nocapture
```

Expected first run: **possibly FAIL** on the first attempt. Read the diagnostic (the diff harness prints up to 10 cell mismatches with row/column/left/right). Likely first-failure modes, in order of probability:

- **Gene-symbol column order disagreement on a single panel.** The spec claims Unicode-codepoint sort. Verified empirically for 4 panels during brainstorming, but only 4 of 8. If one of the other 4 panels uses a different order, the header row mismatches. Fix: look at the first mismatch column pair in the diff output and adjust sort rule.
- **Row order for control samples.** Dube places controls at the tail of each panel file in ingest order, not alphabetical. The pivot preserves ingest order for controls. If that's wrong for some panel, the control rows mismatch.
- **Panel slug transform:** `Cardiometabolic_II` → `cardiometabolic_ii_npx.csv` via plain ASCII lowercase. Confirm `to_ascii_lowercase()` is used, not `to_lowercase()` (the latter is locale-aware and slightly slower but should still match for ASCII panel names).
- **Fold-change missingness rule:** if the numeric diff fails on assays where at least one participant is missing in at least one exposure, the "available participants only" interpretation is correct (plan assumes this). If it fails on assays with full data, investigate the `mean(2^x) / mean(2^y)` path — specifically whether Dube does `log2(mean(2^x)) - log2(mean(2^y))` vs `log2(mean(2^x) / mean(2^y))` — these are mathematically equivalent in f64 except for catastrophic cancellation on huge values.
- **Float formatting of the wide NPX file is NOT a likely failure mode** — the plan writes source strings verbatim, avoiding any parse+reformat round-trip.

**Iteration workflow:**
1. Run with `--nocapture`: `cargo test -p atman --test dube_reproduction -- --nocapture`.
2. Read the first cell mismatch in the panic message. Note `row_key`, `column`, `left`, `right`.
3. If the mismatch is a **header** cell: the pivot assay-sort order is wrong. Inspect the reference header against your output header for that panel to find the boundary position.
4. If the mismatch is a **row key (first column)**: the row partition or sort is wrong. Inspect the sample classification logic.
5. If the mismatch is a **data cell**: trace back through ingest → qc → pivot for that specific `(panel, gene_symbol, sample_id)`. Add a temporary print at the pivot layer to confirm the expected value is present in `npx_source_str`.
6. Form a hypothesis, write a tiny unit test at the `atman-core` level reproducing the bug, fix, re-run.
7. Commit each fix with `fix(matrix): <what>` or `fix(ingest): <what>`.

For speed, run in debug mode during iteration (`cargo test -p atman --test dube_reproduction -- --nocapture`); use `--release` only for the final acceptance check (spec §15 runtime budget).

- [ ] **Step 3: Iterate to green**

This may take several iterations. Each fix:
1. Identify the specific diff (first cell mismatch).
2. Form a hypothesis.
3. Write a small unit test reproducing the bug at `atman-core` level if possible.
4. Fix.
5. Re-run integration test.
6. Commit.

Expected final outcome: all 16 files diff-clean in under 30 seconds.

- [ ] **Step 4: Final commit for the passing test**

Once green:

```bash
git add crates/atman/tests/dube_reproduction.rs
git commit -m "test(atman): end-to-end Dube reproduction — passing"
```

---

## Task 20: README, Lint, Fmt, Final Test Pass

**Files:**
- Create: `README.md`

- [ ] **Step 1: Write `README.md`**

```markdown
# Atman

A local-first Rust engine for proteomics data, extending the karna bioinformatics platform into the protein-abundance modality.

**v0.1 scope:** Ingest, QC-filter, and reproduce published Dube et al. *Scientific Data* 2023 Olink Explore NGS outputs — byte-for-byte on filtered NPX files and within 1e-4 on log2 fold-change files. No statistical testing, no normalization beyond what Olink already applied.

## Reproduction recipe

```bash
cargo build --release
cargo test --workspace --release
```

The workspace integration test `dube_reproduction` runs the full pipeline:
raw NPX CSVs → long TSV → QC → per-panel Dube-wide CSVs → per-panel log2 FC,
then diffs every output against the published reference files at
`example_data/dube_heat_2023/filtered_npx/` and `.../fold_changes/`.

## Manual pipeline

```bash
atman ingest \
    --platform olink-explore-ngs --parser dube \
    --output-dir out/ \
    example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv \
    example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv

atman qc --input-dir out/ --output-dir out/ --rule dube
atman matrix --input-dir out/ --output-dir out/ --format dube-wide --split-by panel
atman fold-change --input-dir out/ --output-dir out/ \
    --groups "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1"
```

## Design

See `docs/superpowers/specs/2026-04-14-atman-olink-reproduction-design.md`.
```

- [ ] **Step 2: Lint and format**

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
```

Fix any clippy findings (usually: unused imports, needless clones, format literals).

- [ ] **Step 3: Final test pass**

```bash
cargo test --workspace --release
```

All tests must pass. Target runtime: <30 seconds for the full suite.

- [ ] **Step 4: Commit**

```bash
git add README.md
git add -u
git commit -m "docs(readme): v0.1 reproduction recipe + final lint/fmt"
```

- [ ] **Step 5: Verify acceptance criteria (spec §15)**

```bash
cargo test --workspace --release          # criterion 1, 2
cargo clippy --workspace --all-targets -- -D warnings  # criterion 3
cargo fmt --check                          # criterion 4
```

All three must be green. That's v0.1 done.

---

## Completion Criteria

v0.1 is done when **all** of these hold:

1. `cargo test --workspace --release` passes on a clean checkout.
2. The `dube_reproduction` integration test passes (all 16 reference files diff-clean).
3. `cargo clippy --workspace --all-targets -- -D warnings` is clean.
4. `cargo fmt --check` is clean.
5. Reproduction test runtime < 30 seconds on a Mac Studio M-series.
6. README has the one-command reproduction recipe.
7. No `unwrap()` on a reachable CLI path produces a panic — errors flow through `Result`.

Each of the 20 tasks above ends with a commit. Total expected commit count: ~25 (some tasks split into multiple commits during the iteration loop on Task 19).

---

**End of plan.**
