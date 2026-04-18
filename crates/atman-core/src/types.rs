//! Core data model. Platform-agnostic. No I/O, no CSV, no filesystem.

use serde::{Deserialize, Serialize};
use std::str::FromStr;

/// A proteomics platform identifier. The current CLI adapter supports
/// `OlinkExploreNgs`.
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

impl FromStr for Platform {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "olink_explore_ngs" | "olink-explore-ngs" => Ok(Self::OlinkExploreNgs),
            "olink_target_qpcr" | "olink-target-qpcr" => Ok(Self::OlinkTargetQpcr),
            "somascan" | "soma_scan" => Ok(Self::SomaScan),
            "maxquant_lfq" | "maxquant-lfq" => Ok(Self::MaxQuantLfq),
            "diann_report" | "diann-report" | "diann_dia" | "dia-nn" | "diann" => {
                Ok(Self::DiannReport)
            }
            "spectronaut_report" | "spectronaut-report" | "spectronaut_dia" => {
                Ok(Self::SpectronautReport)
            }
            other => Err(format!("unsupported platform {other:?}")),
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
/// intensity. The current Olink adapter emits `Log2Npx`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Abundance {
    Log2Npx(f64),
    Log2Intensity(f64),
    Ibaq(f64),
    Raw(f64),
}

impl Abundance {
    /// Returns the underlying f64 regardless of unit.
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
/// Carries BOTH the canonical f64 value (`abundance`) and the original NPX
/// source string (`npx_source_str`) so the Dube-wide CSV writer can emit
/// byte-exact source strings without any f64 parse/reformat round-trip.
/// Carries `gene_symbol` so the pivot and fold-change paths can key on gene
/// symbol without threading the protein catalog through every transformation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MeasurementRecord {
    pub platform: Platform,
    pub assay_id: AssayId,
    pub gene_symbol: Option<String>,
    pub sample_id: String,
    pub abundance: Abundance,
    pub abundance_raw: Abundance,
    pub npx_source_str: String,
    pub qc_sample: QcFlag,
    pub qc_assay: QcFlag,
    pub detection_limit: DetectionLimit,
    pub below_lod: bool,
    pub batch: Batch,
    pub dropped_by_qc: bool,
    pub ingest_order: u64,
    pub panel: Option<String>,
}

impl MeasurementRecord {
    /// Returns `None` when the record has been masked by QC, otherwise the
    /// f64 abundance. Every downstream consumer that needs the *effective*
    /// value should go through this accessor.
    pub fn effective_abundance(&self) -> Option<f64> {
        if self.dropped_by_qc {
            None
        } else {
            Some(self.abundance.as_f64())
        }
    }
}

/// A biological or control sample.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Sample {
    pub sample_id: String,
    pub subject_id: Option<String>,
    pub condition: Option<String>,
    pub is_control: bool,
    pub sample_type: Option<String>,
    pub ingest_order: u64,
}
