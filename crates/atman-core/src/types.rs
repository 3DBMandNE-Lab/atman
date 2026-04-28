//! Core data model. Platform-agnostic. No I/O, no CSV, no filesystem.

use serde::{Deserialize, Serialize};
use std::str::FromStr;

/// A proteomics platform identifier. The current CLI adapter supports
/// `OlinkExploreNgs`.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Platform {
    OlinkExploreNgs,
    OlinkTargetQpcr,
    SomaScan,
    MaxQuantLfq,
    MaxQuantTmt,
    DiannReport,
    SpectronautReport,
    /// Any platform identifier not in atman's well-known list. Adapters
    /// declare their own platform string (e.g. `"cptac_tmt_proteome"`,
    /// `"thermo_proteome_discoverer"`); atman accepts and round-trips it
    /// without requiring a recompile to add the platform. The well-known
    /// variants exist only because some platforms (Olink NPX QC, MaxQuant
    /// LFQ matrix shape, etc.) have platform-specific code paths
    /// elsewhere; everything else is `Custom`.
    Custom(String),
}

impl Platform {
    /// Render the platform as it appears in the canonical TSV `platform`
    /// column. Returns an owned string because `Custom(String)` cannot be
    /// borrowed as `&'static str`.
    pub fn as_str(&self) -> &str {
        match self {
            Self::OlinkExploreNgs => "olink_explore_ngs",
            Self::OlinkTargetQpcr => "olink_target_qpcr",
            Self::SomaScan => "somascan",
            Self::MaxQuantLfq => "maxquant_lfq",
            Self::MaxQuantTmt => "maxquant_tmt",
            Self::DiannReport => "diann_report",
            Self::SpectronautReport => "spectronaut_report",
            Self::Custom(s) => s.as_str(),
        }
    }
}

impl FromStr for Platform {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let trimmed = s.trim();
        if trimmed.is_empty() {
            return Err("empty platform string".to_string());
        }
        match trimmed.to_ascii_lowercase().as_str() {
            "olink_explore_ngs" | "olink-explore-ngs" => Ok(Self::OlinkExploreNgs),
            "olink_target_qpcr" | "olink-target-qpcr" => Ok(Self::OlinkTargetQpcr),
            "somascan" | "soma_scan" => Ok(Self::SomaScan),
            "maxquant_lfq" | "maxquant-lfq" => Ok(Self::MaxQuantLfq),
            "maxquant_tmt" | "maxquant-tmt" | "tmt" => Ok(Self::MaxQuantTmt),
            "diann_report" | "diann-report" | "diann_dia" | "dia-nn" | "diann" => {
                Ok(Self::DiannReport)
            }
            "spectronaut_report" | "spectronaut-report" | "spectronaut_dia" => {
                Ok(Self::SpectronautReport)
            }
            // Anything else is accepted as a custom platform; the canonical
            // form preserves the operator's exact lowercased string so
            // round-trips through `as_str()` are stable.
            _ => Ok(Self::Custom(trimmed.to_ascii_lowercase())),
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

/// Peptide catalog entry for MS-based proteomics, keyed by
/// `(AssayId, peptide_id)`. `assay_id` points at the parent protein in
/// `proteins.tsv`; peptide-level differential abundance models aggregate
/// `PeptideMeasurementRecord` rows grouped by `assay_id`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PeptideIdentity {
    pub peptide_id: String,
    pub assay_id: AssayId,
    pub sequence: Option<String>,
    pub charge: Option<i32>,
    pub modifications: Option<String>,
    pub missed_cleavages: Option<u32>,
}

/// One sample × one peptide abundance measurement. Primary key is
/// `(sample_id, peptide_id)`. Abundance units are recorded as a free-form
/// string (for example `log2_intensity` or `log2_maxlfq`) so downstream
/// commands can reason about log-scale assumptions without the platform
/// enum churn required by `MeasurementRecord`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PeptideMeasurementRecord {
    pub sample_id: String,
    pub peptide_id: String,
    pub abundance: f64,
    pub abundance_unit: String,
    pub dropped_by_qc: bool,
    pub below_lod: bool,
}

impl PeptideMeasurementRecord {
    /// Returns `None` when the record is masked by QC, otherwise the
    /// f64 abundance. Mirrors `MeasurementRecord::effective_abundance`.
    pub fn effective_abundance(&self) -> Option<f64> {
        if self.dropped_by_qc || !self.abundance.is_finite() {
            None
        } else {
            Some(self.abundance)
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

#[cfg(test)]
mod tests {
    use super::Platform;

    #[test]
    fn platform_round_trips_through_string() {
        for variant in [
            Platform::OlinkExploreNgs,
            Platform::OlinkTargetQpcr,
            Platform::SomaScan,
            Platform::MaxQuantLfq,
            Platform::MaxQuantTmt,
            Platform::DiannReport,
            Platform::SpectronautReport,
            Platform::Custom("cptac_tmt_proteome".to_string()),
            Platform::Custom("thermo_proteome_discoverer".to_string()),
        ] {
            let s = variant.as_str().to_string();
            let parsed: Platform = s.parse().expect("canonical string parses back");
            assert_eq!(variant, parsed, "round-trip for {s}");
        }
    }

    #[test]
    fn unknown_platform_lands_in_custom_variant_without_recompile() {
        // Adapters declare their own platform string; atman accepts any
        // non-empty identifier without requiring it to be in the well-known
        // list. This is the extensibility seam — no recompile to add a
        // new platform.
        let parsed: Platform = "totally_made_up_platform".parse().unwrap();
        assert_eq!(
            parsed,
            Platform::Custom("totally_made_up_platform".to_string())
        );
        assert_eq!(parsed.as_str(), "totally_made_up_platform");
    }

    #[test]
    fn empty_platform_string_still_errors() {
        // The one rejection: an empty platform field is a schema error,
        // not a custom platform. Adapters must declare *something*.
        let err = "".parse::<Platform>().unwrap_err();
        assert!(err.contains("empty"), "got error: {err}");
        let err2 = "   ".parse::<Platform>().unwrap_err();
        assert!(err2.contains("empty"), "got error: {err2}");
    }

    #[test]
    fn maxquant_tmt_accepts_aliases() {
        for alias in [
            "maxquant_tmt",
            "maxquant-tmt",
            "tmt",
            "TMT",
            " MaxQuant_TMT ",
        ] {
            let parsed: Platform = alias.parse().unwrap_or_else(|e| panic!("{alias}: {e}"));
            assert_eq!(parsed, Platform::MaxQuantTmt);
        }
    }
}
