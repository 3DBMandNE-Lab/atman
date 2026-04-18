//! Error types for atman-core.
use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum IngestError {
    #[error(
        "schema mismatch in {file}: missing columns {missing:?}, unexpected columns {extra:?}",
        file = file.display()
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

    #[error("invalid abundance value {raw:?} at {}:{line}", file.display())]
    InvalidAbundance {
        file: PathBuf,
        line: usize,
        raw: String,
    },

    #[error("missing required column {column} in {}", file.display())]
    MissingRequiredColumn { file: PathBuf, column: String },

    #[error("unexpected normalization value {value:?} at {}:{line} (expected 'Plate control')", file.display())]
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

    #[error("csv parse error in {}: {source}", file.display())]
    CsvError {
        file: PathBuf,
        #[source]
        source: csv::Error,
    },

    #[error("io error in {}: {source}", file.display())]
    IoError {
        file: PathBuf,
        #[source]
        source: std::io::Error,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

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
