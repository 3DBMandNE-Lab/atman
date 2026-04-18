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
    let content = format!(
        "{}\nSSNA-001B-PR1;1;OID1;P1;X;0.1;P;L;PL;PASS;0.1;0.2;Plate control\n",
        bad_header
    );
    let p = write(&dir, "f.csv", &content);
    let err = OlinkExploreLongCsv
        .read(&[p], &DubeSampleIdParser)
        .unwrap_err();
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
    let err = OlinkExploreLongCsv
        .read(&[p], &DubeSampleIdParser)
        .unwrap_err();
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
    let err = OlinkExploreLongCsv
        .read(&[p], &DubeSampleIdParser)
        .unwrap_err();
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
    let err = OlinkExploreLongCsv
        .read(&[p], &DubeSampleIdParser)
        .unwrap_err();
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
    let err = OlinkExploreLongCsv
        .read(&[p], &DubeSampleIdParser)
        .unwrap_err();
    assert!(matches!(err, IngestError::UnparseableSampleId { .. }));
}
