use proteome_core::ingest::olink_explore::OlinkExploreLongCsv;
use proteome_core::{Abundance, DubeSampleIdParser, Platform, ProteomeIngest};
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

fn write_fixture(dir: &TempDir, name: &str, content: &str) -> PathBuf {
    let p = dir.path().join(name);
    fs::write(&p, content).unwrap();
    p
}

const HEADER: &str =
    "SampleID;Index;OlinkID;UniProt;Assay;MissingFreq;Panel;Panel_Lot_Nr;PlateID;QC_Warning;LOD;NPX;Normalization;Assay_Warning";

#[test]
fn ingest_minimal_happy_path() {
    let dir = TempDir::new().unwrap();
    let rows = vec![
        "SSNA-001B-PR1;1;OID20838;Q9UJY5;GGA1;0.1;Neurology;B04414;PLATE1;PASS;0.20;0.15;Plate control;PASS",
        "SSNA-001B-PR2;2;OID20838;Q9UJY5;GGA1;0.1;Neurology;B04414;PLATE1;WARN;0.20;0.10;Plate control;PASS",
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
    assert!(out.measurements.iter().all(|m| m.platform == Platform::OlinkExploreNgs));

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
