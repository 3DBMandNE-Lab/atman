//! End-to-end test for the CPTAC TMT proteome Python adapter.
//!
//! Runs `adapters/cptac/cptac_tmt_proteome_to_atman.py` against a small
//! real-data subset of the CPTAC HCC proteome (`cptac_hcc_subset.tsv`) and
//! checks that the canonical TSVs land correctly and pass `atman validate`
//! with the recovered `tumor` / `paired_non_tumor` contrast.
//!
//! The fixture is a 4-gene × 4-sample slice of the published HCC proteome
//! (T112, P111, T113, P114; A1BG, A1CF, A2M, AAAS). It exercises every
//! load-bearing path in the adapter: leading summary rows (Mean, Median,
//! StdDev) skipped, `Log Ratio` / `Unshared Log Ratio` discrimination,
//! NCBIGeneID promotion to `assay_id`, and the `T<N>` / `P<N>` prefix
//! recovery into the `condition` column.

use std::path::PathBuf;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

#[test]
fn cptac_tmt_adapter_recovers_paired_condition_from_real_subset() {
    let root = workspace_root();
    let adapter = root.join("adapters/cptac/cptac_tmt_proteome_to_atman.py");
    let fixture =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/cptac_hcc_subset.tsv");
    assert!(adapter.exists(), "adapter missing: {adapter:?}");
    assert!(fixture.exists(), "fixture missing: {fixture:?}");

    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path();

    let py_status = Command::new("python3")
        .arg(&adapter)
        .arg("--proteome")
        .arg(&fixture)
        .arg("--tumor-tag")
        .arg("HCC")
        .arg("--output-dir")
        .arg(out)
        .arg("--sample-id-regex")
        .arg(r"^(?P<condition_key>[TP])(?P<subject>\d+)$")
        .arg("--condition-key-map")
        .arg("T=tumor,P=paired_non_tumor")
        .arg("--quiet")
        .status()
        .expect("run cptac_tmt_proteome_to_atman.py");
    assert!(py_status.success(), "adapter exited non-zero");

    // Canonical TSVs landed.
    let samples_path = out.join("samples.tsv");
    let proteins_path = out.join("proteins.tsv");
    let measurements_path = out.join("measurements.tsv");
    for f in [&samples_path, &proteins_path, &measurements_path] {
        assert!(f.exists(), "adapter did not write {f:?}");
    }

    // samples.tsv: 4 biological rows; T -> tumor, P -> paired_non_tumor;
    // subject_id is the numeric tail.
    let samples = std::fs::read_to_string(&samples_path).unwrap();
    let sample_lines: Vec<&str> = samples.lines().collect();
    assert_eq!(
        sample_lines[0],
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order"
    );
    assert_eq!(sample_lines.len(), 5, "expected 1 header + 4 sample rows");
    let by_sample: std::collections::BTreeMap<&str, Vec<&str>> = sample_lines[1..]
        .iter()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0], f)
        })
        .collect();
    assert_eq!(by_sample["T112"][1], "112", "subject_id for T112");
    assert_eq!(by_sample["T112"][2], "tumor", "T -> tumor");
    assert_eq!(by_sample["P111"][1], "111", "subject_id for P111");
    assert_eq!(
        by_sample["P111"][2], "paired_non_tumor",
        "P -> paired_non_tumor"
    );
    assert_eq!(by_sample["T113"][2], "tumor");
    assert_eq!(by_sample["P114"][2], "paired_non_tumor");

    // proteins.tsv: 4 distinct proteins; assay_id = NCBIGeneID.
    let proteins = std::fs::read_to_string(&proteins_path).unwrap();
    let protein_lines: Vec<&str> = proteins.lines().collect();
    assert_eq!(
        protein_lines[0],
        "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot"
    );
    assert_eq!(protein_lines.len(), 5, "expected 1 header + 4 protein rows");
    let by_gene: std::collections::BTreeMap<&str, Vec<&str>> = protein_lines[1..]
        .iter()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[3], f) // gene_symbol -> row
        })
        .collect();
    assert_eq!(by_gene["A1BG"][1], "1", "A1BG NCBIGeneID");
    assert_eq!(by_gene["A2M"][1], "2", "A2M NCBIGeneID");
    assert_eq!(by_gene["A1CF"][1], "29974", "A1CF NCBIGeneID");
    assert_eq!(by_gene["AAAS"][1], "8086", "AAAS NCBIGeneID");
    for (_, row) in &by_gene {
        assert_eq!(row[0], "cptac_tmt_proteome");
        assert_eq!(row[4], "CPTAC_HCC", "panel = CPTAC_<tumor-tag>");
    }

    // measurements.tsv: 4 samples × 4 proteins = 16 measurements.
    let measurements = std::fs::read_to_string(&measurements_path).unwrap();
    let m_lines: Vec<&str> = measurements.lines().collect();
    assert_eq!(m_lines.len(), 17, "expected 1 header + 16 measurement rows");

    // Spot-check one cell against the raw fixture: A1BG / T112 = 0.82994022391227695
    let a1bg_t112 = m_lines[1..]
        .iter()
        .find(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f[1] == "T112" && f[3] == "A1BG"
        })
        .expect("A1BG / T112 measurement missing");
    let f: Vec<&str> = a1bg_t112.split('\t').collect();
    assert_eq!(
        f[5], "0.82994022391227695",
        "npx_source_str should preserve full raw precision"
    );
    let abundance: f64 = f[6].parse().unwrap();
    assert!(
        (abundance - 0.82994022391227695).abs() < 1e-9,
        "abundance value drifted from raw"
    );

    // End-to-end: validate accepts the canonical output for the recovered
    // tumor vs. paired_non_tumor contrast. This is the auditability story
    // in test form: the schema-first stage does not silently drop the
    // condition information that downstream DE depends on.
    let v = run_atman(&[
        "validate",
        "--input-dir",
        out.to_str().unwrap(),
        "--groups",
        "tumor-paired_non_tumor",
        "--min-pairs",
        "2",
    ]);
    assert!(
        v.status.success(),
        "atman validate should accept the adapter output: stdout={}\nstderr={}",
        String::from_utf8_lossy(&v.stdout),
        String::from_utf8_lossy(&v.stderr),
    );
}

#[test]
fn cptac_tmt_adapter_default_condition_for_opaque_ids() {
    // For cohorts without inline condition encoding (GBM, BRCA, COAD,
    // HNSCC), the adapter falls back to --default-condition (tumor) with
    // subject_id == sample_id. Use the same HCC subset but omit the
    // sample-id regex; everything should be labeled `tumor`.
    let root = workspace_root();
    let adapter = root.join("adapters/cptac/cptac_tmt_proteome_to_atman.py");
    let fixture =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/cptac_hcc_subset.tsv");

    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path();

    let py_status = Command::new("python3")
        .arg(&adapter)
        .arg("--proteome")
        .arg(&fixture)
        .arg("--tumor-tag")
        .arg("OPAQUE")
        .arg("--output-dir")
        .arg(out)
        .arg("--quiet")
        .status()
        .expect("run cptac_tmt_proteome_to_atman.py");
    assert!(py_status.success());

    let samples = std::fs::read_to_string(out.join("samples.tsv")).unwrap();
    for line in samples.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            f[1], f[0],
            "subject_id should equal sample_id with no regex"
        );
        assert_eq!(f[2], "tumor", "default condition should apply");
    }
}
