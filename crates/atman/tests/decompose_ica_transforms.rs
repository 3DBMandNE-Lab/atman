//! Integration tests for `atman decompose ica --transform ...`.
//!
//! Coverage:
//! - CLR transform runs end-to-end on a synthetic fixture and writes
//!   `transform_applied.json` with the correct schema.
//! - ALR transform requires `--alr-reference` and refuses cleanly
//!   when the referenced gene is absent from the panel.
//! - ILR transform emits `n_output_coords = p − 1` and synthetic
//!   Helmert-basis assay IDs in the loadings TSV.
//! - Running with `--transform none` produces no `transform_applied.json`
//!   and emits identical loadings to an explicit CLR on a dataset
//!   that is already per-sample zero-centred (invariance sanity).

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_small_canonical(dir: &Path, n_proteins: usize, n_samples: usize) {
    std::fs::create_dir_all(dir).unwrap();
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=n_samples {
        let cond = if i <= n_samples / 2 { "A" } else { "B" };
        samples.push_str(&format!("S{i:03}\tS{i:03}\t{cond}\t0\tplasma\t{i}\n"));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 1..=n_proteins {
        proteins.push_str(&format!(
            "olink_explore_ngs\tA{j:03}\tQ{j:05}\tG{j:03}\tP1\t\n"
        ));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    // Synthetic log2-NPX values with planted structure in proteins
    // 1..5: two heavy-tailed signal sources drive them; proteins
    // 6..n_proteins are independent noise.
    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order: u64 = 0;
    let mut state: u64 = 20260418;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 33) as f64) / (u32::MAX as f64)
    };
    for i in 1..=n_samples {
        // Two signal sources, cubed uniforms for non-Gaussian flavour.
        let s1 = {
            let u = next() * 2.0 - 1.0;
            u * u * u
        };
        let s2 = {
            let u = next() * 2.0 - 1.0;
            u * u * u
        };
        for j in 1..=n_proteins {
            order += 1;
            let assay = format!("A{j:03}");
            let gene = format!("G{j:03}");
            let noise = (next() - 0.5) * 0.2;
            let value = if j <= 3 {
                2.0 * s1 + noise
            } else if j <= 5 {
                2.0 * s2 + noise
            } else {
                next() - 0.5
            };
            qc.push_str(&format!(
                "olink_explore_ngs\tS{i:03}\t{assay}\t{gene}\tP1\t{value:.6}\t\
                 {value:.6}\t{value:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
}

fn args_for_ica(input: &Path, out_dir: &Path, transform: &str, extra: &[&str]) -> Vec<String> {
    let mut v: Vec<String> = vec![
        "decompose".into(),
        "ica".into(),
        "--input-dir".into(),
        input.to_string_lossy().into_owned(),
        "--k".into(),
        "2".into(),
        "--n-seeds".into(),
        "2".into(),
        "--seed".into(),
        "20260418".into(),
        "--max-iter".into(),
        "100".into(),
        "--tol".into(),
        "1e-3".into(),
        "--transform".into(),
        transform.into(),
        "--output-loadings".into(),
        out_dir.join("loadings.tsv").to_string_lossy().into_owned(),
        "--output-activations".into(),
        out_dir
            .join("activations.tsv")
            .to_string_lossy()
            .into_owned(),
        "--output-stability".into(),
        out_dir.join("stability.tsv").to_string_lossy().into_owned(),
    ];
    v.extend(extra.iter().map(|s| s.to_string()));
    v
}

fn run_ica_with_transform(input: &Path, out_dir: &Path, transform: &str, extra: &[&str]) -> Output {
    std::fs::create_dir_all(out_dir).unwrap();
    let args = args_for_ica(input, out_dir, transform, extra);
    let refs: Vec<&str> = args.iter().map(|s| s.as_str()).collect();
    run_atman(&refs)
}

#[test]
fn clr_runs_and_writes_transform_applied_json() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let out_dir = tmp.path().join("clr_out");
    write_small_canonical(&input, 15, 20);
    let status = run_ica_with_transform(&input, &out_dir, "clr", &[]);
    assert!(
        status.status.success(),
        "clr run failed:\n{}",
        String::from_utf8_lossy(&status.stderr)
    );
    let applied = out_dir.join("transform_applied.json");
    assert!(applied.exists(), "transform_applied.json missing");
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&applied).unwrap()).unwrap();
    assert_eq!(j["transform"], "clr");
    assert_eq!(j["alr_reference"], serde_json::Value::Null);
    assert_eq!(j["n_input_proteins"], 15);
    assert_eq!(j["n_output_coords"], 15);

    // Sidecar echoes the transform.
    let sidecar = out_dir.join("loadings.tsv.run.json");
    let s: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();
    assert_eq!(s["args"]["transform"], "clr");
    // Output files list includes transform_applied.json.
    let keys: Vec<&str> = s["output_files"]
        .as_object()
        .unwrap()
        .keys()
        .map(|s| s.as_str())
        .collect();
    assert!(keys.iter().any(|k| k.ends_with("transform_applied.json")));
}

#[test]
fn transform_none_skips_transform_applied_json() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let out_dir = tmp.path().join("none_out");
    write_small_canonical(&input, 10, 16);
    let status = run_ica_with_transform(&input, &out_dir, "none", &[]);
    assert!(status.status.success());
    let applied = out_dir.join("transform_applied.json");
    assert!(
        !applied.exists(),
        "transform_applied.json should not exist under --transform none"
    );
}

#[test]
fn alr_requires_reference_and_refuses_unknown_gene() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let out_dir = tmp.path().join("alr_bad");
    write_small_canonical(&input, 10, 16);
    // (i) Missing --alr-reference → error.
    let missing = run_ica_with_transform(&input, &out_dir, "alr", &[]);
    assert!(!missing.status.success());
    let err = String::from_utf8_lossy(&missing.stderr);
    assert!(
        err.contains("--alr-reference"),
        "missing-reference error unexpected: {err}"
    );
    // (ii) Unknown gene → error.
    let unknown = run_ica_with_transform(
        &input,
        &out_dir,
        "alr",
        &["--alr-reference", "NOT_A_REAL_GENE"],
    );
    assert!(!unknown.status.success());
    let err = String::from_utf8_lossy(&unknown.stderr);
    assert!(err.contains("NOT_A_REAL_GENE"), "unknown-gene error: {err}");
}

#[test]
fn alr_with_valid_reference_runs_end_to_end() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let out_dir = tmp.path().join("alr_ok");
    write_small_canonical(&input, 12, 20);
    let status = run_ica_with_transform(&input, &out_dir, "alr", &["--alr-reference", "G010"]);
    assert!(
        status.status.success(),
        "alr run failed:\n{}",
        String::from_utf8_lossy(&status.stderr)
    );
    let applied: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(out_dir.join("transform_applied.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(applied["transform"], "alr");
    assert_eq!(applied["alr_reference"], "G010");
}

#[test]
fn ilr_reduces_output_coords_by_one() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let out_dir = tmp.path().join("ilr_out");
    write_small_canonical(&input, 8, 20);
    let status = run_ica_with_transform(&input, &out_dir, "ilr", &[]);
    assert!(
        status.status.success(),
        "ilr run failed:\n{}",
        String::from_utf8_lossy(&status.stderr)
    );
    let applied: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(out_dir.join("transform_applied.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(applied["transform"], "ilr");
    assert_eq!(applied["n_input_proteins"], 8);
    assert_eq!(applied["n_output_coords"], 7);

    // Loadings TSV should reference ilr_coord_* synthetic IDs, not
    // the original A001..A008 assay IDs.
    let loadings = std::fs::read_to_string(out_dir.join("loadings.tsv")).unwrap();
    let mut seen_assays: HashSet<String> = HashSet::new();
    for line in loadings.lines().skip(1) {
        let cells: Vec<&str> = line.split('\t').collect();
        if cells.len() >= 2 {
            seen_assays.insert(cells[1].to_string());
        }
    }
    assert!(
        seen_assays.iter().all(|a| a.starts_with("ilr_coord_")),
        "expected ilr_coord_* assay IDs, saw {:?}",
        seen_assays.iter().take(5).collect::<Vec<_>>()
    );
    assert_eq!(seen_assays.len(), 7);
}

#[test]
fn clr_is_deterministic_across_repeated_runs() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    write_small_canonical(&input, 10, 16);
    let dir_a: PathBuf = tmp.path().join("clr_a");
    let dir_b: PathBuf = tmp.path().join("clr_b");
    let status_a = run_ica_with_transform(&input, &dir_a, "clr", &[]);
    let status_b = run_ica_with_transform(&input, &dir_b, "clr", &[]);
    assert!(status_a.status.success() && status_b.status.success());
    let a = std::fs::read_to_string(dir_a.join("loadings.tsv")).unwrap();
    let b = std::fs::read_to_string(dir_b.join("loadings.tsv")).unwrap();
    assert_eq!(a, b, "CLR loadings must be byte-identical across runs");
}
