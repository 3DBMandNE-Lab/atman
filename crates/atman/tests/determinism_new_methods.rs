//! Integration test: Phase M — byte-identity determinism for three new methods.
//!
//! Verifies that `atman decompose nmf`, `atman decompose ica --missingness-model`,
//! and `atman de --test limma --adjust-for` produce byte-identical output across
//! two runs with pinned BLAS threading (`BLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`).

use std::path::Path;
use std::process::Command;

fn run_atman(args: &[&str]) -> std::process::Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    let mut cmd = Command::new(bin);
    // Pin BLAS threading to 1 to ensure determinism.
    cmd.env("BLAS_NUM_THREADS", "1")
        .env("OMP_NUM_THREADS", "1")
        .env("OPENBLAS_NUM_THREADS", "1");
    cmd.args(args).output().expect("run atman")
}

fn manifest_dir() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).to_path_buf()
}

/// Compare two TSV files byte-for-byte.
fn assert_tsv_byte_identical(path_a: &Path, path_b: &Path, label: &str) {
    let bytes_a = std::fs::read(path_a).unwrap_or_else(|_| panic!("read {}", path_a.display()));
    let bytes_b = std::fs::read(path_b).unwrap_or_else(|_| panic!("read {}", path_b.display()));
    assert_eq!(
        bytes_a,
        bytes_b,
        "{} TSV byte mismatch:\n  A: {}\n  B: {}",
        label,
        path_a.display(),
        path_b.display()
    );
}

/// Parse .run.json sidecar and extract input/output SHA fields for comparison.
fn extract_run_json_shas(path: &Path) -> (serde_json::Value, serde_json::Value) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|_| panic!("read {}", path.display()));
    let j: serde_json::Value =
        serde_json::from_str(&text).unwrap_or_else(|_| panic!("parse {}", path.display()));

    let input_sha = j["inputs_sha256"].clone();
    let output_files = j["output_files"].clone();

    (input_sha, output_files)
}

/// Compare .run.json sidecar SHA fields (but not timestamps or other mutable fields).
fn assert_run_json_shas_match(path_a: &Path, path_b: &Path, label: &str) {
    let (input_sha_a, outputs_a) = extract_run_json_shas(path_a);
    let (input_sha_b, outputs_b) = extract_run_json_shas(path_b);

    assert_eq!(
        input_sha_a, input_sha_b,
        "{} .run.json inputs_sha256 mismatch:\n  A: {}\n  B: {}",
        label, input_sha_a, input_sha_b
    );

    // Extract just the filenames and SHAs, not full paths (which differ between runs)
    let extract_filename_shas = |obj: &serde_json::Value| -> Vec<(String, String)> {
        let mut result = Vec::new();
        if let Some(map) = obj.as_object() {
            for (full_path, sha_val) in map {
                // Extract just the filename from the full path
                let filename = std::path::Path::new(full_path)
                    .file_name()
                    .map(|n| n.to_string_lossy().to_string())
                    .unwrap_or_else(|| full_path.clone());
                let sha = sha_val.as_str().unwrap_or("").to_string();
                result.push((filename, sha));
            }
        }
        result.sort();
        result
    };

    let shas_a = extract_filename_shas(&outputs_a);
    let shas_b = extract_filename_shas(&outputs_b);

    assert_eq!(
        shas_a, shas_b,
        "{} .run.json output file SHAs mismatch:\n  A: {:?}\n  B: {:?}",
        label, shas_a, shas_b
    );
}

/// Convert matrix TSV (sample_id, gene_symbol, abundance) to canonical Atman input dir
fn build_canonical_input_from_matrix(matrix_path: &Path, output_dir: &Path) {
    std::fs::create_dir_all(output_dir).unwrap();

    let matrix_text = std::fs::read_to_string(matrix_path).unwrap();
    let mut samples = std::collections::BTreeSet::new();
    let mut genes = Vec::new();
    let mut seen_genes = std::collections::BTreeSet::new();
    let mut measurements = Vec::new();

    let lines: Vec<&str> = matrix_text.lines().collect();
    let data_start = if lines[0].starts_with('#') { 2 } else { 1 };

    for line in &lines[data_start..] {
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        let sample_id = parts[0].to_string();
        let gene_symbol = parts[1].to_string();

        // Parse abundance, skip if parse fails
        let abundance = match parts[2].parse::<f64>() {
            Ok(v) => v,
            Err(_) => continue,
        };

        samples.insert(sample_id.clone());
        if seen_genes.insert(gene_symbol.clone()) {
            genes.push(gene_symbol.clone());
        }
        measurements.push((sample_id, gene_symbol, abundance));
    }

    // Write measurements.tsv
    let headers = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
                   abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
                   detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order";
    let mut measurements_tsv = String::from(headers);
    measurements_tsv.push('\n');
    for (ingest_order, (sample_id, gene_symbol, abundance)) in measurements.into_iter().enumerate()
    {
        let gene_idx = genes.iter().position(|g| g == &gene_symbol).unwrap();
        let assay_id = format!("A{:03}", gene_idx);
        let src = format!("{:.10}", abundance);
        measurements_tsv.push_str(&format!(
            "proteomics\t{}\t{}\t{}\tP1\t{}\t{:.10}\t{:.10}\tlog2_scale\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            sample_id, assay_id, gene_symbol, src, abundance, abundance, ingest_order + 1
        ));
    }
    std::fs::write(output_dir.join("measurements.tsv"), measurements_tsv).unwrap();

    // Write samples.tsv
    let mut samples_tsv =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for (i, sample_id) in samples.iter().enumerate() {
        samples_tsv.push_str(&format!(
            "{}\t{}\tcase\t0\tcsf\t{}\n",
            sample_id,
            sample_id,
            i + 1
        ));
    }
    std::fs::write(output_dir.join("samples.tsv"), samples_tsv).unwrap();

    // Write proteins.tsv
    let mut proteins_tsv =
        String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for (gene_idx, gene) in genes.iter().enumerate() {
        let assay_id = format!("A{:03}", gene_idx);
        proteins_tsv.push_str(&format!("proteomics\t{}\t\t{}\tP1\t\n", assay_id, gene));
    }
    std::fs::write(output_dir.join("proteins.tsv"), proteins_tsv).unwrap();
}

#[test]
fn decompose_nmf_is_byte_deterministic() {
    let fixture_path = manifest_dir().join("tests/fixtures/nmf_frobenius_input.tsv");
    assert!(
        fixture_path.exists(),
        "fixture not found: {}",
        fixture_path.display()
    );

    let tmp1 = tempfile::tempdir().unwrap();
    let tmp2 = tempfile::tempdir().unwrap();

    let input_1 = tmp1.path().join("input");
    let input_2 = tmp2.path().join("input");
    let output_loadings_1 = tmp1.path().join("nmf_loadings.tsv");
    let output_loadings_2 = tmp2.path().join("nmf_loadings.tsv");
    let output_activations_1 = tmp1.path().join("nmf_activations.tsv");
    let output_activations_2 = tmp2.path().join("nmf_activations.tsv");

    build_canonical_input_from_matrix(&fixture_path, &input_1);
    build_canonical_input_from_matrix(&fixture_path, &input_2);

    // Run 1
    let result_1 = run_atman(&[
        "decompose",
        "nmf",
        "--input-dir",
        input_1.to_str().unwrap(),
        "--k",
        "4",
        "--output-loadings",
        output_loadings_1.to_str().unwrap(),
        "--output-activations",
        output_activations_1.to_str().unwrap(),
    ]);
    assert!(
        result_1.status.success(),
        "run 1 failed:\n{}",
        String::from_utf8_lossy(&result_1.stderr)
    );

    // Run 2
    let result_2 = run_atman(&[
        "decompose",
        "nmf",
        "--input-dir",
        input_2.to_str().unwrap(),
        "--k",
        "4",
        "--output-loadings",
        output_loadings_2.to_str().unwrap(),
        "--output-activations",
        output_activations_2.to_str().unwrap(),
    ]);
    assert!(
        result_2.status.success(),
        "run 2 failed:\n{}",
        String::from_utf8_lossy(&result_2.stderr)
    );

    // Compare outputs byte-for-byte
    assert_tsv_byte_identical(&output_loadings_1, &output_loadings_2, "nmf_loadings");
    assert_tsv_byte_identical(
        &output_activations_1,
        &output_activations_2,
        "nmf_activations",
    );

    // Compare .run.json SHA fields
    let sidecar_1 = {
        let mut s = output_loadings_1.as_os_str().to_os_string();
        s.push(".run.json");
        std::path::PathBuf::from(s)
    };
    let sidecar_2 = {
        let mut s = output_loadings_2.as_os_str().to_os_string();
        s.push(".run.json");
        std::path::PathBuf::from(s)
    };
    assert_run_json_shas_match(&sidecar_1, &sidecar_2, "nmf_loadings sidecar");
}

#[test]
fn decompose_ica_missingness_is_byte_deterministic() {
    let fixture_dir = manifest_dir().join("tests/fixtures");

    let tmp1 = tempfile::tempdir().unwrap();
    let tmp2 = tempfile::tempdir().unwrap();

    let input_1 = tmp1.path().join("input");
    let input_2 = tmp2.path().join("input");
    let output_loadings_1 = tmp1.path().join("ica_loadings.tsv");
    let output_loadings_2 = tmp2.path().join("ica_loadings.tsv");
    let output_activations_1 = tmp1.path().join("ica_activations.tsv");
    let output_activations_2 = tmp2.path().join("ica_activations.tsv");
    let output_stability_1 = tmp1.path().join("ica_stability.tsv");
    let output_stability_2 = tmp2.path().join("ica_stability.tsv");

    build_canonical_input_from_matrix(&fixture_dir.join("mnar_observed_abundance.tsv"), &input_1);
    build_canonical_input_from_matrix(&fixture_dir.join("mnar_observed_abundance.tsv"), &input_2);

    // Run 1
    let result_1 = run_atman(&[
        "decompose",
        "ica",
        "--input-dir",
        input_1.to_str().unwrap(),
        "--missingness-model",
        "abundance-conditional",
        "--k",
        "4",
        "--output-loadings",
        output_loadings_1.to_str().unwrap(),
        "--output-activations",
        output_activations_1.to_str().unwrap(),
        "--output-stability",
        output_stability_1.to_str().unwrap(),
    ]);
    assert!(
        result_1.status.success(),
        "run 1 failed:\n{}",
        String::from_utf8_lossy(&result_1.stderr)
    );

    // Run 2
    let result_2 = run_atman(&[
        "decompose",
        "ica",
        "--input-dir",
        input_2.to_str().unwrap(),
        "--missingness-model",
        "abundance-conditional",
        "--k",
        "4",
        "--output-loadings",
        output_loadings_2.to_str().unwrap(),
        "--output-activations",
        output_activations_2.to_str().unwrap(),
        "--output-stability",
        output_stability_2.to_str().unwrap(),
    ]);
    assert!(
        result_2.status.success(),
        "run 2 failed:\n{}",
        String::from_utf8_lossy(&result_2.stderr)
    );

    // Compare outputs byte-for-byte
    assert_tsv_byte_identical(&output_loadings_1, &output_loadings_2, "ica_loadings");
    assert_tsv_byte_identical(
        &output_activations_1,
        &output_activations_2,
        "ica_activations",
    );

    // Compare .run.json SHA fields
    let sidecar_1 = {
        let mut s = output_loadings_1.as_os_str().to_os_string();
        s.push(".run.json");
        std::path::PathBuf::from(s)
    };
    let sidecar_2 = {
        let mut s = output_loadings_2.as_os_str().to_os_string();
        s.push(".run.json");
        std::path::PathBuf::from(s)
    };
    assert_run_json_shas_match(&sidecar_1, &sidecar_2, "ica_loadings sidecar");
}

#[test]
fn de_adjust_for_is_byte_deterministic() {
    let fixture_dir = manifest_dir().join("tests/fixtures");
    let input_dir = fixture_dir.join("de_adjust_for_input");
    let covariates_path = fixture_dir.join("de_adjust_for_covariates.tsv");

    assert!(
        input_dir.exists(),
        "input dir not found: {}",
        input_dir.display()
    );
    assert!(
        covariates_path.exists(),
        "covariates file not found: {}",
        covariates_path.display()
    );

    let tmp1 = tempfile::tempdir().unwrap();
    let tmp2 = tempfile::tempdir().unwrap();

    let output_dir_1 = tmp1.path();
    let output_dir_2 = tmp2.path();

    // Run 1
    let result_1 = run_atman(&[
        "de",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--output-dir",
        output_dir_1.to_str().unwrap(),
        "--test",
        "limma",
        "--groups",
        "b-a",
        "--adjust-for",
        covariates_path.to_str().unwrap(),
    ]);
    assert!(
        result_1.status.success(),
        "run 1 failed:\n{}",
        String::from_utf8_lossy(&result_1.stderr)
    );

    // Run 2
    let result_2 = run_atman(&[
        "de",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--output-dir",
        output_dir_2.to_str().unwrap(),
        "--test",
        "limma",
        "--groups",
        "b-a",
        "--adjust-for",
        covariates_path.to_str().unwrap(),
    ]);
    assert!(
        result_2.status.success(),
        "run 2 failed:\n{}",
        String::from_utf8_lossy(&result_2.stderr)
    );

    // Find the actual de_results TSV (it may be named de_results_limma.tsv or similar)
    let find_de_results = |dir: &Path| -> std::path::PathBuf {
        for entry in std::fs::read_dir(dir).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            if path.is_file()
                && path
                    .file_name()
                    .unwrap()
                    .to_string_lossy()
                    .starts_with("de_results")
                && path.extension().is_some_and(|ext| ext == "tsv")
            {
                return path;
            }
        }
        panic!("de_results*.tsv not found in {}", dir.display())
    };

    let actual_output_1 = find_de_results(output_dir_1);
    let actual_output_2 = find_de_results(output_dir_2);

    // Compare output byte-for-byte
    assert_tsv_byte_identical(&actual_output_1, &actual_output_2, "de_results");

    // Compare .run.json SHA fields
    let sidecar_1 = {
        let mut s = actual_output_1.as_os_str().to_os_string();
        s.push(".run.json");
        std::path::PathBuf::from(s)
    };
    let sidecar_2 = {
        let mut s = actual_output_2.as_os_str().to_os_string();
        s.push(".run.json");
        std::path::PathBuf::from(s)
    };
    assert_run_json_shas_match(&sidecar_1, &sidecar_2, "de_results sidecar");
}
