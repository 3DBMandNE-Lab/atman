use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

/// Build a tiny canonical Atman input directory from synthetic data with two
/// recoverable non-gaussian sources. Returns the input dir path.
fn build_synthetic_canonical(tmp: &std::path::Path) -> std::path::PathBuf {
    let input_dir = tmp.join("canonical");
    std::fs::create_dir(&input_dir).unwrap();

    // Two independent uniform sources on [-1, 1].
    // Linear mix into 40 assays across 30 samples; ICA should recover both.
    let n_samples = 30;
    let n_assays = 40;
    let mut xs = Vec::with_capacity(n_samples);
    // Deterministic LCG for synthetic data (independent from atman's Xoshiro).
    let mut state: u64 = 1234567;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 33) as f64) / (u32::MAX as f64)
    };
    for _ in 0..n_samples {
        let s1 = (next() - 0.5) * 2.0;
        let s2 = (next() - 0.5) * 2.0;
        let mut row = Vec::with_capacity(n_assays);
        for j in 0..n_assays {
            let a1 = if j < n_assays / 2 { 1.0 } else { 0.1 };
            let a2 = if j >= n_assays / 2 { 1.0 } else { 0.1 };
            let noise = (next() - 0.5) * 0.05;
            row.push(a1 * s1 + a2 * s2 + noise);
        }
        xs.push(row);
    }

    // Write measurements.tsv in canonical long format.
    let headers = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
                    abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
                    detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order";
    let mut measurements = String::from(headers);
    measurements.push('\n');
    let mut ingest = 0u64;
    for (i, row) in xs.iter().enumerate() {
        let sample_id = format!("S{i:02}");
        for (j, v) in row.iter().enumerate() {
            let assay_id = format!("A{j:03}");
            let gene = format!("GENE{j:03}");
            ingest += 1;
            let src = format!("{v:.6}");
            measurements.push_str(&format!(
                "olink_explore_ngs\t{sample_id}\t{assay_id}\t{gene}\tP1\t{src}\t{v:.6}\t{v:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{ingest}\n",
            ));
        }
    }
    std::fs::write(input_dir.join("measurements.tsv"), &measurements).unwrap();

    // samples.tsv with subject_id == sample_id.
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 0..n_samples {
        samples.push_str(&format!("S{i:02}\tS{i:02}\tcase\t0\tcsf\t{}\n", i + 1));
    }
    std::fs::write(input_dir.join("samples.tsv"), samples).unwrap();

    // proteins.tsv (minimal).
    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 0..n_assays {
        proteins.push_str(&format!("olink_explore_ngs\tA{j:03}\t\tGENE{j:03}\tP1\t\n",));
    }
    std::fs::write(input_dir.join("proteins.tsv"), proteins).unwrap();

    input_dir
}

#[test]
fn decompose_ica_emits_loadings_activations_and_stability() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = build_synthetic_canonical(tmp.path());
    let loadings = tmp.path().join("loadings.tsv");
    let activations = tmp.path().join("activations.tsv");
    let stability = tmp.path().join("stability.tsv");

    let output = run_atman(&[
        "decompose",
        "ica",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--k",
        "2",
        "--n-seeds",
        "5",
        "--seed",
        "20260418",
        "--seed-stability-threshold",
        "0.5",
        "--stability-top-n",
        "10",
        "--output-loadings",
        loadings.to_str().unwrap(),
        "--output-activations",
        activations.to_str().unwrap(),
        "--output-stability",
        stability.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let loadings_text = std::fs::read_to_string(&loadings).unwrap();
    assert!(loadings_text.starts_with("program\tassay_id\tgene_symbol\tloading\n"));
    let loading_rows: Vec<&str> = loadings_text.lines().skip(1).collect();
    // 2 programs * 40 assays = 80 rows.
    assert_eq!(loading_rows.len(), 80);
    assert!(loading_rows.iter().any(|r| r.starts_with("program_01\t")));
    assert!(loading_rows.iter().any(|r| r.starts_with("program_02\t")));

    let activations_text = std::fs::read_to_string(&activations).unwrap();
    assert!(activations_text.starts_with("cohort\tsubject_id\tsample_id\tprogram\tactivation\n"));
    let activation_rows: Vec<&str> = activations_text.lines().skip(1).collect();
    // 2 programs * 30 samples = 60 rows.
    assert_eq!(activation_rows.len(), 60);

    let stability_text = std::fs::read_to_string(&stability).unwrap();
    assert!(stability_text.starts_with(
        "program\treference_seed\tn_alt_seeds\tn_stable_runs\tseed_stability_fraction\tmean_best_jaccard\tflag_below_threshold\n"
    ));
    let stability_rows: Vec<Vec<&str>> = stability_text
        .lines()
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(stability_rows.len(), 2);
    for row in &stability_rows {
        assert_eq!(row[1], "20260418");
        assert_eq!(row[2], "4");
        let frac: f64 = row[4].parse().unwrap();
        // With a clean linear mix and only 2 components, stability should be high
        // across seeds; require >= 0.5 recovered with top-10 Jaccard >= 0.5.
        assert!(frac >= 0.5, "program {} stability={frac}", row[0]);
    }

    // Sidecar: one-per-invocation provenance JSON next to the loadings TSV.
    let sidecar_path = tmp.path().join("loadings.tsv.run.json");
    let sidecar_text = std::fs::read_to_string(&sidecar_path)
        .expect("decompose ica should write <primary_output>.run.json");
    let sidecar: serde_json::Value = serde_json::from_str(&sidecar_text).unwrap();
    assert_eq!(sidecar["command"], "decompose ica");
    assert_eq!(sidecar["exit_code"], 0);
    assert_eq!(sidecar["args"]["k"], 2);
    assert_eq!(sidecar["args"]["n-seeds"], 5);
    assert_eq!(sidecar["args"]["seed"], 20260418);
    assert_eq!(sidecar["args"]["seed-stability-threshold"], 0.5);
    assert_eq!(sidecar["args"]["min-stable-seed-fraction"], 0.9);
    assert_eq!(sidecar["args"]["stability-metric"], "jaccard-top20");
    assert_eq!(sidecar["args"]["stability-top-n"], 10);
    assert!(sidecar["atman_version"].is_string());
    assert!(sidecar["atman_git_sha"].is_string());
    assert!(sidecar["os_arch"].is_string());
    // v1: inputs_sha256 is a per-file dict, `{path: "sha256:<hex>"}`.
    let inputs = sidecar["inputs_sha256"].as_object().expect("inputs_sha256");
    assert!(!inputs.is_empty(), "expected per-file input hashes");
    for (_k, v) in inputs {
        let s = v.as_str().unwrap();
        assert!(s.starts_with("sha256:"));
        assert_eq!(s.len(), "sha256:".len() + 64);
    }
    assert_eq!(sidecar["schema_version"], 1);
    assert!(sidecar["run_uuid"].as_str().unwrap().len() == 36);
    assert!(sidecar["reinvoke"]
        .as_str()
        .unwrap()
        .starts_with("atman decompose ica"));
    assert!(sidecar["build_env"]["rustc_version"].is_string());
    assert!(sidecar["build_env"]["cargo_lock_sha256"].is_string());
    assert!(sidecar["build_env"]["profile"].is_string());
    assert!(sidecar["cwd_at_start"].is_string());
    let started = sidecar["started_at"].as_str().unwrap();
    let finished = sidecar["finished_at"].as_str().unwrap();
    assert!(started.ends_with('Z') && started.len() == 20);
    assert!(finished.ends_with('Z') && finished.len() == 20);
    assert!(finished.as_bytes() >= started.as_bytes());
    // Every output file we declared should appear in output_files with a
    // sha256 prefix.
    let outputs = sidecar["output_files"].as_object().unwrap();
    assert_eq!(outputs.len(), 3);
    for (_, v) in outputs {
        assert!(v.as_str().unwrap().starts_with("sha256:"));
    }
}

#[test]
fn decompose_ica_is_deterministic_across_runs() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = build_synthetic_canonical(tmp.path());
    let loadings_a = tmp.path().join("loadings_a.tsv");
    let activations_a = tmp.path().join("activations_a.tsv");
    let stability_a = tmp.path().join("stability_a.tsv");
    let loadings_b = tmp.path().join("loadings_b.tsv");
    let activations_b = tmp.path().join("activations_b.tsv");
    let stability_b = tmp.path().join("stability_b.tsv");

    let args = |l: &str, a: &str, s: &str| {
        vec![
            "decompose".to_string(),
            "ica".to_string(),
            "--input-dir".to_string(),
            input_dir.to_str().unwrap().to_string(),
            "--k".to_string(),
            "2".to_string(),
            "--n-seeds".to_string(),
            "3".to_string(),
            "--seed".to_string(),
            "20260418".to_string(),
            "--stability-top-n".to_string(),
            "10".to_string(),
            "--output-loadings".to_string(),
            l.to_string(),
            "--output-activations".to_string(),
            a.to_string(),
            "--output-stability".to_string(),
            s.to_string(),
        ]
    };

    let aa: Vec<&str> = vec![
        "decompose",
        "ica",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--k",
        "2",
        "--n-seeds",
        "3",
        "--seed",
        "20260418",
        "--stability-top-n",
        "10",
        "--output-loadings",
        loadings_a.to_str().unwrap(),
        "--output-activations",
        activations_a.to_str().unwrap(),
        "--output-stability",
        stability_a.to_str().unwrap(),
    ];
    let bb: Vec<String> = args(
        loadings_b.to_str().unwrap(),
        activations_b.to_str().unwrap(),
        stability_b.to_str().unwrap(),
    );
    let bb_refs: Vec<&str> = bb.iter().map(|s| s.as_str()).collect();

    let out_a = run_atman(&aa);
    assert!(out_a.status.success());
    let out_b = run_atman(&bb_refs);
    assert!(out_b.status.success());

    assert_eq!(
        std::fs::read_to_string(&loadings_a).unwrap(),
        std::fs::read_to_string(&loadings_b).unwrap()
    );
    assert_eq!(
        std::fs::read_to_string(&activations_a).unwrap(),
        std::fs::read_to_string(&activations_b).unwrap()
    );
    assert_eq!(
        std::fs::read_to_string(&stability_a).unwrap(),
        std::fs::read_to_string(&stability_b).unwrap()
    );
}

#[test]
fn decompose_ica_min_stable_seed_fraction_zero_disables_flag() {
    // Lower bound regression: floor=0.0 makes the strict-less-than
    // comparison impossible, so flag_below_threshold must always read 0.
    // (Behavior under non-zero floors is data-dependent and asserted at
    // the sidecar level in the main test.)
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = build_synthetic_canonical(tmp.path());
    let stability = tmp.path().join("stability.tsv");
    let loadings = tmp.path().join("loadings.tsv");
    let activations = tmp.path().join("activations.tsv");

    let result = run_atman(&[
        "decompose",
        "ica",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--k",
        "2",
        "--n-seeds",
        "5",
        "--seed",
        "20260418",
        "--seed-stability-threshold",
        "0.999",
        "--min-stable-seed-fraction",
        "0.0",
        "--stability-top-n",
        "10",
        "--output-loadings",
        loadings.to_str().unwrap(),
        "--output-activations",
        activations.to_str().unwrap(),
        "--output-stability",
        stability.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let text = std::fs::read_to_string(&stability).unwrap();
    let flags: Vec<&str> = text
        .lines()
        .skip(1)
        .map(|line| line.split('\t').last().unwrap())
        .collect();
    assert!(!flags.is_empty());
    assert!(
        flags.iter().all(|f| *f == "0"),
        "expected no flag at floor=0.0, got {flags:?}"
    );
}

#[test]
fn decompose_ica_rejects_min_stable_seed_fraction_out_of_range() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = build_synthetic_canonical(tmp.path());
    let stability = tmp.path().join("stability.tsv");
    let loadings = tmp.path().join("loadings.tsv");
    let activations = tmp.path().join("activations.tsv");

    let result = run_atman(&[
        "decompose",
        "ica",
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--k",
        "2",
        "--n-seeds",
        "3",
        "--seed",
        "20260418",
        "--min-stable-seed-fraction",
        "1.5",
        "--output-loadings",
        loadings.to_str().unwrap(),
        "--output-activations",
        activations.to_str().unwrap(),
        "--output-stability",
        stability.to_str().unwrap(),
    ]);
    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("--min-stable-seed-fraction must be in [0, 1]"),
        "stderr was:\n{stderr}"
    );
}
