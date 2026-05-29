use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

/// Write a fixture with `n_plexes` plexes of `plex_size` samples each, and
/// `proteins_per_plex` proteins absent in each plex (and observed in all others).
/// The remaining `shared_proteins` proteins are observed in every sample.
fn write_plex_fixture(
    dir: &std::path::Path,
    n_plexes: usize,
    plex_size: usize,
    proteins_per_plex: usize,
    shared_proteins: usize,
    scramble: bool,
) {
    std::fs::create_dir_all(dir).unwrap();

    let total_samples = n_plexes * plex_size;
    let total_proteins = n_plexes * proteins_per_plex + shared_proteins;

    let mut samples_buf =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    let mut sample_ids: Vec<String> = Vec::new();
    let mut order = 1u64;
    for p in 0..n_plexes {
        for k in 0..plex_size {
            let sid = format!("S_{:02}_{:02}", p, k);
            sample_ids.push(sid.clone());
            samples_buf.push_str(&format!("{sid}\tsubj_{sid}\tcase\t0\ttumor\t{order}\n"));
            order += 1;
        }
    }
    std::fs::write(dir.join("samples.tsv"), samples_buf).unwrap();

    let mut proteins_buf =
        String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for i in 0..total_proteins {
        proteins_buf.push_str(&format!("diann_report\tA{:05}\t\tG{:05}\tms\t\n", i, i));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins_buf).unwrap();

    // Build per-sample absence sets.
    // Plex p's absence set = proteins indexed [p*proteins_per_plex .. (p+1)*proteins_per_plex).
    // Shared proteins are observed in every sample.
    let mut absent_for: Vec<std::collections::HashSet<usize>> =
        vec![std::collections::HashSet::new(); total_samples];
    for p in 0..n_plexes {
        let lo = p * proteins_per_plex;
        let hi = lo + proteins_per_plex;
        for k in 0..plex_size {
            let s = p * plex_size + k;
            for prot in lo..hi {
                absent_for[s].insert(prot);
            }
        }
    }
    if scramble {
        // Preserve per-sample absence count but reassign to random proteins.
        // Deterministic shuffle: linear-congruential pseudo-random with fixed seed.
        let mut rng_state: u64 = 0xDEADBEEFCAFEBABE;
        for s in 0..total_samples {
            let count = absent_for[s].len();
            absent_for[s].clear();
            while absent_for[s].len() < count {
                rng_state = rng_state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let pick = (rng_state >> 32) as usize % total_proteins;
                absent_for[s].insert(pick);
            }
        }
    }

    let mut buf = String::from("platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n");
    let mut row_order = 1u64;
    for (s_idx, sid) in sample_ids.iter().enumerate() {
        for prot in 0..total_proteins {
            if absent_for[s_idx].contains(&prot) {
                continue;
            }
            buf.push_str(&format!(
                "diann_report\t{sid}\tA{:05}\tG{:05}\tms\t\t10.0\t10.0\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{row_order}\n",
                prot, prot
            ));
            row_order += 1;
        }
    }
    std::fs::write(dir.join("measurements.tsv"), buf).unwrap();
}

fn parse_quality(path: &std::path::Path) -> std::collections::HashMap<String, String> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    let row = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    header
        .iter()
        .zip(row.iter())
        .map(|(h, v)| (h.to_string(), v.to_string()))
        .collect()
}

#[test]
fn recover_plex_finds_synthetic_plexes() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output_path = tmp.path().join("recovered.tsv");
    write_plex_fixture(&input, 4, 8, 50, 200, false);

    let result = run_atman(&[
        "recover-plex",
        "--input-dir",
        input.to_str().unwrap(),
        "--output",
        output_path.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    let body = std::fs::read_to_string(&output_path).unwrap();
    let n_data_rows = body.lines().count() - 1;
    assert_eq!(
        n_data_rows, 32,
        "expected 32 sample rows, got {n_data_rows}"
    );

    let quality_path = output_path.with_file_name("recovered_quality.tsv");
    let q = parse_quality(&quality_path);
    assert_eq!(q.get("n_clusters").map(String::as_str), Some("4"));
    assert_eq!(q.get("n_singletons").map(String::as_str), Some("0"));
    assert_eq!(
        q.get("median_cluster_size").map(String::as_str),
        Some("8.000000")
    );
    assert_eq!(
        q.get("diagnostic").map(String::as_str),
        Some("plex_coherent")
    );
    let separation: f64 = q.get("separation_ratio").unwrap().parse().unwrap();
    assert!(separation > 5.0, "separation ratio too low: {separation}");
}

#[test]
fn recover_plex_infers_pairs_within_plex() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    std::fs::create_dir_all(&input).unwrap();

    // Two plexes of 4 samples each: 2 tumor + 2 paired_non_tumor per plex,
    // pairs at consecutive numeric stems within each plex.
    // Plex 1 stems: T100, P101, T200, P201
    // Plex 2 stems: T500, P501, T600, P601
    let plex_assignments: Vec<(&str, &str, usize, &str)> = vec![
        ("T100", "tumor", 0, "100"),
        ("P101", "paired_non_tumor", 0, "101"),
        ("T200", "tumor", 0, "200"),
        ("P201", "paired_non_tumor", 0, "201"),
        ("T500", "tumor", 1, "500"),
        ("P501", "paired_non_tumor", 1, "501"),
        ("T600", "tumor", 1, "600"),
        ("P601", "paired_non_tumor", 1, "601"),
    ];

    let mut samples_buf =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for (i, (sid, cond, _plex, subj)) in plex_assignments.iter().enumerate() {
        samples_buf.push_str(&format!("{sid}\t{subj}\t{cond}\t0\ttumor\t{}\n", i + 1));
    }
    std::fs::write(input.join("samples.tsv"), samples_buf).unwrap();

    // 50 proteins absent in plex 0, 50 different proteins absent in plex 1, 200 shared.
    let total_proteins = 300;
    let mut proteins_buf =
        String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for i in 0..total_proteins {
        proteins_buf.push_str(&format!("diann_report\tA{:05}\t\tG{:05}\tms\t\n", i, i));
    }
    std::fs::write(input.join("proteins.tsv"), proteins_buf).unwrap();

    let mut buf = String::from("platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n");
    let mut row_order = 1u64;
    for (sid, _cond, plex, _subj) in &plex_assignments {
        for prot in 0..total_proteins {
            // Plex p's absent set: [p*50, (p+1)*50)
            if prot >= plex * 50 && prot < (plex + 1) * 50 {
                continue;
            }
            buf.push_str(&format!(
                "diann_report\t{sid}\tA{:05}\tG{:05}\tms\t\t10.0\t10.0\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{row_order}\n",
                prot, prot
            ));
            row_order += 1;
        }
    }
    std::fs::write(input.join("measurements.tsv"), buf).unwrap();

    let output_path = tmp.path().join("recovered.tsv");
    let augmented_path = tmp.path().join("samples_aug.tsv");
    let result = run_atman(&[
        "recover-plex",
        "--input-dir",
        input.to_str().unwrap(),
        "--output",
        output_path.to_str().unwrap(),
        "--augmented-samples-output",
        augmented_path.to_str().unwrap(),
        "--infer-pairs",
        "tumor-paired_non_tumor",
        "--expected-cluster-size",
        "4",
    ]);
    assert!(
        result.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    let pairs_path = output_path.with_file_name("recovered_pairs.tsv");
    let pairs_body = std::fs::read_to_string(&pairs_path).unwrap();
    let n_pairs = pairs_body.lines().count() - 1;
    assert_eq!(n_pairs, 4, "expected 4 patient pairs, got {n_pairs}");
    // Both members of every pair should share patient_id in augmented samples.
    let aug = std::fs::read_to_string(&augmented_path).unwrap();
    let mut header_iter = aug.lines();
    let header: Vec<&str> = header_iter.next().unwrap().split('\t').collect();
    let pid_col = header.iter().position(|h| *h == "patient_id").unwrap();
    let sid_col = header.iter().position(|h| *h == "sample_id").unwrap();
    let mut by_patient: std::collections::HashMap<String, Vec<String>> =
        std::collections::HashMap::new();
    for line in header_iter {
        let cols: Vec<&str> = line.split('\t').collect();
        by_patient
            .entry(cols[pid_col].to_string())
            .or_default()
            .push(cols[sid_col].to_string());
    }
    assert_eq!(by_patient.len(), 4);
    for (_pid, members) in &by_patient {
        assert_eq!(members.len(), 2, "every patient should have 2 samples");
    }
    // Specifically check that T100 and P101 are paired.
    let t100_pid = by_patient
        .iter()
        .find(|(_, v)| v.iter().any(|s| s == "T100"))
        .map(|(k, _)| k.clone())
        .unwrap();
    assert!(by_patient[&t100_pid].iter().any(|s| s == "P101"));

    let sidecar: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(output_path.with_extension("tsv.run.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(sidecar["command"], "recover-plex");
    assert!(sidecar["reinvoke"]
        .as_str()
        .unwrap()
        .starts_with("atman recover-plex "));
    let outputs = sidecar["output_files"].as_object().unwrap();
    assert!(outputs.contains_key(&pairs_path.display().to_string()));
    assert!(outputs.contains_key(&augmented_path.display().to_string()));
}

#[test]
fn recover_plex_flags_scrambled_absence() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output_path = tmp.path().join("recovered.tsv");
    write_plex_fixture(&input, 4, 8, 50, 200, true);

    let result = run_atman(&[
        "recover-plex",
        "--input-dir",
        input.to_str().unwrap(),
        "--output",
        output_path.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    let quality_path = output_path.with_file_name("recovered_quality.tsv");
    let q = parse_quality(&quality_path);
    let diagnostic = q.get("diagnostic").map(String::as_str).unwrap_or("");
    assert_ne!(
        diagnostic, "plex_coherent",
        "scrambled absence should not be classified plex_coherent (got {diagnostic})"
    );
}
