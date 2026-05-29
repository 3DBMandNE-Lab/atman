use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

/// Build a fixture with `n_pairs` patient pairs, two proteins:
/// - "STRONG": tumor consistently 1.0 unit above non-tumor
/// - "FRAGILE": all pairs near zero except one outlier driving the effect
fn write_fixture(dir: &std::path::Path, n_pairs: usize) {
    std::fs::create_dir_all(dir).unwrap();
    let mut samples_buf = String::from(
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tpatient_id\n",
    );
    let mut order = 1u64;
    for k in 0..n_pairs {
        let pid = format!("pair_{:03}", k);
        samples_buf.push_str(&format!(
            "T{:03}\t{}\ttumor\t0\ttumor\t{}\t{}\n",
            k, k, order, pid
        ));
        order += 1;
        samples_buf.push_str(&format!(
            "P{:03}\t{}\tpaired_non_tumor\t0\tnormal\t{}\t{}\n",
            k, k, order, pid
        ));
        order += 1;
    }
    std::fs::write(dir.join("samples.tsv"), samples_buf).unwrap();

    let proteins = "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
diann_report\tA_STRONG\t\tSTRONG\tms\t\n\
diann_report\tA_FRAGILE\t\tFRAGILE\tms\t\n";
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    let mut buf = String::from("platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n");
    let mut row = 1u64;
    for k in 0..n_pairs {
        // STRONG: tumor 11.0 / normal 10.0 (Δ ≈ +1.0 with small per-pair jitter so
        // paired-t doesn't bail on zero variance).
        let jitter = ((k as f64) * 0.013).sin() * 0.05;
        let t_val = 11.0 + jitter;
        let n_val = 10.0 - jitter;
        buf.push_str(&format!(
            "diann_report\tT{:03}\tA_STRONG\tSTRONG\tms\t\t{}\t{}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            k, t_val, t_val, row
        ));
        row += 1;
        buf.push_str(&format!(
            "diann_report\tP{:03}\tA_STRONG\tSTRONG\tms\t\t{}\t{}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            k, n_val, n_val, row
        ));
        row += 1;
        // FRAGILE: pair 0 carries a Δ ≈ +20 outlier; every other pair has Δ ≈ 0
        // with small jitter so paired-t can compute when the outlier is dropped.
        let f_jitter = ((k as f64) * 0.27).cos() * 0.05;
        let (tf, nf) = if k == 0 {
            (30.0, 10.0)
        } else {
            (10.0 + f_jitter, 10.0 - f_jitter)
        };
        buf.push_str(&format!(
            "diann_report\tT{:03}\tA_FRAGILE\tFRAGILE\tms\t\t{}\t{}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            k, tf, tf, row
        ));
        row += 1;
        buf.push_str(&format!(
            "diann_report\tP{:03}\tA_FRAGILE\tFRAGILE\tms\t\t{}\t{}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            k, nf, nf, row
        ));
        row += 1;
    }
    std::fs::write(dir.join("measurements.tsv"), buf).unwrap();
}

#[test]
fn robust_paired_distinguishes_robust_from_fragile() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output_path = tmp.path().join("robust.tsv");
    write_fixture(&input, 30);

    let out = run_atman(&[
        "robust-paired",
        "--input-dir",
        input.to_str().unwrap(),
        "--groups",
        "tumor-paired_non_tumor",
        "--paired-by",
        "patient_id",
        "--min-pairs",
        "5",
        "--output",
        output_path.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let body = std::fs::read_to_string(&output_path).unwrap();
    let mut strong: Option<Vec<String>> = None;
    let mut fragile: Option<Vec<String>> = None;
    let header: Vec<String> = body
        .lines()
        .next()
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect();
    let col = |name: &str| header.iter().position(|h| h == name).unwrap();
    for line in body.lines().skip(1) {
        let cols: Vec<String> = line.split('\t').map(String::from).collect();
        match cols[col("gene_symbol")].as_str() {
            "STRONG" => strong = Some(cols),
            "FRAGILE" => fragile = Some(cols),
            _ => {}
        }
    }
    let strong = strong.expect("STRONG row missing");
    let fragile = fragile.expect("FRAGILE row missing");

    let strong_sign: f64 = strong[col("loso_sign_stability")].parse().unwrap();
    let strong_pstab: f64 = strong[col("loso_p_lt_05_stability")].parse().unwrap();
    assert!(
        strong_sign >= 0.999 && strong_pstab >= 0.999,
        "STRONG should be robust (sign={}, p_stab={})",
        strong_sign,
        strong_pstab
    );

    // FRAGILE: 1 outlier pair drives a Δ=+5/30 pop-out; dropping it should erase the effect.
    // Sign stability is 1 only because the outlier has |Δ|=5 vs others = 0. But p-value will
    // collapse when the outlier is dropped — so p_lt_05_stability should be 0.
    let fragile_pstab: f64 = fragile[col("loso_p_lt_05_stability")].parse().unwrap();
    assert!(
        fragile_pstab < 0.5,
        "FRAGILE should be unstable (p_stab={})",
        fragile_pstab
    );
    let influencer = &fragile[col("most_influential_pair_id")];
    assert_eq!(
        influencer, "pair_000",
        "outlier pair should be most influential, got {influencer}"
    );

    let sidecar: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(output_path.with_extension("tsv.run.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(sidecar["command"], "robust-paired");
    assert!(sidecar["reinvoke"]
        .as_str()
        .unwrap()
        .starts_with("atman robust-paired "));
}

#[test]
fn robust_paired_requires_pairing_column() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    write_fixture(&input, 10);

    // Strip patient_id column so the command must fail.
    let samples = std::fs::read_to_string(input.join("samples.tsv")).unwrap();
    let stripped: String = samples
        .lines()
        .map(|line| {
            let mut cols: Vec<&str> = line.split('\t').collect();
            cols.pop(); // drop last column (patient_id)
            cols.join("\t")
        })
        .collect::<Vec<_>>()
        .join("\n");
    std::fs::write(input.join("samples.tsv"), stripped + "\n").unwrap();

    let output_path = tmp.path().join("robust.tsv");
    let out = run_atman(&[
        "robust-paired",
        "--input-dir",
        input.to_str().unwrap(),
        "--groups",
        "tumor-paired_non_tumor",
        "--paired-by",
        "patient_id",
        "--output",
        output_path.to_str().unwrap(),
    ]);
    assert!(
        !out.status.success(),
        "expected failure when patient_id missing"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("patient_id"),
        "error message should mention patient_id; got: {stderr}"
    );
}
