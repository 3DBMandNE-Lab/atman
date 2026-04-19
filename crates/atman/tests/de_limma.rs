use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_minimal_canonical(dir: &std::path::Path) {
    std::fs::create_dir_all(dir).unwrap();
    // Two groups, five subjects each, two features.
    let samples = "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n".to_string()
        + &(1..=10)
            .map(|i| {
                let cond = if i <= 5 { "A" } else { "B" };
                format!("S{i:02}\tS{i:02}\t{cond}\t0\tplasma\t{i}\n")
            })
            .collect::<String>();
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();
    let proteins = "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
                    olink_explore_ngs\tA001\tQ00001\tGENE1\tP1\t\n\
                    olink_explore_ngs\tA002\tQ00002\tGENE2\tP1\t\n";
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();
    let mut meas = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order = 0;
    for i in 1..=10 {
        for a in ["A001", "A002"] {
            order += 1;
            let effect = if a == "A001" { 1.0 } else { 0.0 };
            let shift = if i <= 5 { 0.0 } else { effect };
            let noise = ((i * 7 + order) as f64).sin() * 0.1;
            let v = shift + noise;
            meas.push_str(&format!(
                "olink_explore_ngs\tS{i:02}\t{a}\tGENE{a_suffix}\tP1\t{v:.6}\t\
                 {v:.6}\t{v:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n",
                a_suffix = &a[1..],
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &meas).unwrap();
    std::fs::write(dir.join("qc_measurements.tsv"), &meas).unwrap();
}

#[test]
fn atman_de_test_limma_smoke() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("de_out");
    write_minimal_canonical(&input);
    std::fs::create_dir_all(&output).unwrap();

    let out = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "limma",
        "--groups",
        "A-B",
        "--min-pairs",
        "3",
        "--trend",
        "false",
        "--robust",
        "false",
    ]);
    assert!(
        out.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let results = std::fs::read_to_string(output.join("de_results.tsv")).unwrap();
    assert!(results.lines().next().unwrap().contains("s2_posterior"));
    assert!(results.lines().next().unwrap().contains("lfc_threshold"));
    // Body rows should have a populated s2_posterior (non-empty column).
    let body = results.lines().nth(1).unwrap();
    let cols: Vec<&str> = body.split('\t').collect();
    let header: Vec<&str> = results.lines().next().unwrap().split('\t').collect();
    let idx = header.iter().position(|c| *c == "s2_posterior").unwrap();
    assert!(!cols[idx].is_empty());

    // Sidecar.
    let sidecar = output.join("de_results.tsv.run.json");
    assert!(sidecar.exists(), "sidecar missing");
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();
    assert_eq!(j["args"]["test"], "limma");
    assert_eq!(j["args"]["trend"], false);
    assert_eq!(j["args"]["robust"], false);
    assert_eq!(j["args"]["lfc-threshold"], 0.0);
}
