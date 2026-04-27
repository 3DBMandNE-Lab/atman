use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_fixture(dir: &std::path::Path) -> std::path::PathBuf {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(
        dir.join("samples.tsv"),
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
C1\tC1\tControl\t0\tcsf\t1\n\
C2\tC2\tControl\t0\tcsf\t2\n\
K1\tK1\tCase\t0\tcsf\t3\n\
K2\tK2\tCase\t0\tcsf\t4\n",
    )
    .unwrap();
    std::fs::write(
        dir.join("proteins.tsv"),
        "\
platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
spectronaut_report\tP00001\tP00001\tGENE1\tspectronaut\t\n\
spectronaut_report\tP00002\tP00002\tGENE2\tspectronaut\t\n",
    )
    .unwrap();
    let rows = [
        ("C1", "P00001", "GENE1", 10.0),
        ("C1", "P00002", "GENE2", 12.0),
        ("C2", "P00001", "GENE1", 11.0),
        ("C2", "P00002", "GENE2", 13.0),
        ("K1", "P00001", "GENE1", 13.0),
        ("K1", "P00002", "GENE2", 15.0),
        ("K2", "P00001", "GENE1", 15.0),
        ("K2", "P00002", "GENE2", 17.0),
    ];
    let mut measurements = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    for (idx, (sample, assay, gene, value)) in rows.iter().enumerate() {
        measurements.push_str(&format!(
            "spectronaut_report\t{sample}\t{assay}\t{gene}\tspectronaut\t{value}\t{value}\t{value}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
            idx + 1
        ));
    }
    std::fs::write(dir.join("measurements.tsv"), measurements).unwrap();
    let modules = dir.join("modules.tsv");
    std::fs::write(
        &modules,
        "\
module\tgene_symbol\n\
heat\tGENE1\n\
heat\tGENE2\n\
heat\tMISSING\n",
    )
    .unwrap();
    modules
}

#[test]
fn score_modules_writes_coverage_and_canonical_de_inputs() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let out = tmp.path().join("module_scores.tsv");
    let canonical = tmp.path().join("canonical");
    let de_out = tmp.path().join("de");
    let modules = write_fixture(&input);

    let output = run_atman(&[
        "score",
        "modules",
        "--input-dir",
        input.to_str().unwrap(),
        "--modules-tsv",
        modules.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--method",
        "mean",
        "--canonical-output-dir",
        canonical.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let scores = std::fs::read_to_string(&out).unwrap();
    assert!(scores.contains("C1\tC1\tControl\theat\tmean\t11\t3\t2\t0.666"));
    assert!(canonical.join("measurements.tsv").exists());

    let de = run_atman(&[
        "de",
        "--input-dir",
        canonical.to_str().unwrap(),
        "--output-dir",
        de_out.to_str().unwrap(),
        "--test",
        "welch-t",
        "--groups",
        "Case-Control",
        "--min-pairs",
        "2",
    ]);
    assert!(
        de.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&de.stderr)
    );
    let de_results = std::fs::read_to_string(de_out.join("de_results.tsv")).unwrap();
    assert!(de_results.contains("module_score\theat\theat"));
}

#[test]
fn score_modules_supports_all_methods() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let modules = write_fixture(&input);
    for method in ["mean", "median", "zscore", "pc1"] {
        let out = tmp.path().join(format!("{method}.tsv"));
        let output = run_atman(&[
            "score",
            "modules",
            "--input-dir",
            input.to_str().unwrap(),
            "--modules-tsv",
            modules.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--method",
            method,
        ]);
        assert!(
            output.status.success(),
            "method={method} stderr:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let scores = std::fs::read_to_string(out).unwrap();
        assert!(scores.contains(&format!("\theat\t{method}\t")));
    }
}
