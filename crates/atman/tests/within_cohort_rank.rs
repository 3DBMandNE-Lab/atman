use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

#[test]
fn within_cohort_rank_writes_ranked_canonical_directory() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output_dir = tmp.path().join("ranked");
    std::fs::create_dir_all(&input).unwrap();
    std::fs::write(
        input.join("samples.tsv"),
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
         S1\tS1\tCase\t0\tbio\t1\n\
         S2\tS2\tCase\t0\tbio\t2\n\
         S3\tS3\tControl\t0\tbio\t3\n\
         S4\tS4\tControl\t0\tbio\t4\n",
    )
    .unwrap();
    std::fs::write(
        input.join("proteins.tsv"),
        "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
         somascan\tA1\tP1\tG1\t\t\n\
         somascan\tA2\tP2\tG2\t\t\n",
    )
    .unwrap();
    let measurements = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n\
        somascan\tS1\tA1\tG1\t\t10\t10\t10\traw\tPASS\tPASS\t\t0\t0\t\t\t1\n\
        somascan\tS2\tA1\tG1\t\t30\t30\t30\traw\tPASS\tPASS\t\t0\t0\t\t\t2\n\
        somascan\tS3\tA1\tG1\t\t20\t20\t20\traw\tPASS\tPASS\t\t0\t0\t\t\t3\n\
        somascan\tS4\tA1\tG1\t\t20\t20\t20\traw\tPASS\tPASS\t\t0\t0\t\t\t4\n\
        somascan\tS1\tA2\tG2\t\t5\t5\t5\traw\tPASS\tPASS\t\t0\t0\t\t\t5\n\
        somascan\tS2\tA2\tG2\t\t6\t6\t6\traw\tPASS\tPASS\t\t0\t1\t\t\t6\n";
    std::fs::write(input.join("measurements.tsv"), measurements).unwrap();
    std::fs::write(input.join("qc_measurements.tsv"), measurements).unwrap();

    let output = run_atman(&[
        "within-cohort-rank",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let ranked = std::fs::read_to_string(output_dir.join("measurements.tsv")).unwrap();
    assert!(ranked.contains("somascan\tS1\tA1\tG1\t\t1\t1\t1\traw"));
    assert!(ranked.contains("somascan\tS2\tA1\tG1\t\t4\t4\t4\traw"));
    assert!(ranked.contains("somascan\tS3\tA1\tG1\t\t2.5\t2.5\t2.5\traw"));
    assert!(ranked.contains("somascan\tS4\tA1\tG1\t\t2.5\t2.5\t2.5\traw"));
    assert!(ranked.contains("somascan\tS2\tA2\tG2\t\t6\t\t6\traw\tPASS\tPASS\t\t0\t1"));
    assert!(output_dir.join("samples.tsv").exists());
    assert!(output_dir.join("proteins.tsv").exists());
    assert!(output_dir.join("qc_measurements.tsv").exists());
}
