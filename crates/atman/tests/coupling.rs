use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

#[test]
fn coupling_reports_spearman_and_sign_consistency() {
    let tmp = tempfile::tempdir().unwrap();
    let activations = tmp.path().join("activations.tsv");
    let pairs = tmp.path().join("pairs.tsv");
    let out = tmp.path().join("coupling.tsv");
    let summary = tmp.path().join("summary.tsv");

    std::fs::write(
        &activations,
        "cohort\tsubject_id\tprogram\tactivation\n\
         c1\tS1\tprogram_01\t1\n\
         c1\tS1\tprogram_10\t4\n\
         c1\tS1\tprogram_20\t6\n\
         c1\tS2\tprogram_01\t2\n\
         c1\tS2\tprogram_10\t3\n\
         c1\tS2\tprogram_20\t5\n\
         c1\tS3\tprogram_01\t3\n\
         c1\tS3\tprogram_10\t2\n\
         c1\tS3\tprogram_20\t4\n\
         c1\tS4\tprogram_01\t4\n\
         c1\tS4\tprogram_10\t1\n\
         c1\tS4\tprogram_20\t3\n\
         c1\tS5\tprogram_01\t5\n\
         c1\tS5\tprogram_10\t0\n\
         c1\tS5\tprogram_20\t2\n\
         c2\tS1\tprogram_01\t1\n\
         c2\tS1\tprogram_04\t2\n\
         c2\tS2\tprogram_01\t2\n\
         c2\tS2\tprogram_04\t3\n\
         c2\tS3\tprogram_01\t3\n\
         c2\tS3\tprogram_04\t4\n\
         c2\tS4\tprogram_01\t4\n\
         c2\tS4\tprogram_04\t5\n",
    )
    .unwrap();
    std::fs::write(
        &pairs,
        "cohort\tpair_label\tprogs_a\tprogs_b\thypothesized_sign\n\
         c1\tsynaptic_vs_A01\tprogram_01\tprogram_10;program_20\tnegative\n\
         c2\tintrathecal_vs_humoral\tprogram_01\tprogram_04\tpositive\n",
    )
    .unwrap();

    let output = run_atman(&[
        "coupling",
        "--activations",
        activations.to_str().unwrap(),
        "--pairs",
        pairs.to_str().unwrap(),
        "--method",
        "spearman",
        "--report",
        "sign-consistency",
        "--output",
        out.to_str().unwrap(),
        "--output-summary",
        summary.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("cohort\tpair_label\tprogs_a\tprogs_b"));
    let rows: Vec<Vec<&str>> = text
        .lines()
        .skip(1)
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0][0], "c1");
    assert_eq!(rows[0][1], "synaptic_vs_A01");
    assert_eq!(rows[0][5], "5");
    assert!((rows[0][6].parse::<f64>().unwrap() + 1.0).abs() < 1e-12);
    assert_eq!(rows[0][10], "1");
    assert_eq!(rows[1][0], "c2");
    assert_eq!(rows[1][1], "intrathecal_vs_humoral");
    assert_eq!(rows[1][5], "4");
    assert!((rows[1][6].parse::<f64>().unwrap() - 1.0).abs() < 1e-12);
    assert_eq!(rows[1][10], "1");

    let summary_text = std::fs::read_to_string(summary).unwrap();
    assert!(summary_text.contains("all\tall\t2\t2\t1\t1\t0.25\n"));
    assert!(summary_text.contains("pair_label\tsynaptic_vs_A01\t1\t1\t0\t1\t0.5\n"));
}

#[test]
fn coupling_requires_complete_subject_support() {
    let tmp = tempfile::tempdir().unwrap();
    let activations = tmp.path().join("activations.tsv");
    let pairs = tmp.path().join("pairs.tsv");
    let out = tmp.path().join("coupling.tsv");
    let summary = tmp.path().join("summary.tsv");
    std::fs::write(
        &activations,
        "cohort\tsubject_id\tprogram\tactivation\n\
         c1\tS1\tA\t1\n\
         c1\tS1\tB\t2\n\
         c1\tS2\tA\t2\n\
         c1\tS2\tB\t1\n\
         c1\tS3\tA\t3\n",
    )
    .unwrap();
    std::fs::write(
        &pairs,
        "cohort\tpair_label\tprogs_a\tprogs_b\thypothesized_sign\n\
         c1\tbad\tA\tB\tnegative\n",
    )
    .unwrap();

    let output = run_atman(&[
        "coupling",
        "--activations",
        activations.to_str().unwrap(),
        "--pairs",
        pairs.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--output-summary",
        summary.to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("need at least 4"));
}
