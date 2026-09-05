use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn parse(text: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(str::to_string)
        .collect();
    let rows = lines
        .map(|l| l.split('\t').map(str::to_string).collect())
        .collect();
    (header, rows)
}

fn col<'a>(header: &[String], rows: &'a [Vec<String>], name: &str) -> Vec<&'a str> {
    let i = header
        .iter()
        .position(|h| h == name)
        .unwrap_or_else(|| panic!("missing {name}"));
    rows.iter().map(|r| r[i].as_str()).collect()
}

const ACTIVATIONS: &str =
    "sample_id\tcohort\tcondition\tis_control\tprimary_A0001\tprimary_A0002\tresid_A0001\n\
S1\tc1\tCase\t0\t1.0\t1.0\t2.0\n\
S2\tc1\tCtrl\t1\t2.0\t2.0\t5.0\n\
S3\tc1\tCtrl\t1\t3.0\t3.0\t7.0\n\
S4\tc2\tCase\t0\t1.0\t4.0\t1.0\n\
S5\tc2\tCtrl\t1\t2.0\t5.0\t3.0\n\
S6\tc2\tCtrl\t1\t3.0\t6.0\t6.0\n";

#[test]
fn axes_build_orthogonalizes_within_cohort_and_zscores_globally() {
    let tmp = tempfile::tempdir().unwrap();
    let act = tmp.path().join("activations.tsv");
    std::fs::write(&act, ACTIVATIONS).unwrap();
    let out = tmp.path().join("axis_scores.tsv");
    let result = run_atman(&[
        "axes",
        "build",
        "--activations",
        act.to_str().unwrap(),
        "--representatives",
        "axis1=primary_A0002,axis2=resid_A0001",
        "--orthogonalize",
        "axis2",
        "--against",
        "axis1",
        "--within",
        "cohort",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let text = std::fs::read_to_string(&out).unwrap();
    let (header, rows) = parse(&text);
    assert_eq!(
        header,
        vec![
            "sample_id",
            "cohort",
            "condition",
            "is_control",
            "primary_A0001",
            "primary_A0002",
            "resid_A0001",
            "axis1_raw",
            "axis1_z",
            "axis2_unorth_raw",
            "axis2_unorth_z",
            "axis2_raw",
            "axis2_z"
        ]
    );
    assert_eq!(rows.len(), 6);
    assert_eq!(
        col(&header, &rows, "axis1_raw"),
        vec!["1", "2", "3", "4", "5", "6"]
    );
    assert_eq!(
        col(&header, &rows, "axis2_unorth_raw"),
        vec!["2", "5", "7", "1", "3", "6"]
    );
    let a2: Vec<f64> = col(&header, &rows, "axis2_raw")
        .iter()
        .map(|v| v.parse().unwrap())
        .collect();
    assert!((a2[0] + 1.0 / 6.0).abs() < 1e-9);
    assert!((a2[1] - 1.0 / 3.0).abs() < 1e-9);
    assert!((a2[2] + 1.0 / 6.0).abs() < 1e-9);
    assert!((a2[3] - 1.0 / 6.0).abs() < 1e-9);
    assert!((a2[4] + 1.0 / 3.0).abs() < 1e-9);
    assert!((a2[5] - 1.0 / 6.0).abs() < 1e-9);
    let z1: Vec<f64> = col(&header, &rows, "axis1_z")
        .iter()
        .map(|v| v.parse().unwrap())
        .collect();
    assert!((z1[0] + 1.3363062095621219).abs() < 1e-9);
    assert!((z1[5] - 1.3363062095621219).abs() < 1e-9);
    let sidecar = std::fs::read_to_string(tmp.path().join("axis_scores.tsv.run.json")).unwrap();
    assert!(sidecar.contains("\"command\": \"axes build\""));
    assert!(sidecar.contains("\"representatives\""));
}

#[test]
fn axes_build_joins_align_project_outputs_with_cohort_dirs() {
    let tmp = tempfile::tempdir().unwrap();
    let c1 = tmp.path().join("c1");
    let c2 = tmp.path().join("c2");
    std::fs::create_dir_all(&c1).unwrap();
    std::fs::create_dir_all(&c2).unwrap();
    std::fs::write(c1.join("samples.tsv"), "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\nS1\tS1\tCase\t0\tbio\t1\nS2\tS2\tCtrl\t1\tbio\t2\n").unwrap();
    std::fs::write(c2.join("samples.tsv"), "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\nS3\tS3\tCase\t0\tbio\t1\nS4\tS4\tCtrl\t1\tbio\t2\n").unwrap();
    let p1 = tmp.path().join("primary.tsv");
    let p2 = tmp.path().join("resid.tsv");
    std::fs::write(
        &p1,
        "sample_id\tA0001\tA0002\nS1\t1\t2\nS2\t3\t4\nS3\t5\t6\nS4\t7\t8\n",
    )
    .unwrap();
    std::fs::write(&p2, "sample_id\tA0001\nS1\t9\nS2\t8\nS3\t7\nS4\t6\n").unwrap();
    let out = tmp.path().join("axis_scores.tsv");
    let result = run_atman(&[
        "axes",
        "build",
        "--activations",
        &format!("primary={},resid={}", p1.display(), p2.display()),
        "--cohort-dirs",
        &format!("c1={},c2={}", c1.display(), c2.display()),
        "--representatives",
        "axis1=primary_A0002,axis3=resid_A0001",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let (header, rows) = parse(&std::fs::read_to_string(&out).unwrap());
    assert_eq!(
        &header[..7],
        &[
            "sample_id",
            "cohort",
            "condition",
            "is_control",
            "primary_A0001",
            "primary_A0002",
            "resid_A0001"
        ]
    );
    assert_eq!(col(&header, &rows, "cohort"), vec!["c1", "c1", "c2", "c2"]);
    assert_eq!(col(&header, &rows, "is_control"), vec!["0", "1", "0", "1"]);
    assert_eq!(col(&header, &rows, "axis3_raw"), vec!["9", "8", "7", "6"]);
    assert!(!header.iter().any(|h| h == "axis3_unorth_raw"));
}

#[test]
fn axes_build_rejects_unknown_representative_and_missing_against() {
    let tmp = tempfile::tempdir().unwrap();
    let act = tmp.path().join("activations.tsv");
    std::fs::write(&act, ACTIVATIONS).unwrap();
    let out = tmp.path().join("o.tsv");
    let r = run_atman(&[
        "axes",
        "build",
        "--activations",
        act.to_str().unwrap(),
        "--representatives",
        "axis1=nope",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("nope"));
    let r = run_atman(&[
        "axes",
        "build",
        "--activations",
        act.to_str().unwrap(),
        "--representatives",
        "axis1=primary_A0002,axis2=resid_A0001",
        "--orthogonalize",
        "axis2",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("--against"));
}
