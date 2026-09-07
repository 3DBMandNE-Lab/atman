use std::collections::HashMap;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn read_rows(path: &std::path::Path) -> Vec<HashMap<String, String>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|l| {
            header
                .iter()
                .zip(l.split('\t'))
                .map(|(h, v)| (h.to_string(), v.to_string()))
                .collect()
        })
        .collect()
}

#[test]
fn axes_icc_from_cohort_dir_subject_ids() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().join("healthy");
    std::fs::create_dir_all(&dir).unwrap();
    let mut scores = String::from("sample_id\tcohort\tcondition\tis_control\ts1\n");
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    // 3 subjects x 2 samples with subject means 1, 3, 5 and ±0.5 within;
    // plus one subject with a single sample (excluded).
    let data = [
        ("A", 0.5),
        ("A", 1.5),
        ("B", 2.5),
        ("B", 3.5),
        ("C", 4.5),
        ("C", 5.5),
        ("D", 9.0),
    ];
    for (i, (subject, v)) in data.iter().enumerate() {
        scores.push_str(&format!("R{i}\thealthy\thealthy\t1\t{v}\n"));
        samples.push_str(&format!("R{i}\t{subject}\thealthy\t1\tbio\t{}\n", i + 1));
    }
    let scores_path = tmp.path().join("scores.tsv");
    std::fs::write(&scores_path, scores).unwrap();
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();
    let out = tmp.path().join("icc.tsv");
    let r = run_atman(&[
        "axes",
        "icc",
        "--scores",
        scores_path.to_str().unwrap(),
        "--cohort-dirs",
        &format!("healthy={}", dir.display()),
        "--score-cols",
        "s1",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&out);
    assert_eq!(rows.len(), 1);
    let row = &rows[0];
    assert_eq!(row["cohort"], "healthy");
    assert_eq!(row["n_subjects"], "3");
    assert_eq!(row["n_samples"], "6");
    assert_eq!(row["k_mean"], "2");
    assert_eq!(row["k0"], "2");
    assert_eq!(row["n_singletons_excluded"], "1");
    // These two tolerances are tighter than six-decimal output can carry:
    // 7.5/8.5 written at `{:.6}` parses back 5.9e-8 away, which would fail
    // this assertion. They pass because the `axes` writers use
    // `io::format_f64` (`{}`, shortest round-trip), not `io::format_float`
    // (`{:.6}`). If an `axes` column is ever switched to the narrow writer,
    // this test fails and the fix is the tolerance, NOT the writer.
    assert!((row["icc1"].parse::<f64>().unwrap() - 7.5 / 8.5).abs() < 1e-12);
    assert!((row["within_sd"].parse::<f64>().unwrap() - 0.5f64.sqrt()).abs() < 1e-12);
    assert!(tmp.path().join("icc.tsv.run.json").exists());
}
