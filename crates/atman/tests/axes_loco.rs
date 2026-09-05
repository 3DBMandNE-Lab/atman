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

fn num(r: &HashMap<String, String>, c: &str) -> f64 {
    r[c].parse().unwrap_or_else(|_| panic!("{c} = {:?}", r[c]))
}

#[test]
fn axes_loco_recomputes_z_without_each_cohort() {
    let tmp = tempfile::tempdir().unwrap();
    let scores = tmp.path().join("scores.tsv");
    std::fs::write(
        &scores,
        "sample_id\tcohort\tcondition\tis_control\ta\n\
         P1\tc1\tCase\t0\t1\nP2\tc1\tCtrl\t1\t2\nP3\tc2\tCase\t0\t3\nP4\tc2\tCtrl\t1\t4\nP5\tc3\tCtrl\t1\t5\nP6\tc3\tCtrl\t1\t6\n",
    )
    .unwrap();
    let out = tmp.path().join("loco.tsv");
    let r = run_atman(&[
        "axes",
        "loco",
        "--scores",
        scores.to_str().unwrap(),
        "--score-cols",
        "a",
        "--groups",
        "Controls=Ctrl;Disease=Case",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&out);
    assert_eq!(rows.len(), 6);
    let full_z = |v: f64| (v - 3.5) / 1.8708286933869707;
    let d_c3 = rows
        .iter()
        .find(|r| r["group"] == "Disease" && r["dropped_cohort"] == "c3")
        .unwrap();
    assert_eq!(d_c3["n_remaining"], "2");
    assert!((num(d_c3, "full_z_mean") - (full_z(1.0) + full_z(3.0)) / 2.0).abs() < 1e-9);
    assert!((num(d_c3, "loo_z_mean") + 0.3872983346207417).abs() < 1e-9);
    assert!(
        (num(d_c3, "delta") - (num(d_c3, "loo_z_mean") - num(d_c3, "full_z_mean"))).abs() < 1e-12
    );
    let d_c1 = rows
        .iter()
        .find(|r| r["group"] == "Disease" && r["dropped_cohort"] == "c1")
        .unwrap();
    assert_eq!(d_c1["n_remaining"], "1");
    assert_eq!(
        rows.iter()
            .map(|r| r["dropped_cohort"].clone())
            .collect::<Vec<_>>()[..3],
        ["c1", "c2", "c3"]
    );
    assert!(tmp.path().join("loco.tsv.run.json").exists());
}
