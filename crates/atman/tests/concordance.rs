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

fn write_tables(dir: &std::path::Path) -> std::path::PathBuf {
    std::fs::write(dir.join("x_primary.tsv"), "gene_symbol\tcohen_d\tfdr\nG1\t1\t0.01\nG2\t2\t0.01\nG3\t3\t0.5\nG4\t4\t0.01\nG5\t5\t0.01\nG6\t9\t0.01\n").unwrap();
    std::fs::write(dir.join("y_primary.tsv"), "gene_symbol\tcohen_d\tfdr\nG1\t1\t0.01\nG2\t3\t0.5\nG3\t2\t0.01\nG4\t-5\t0.01\nG5\t4\t0.01\n").unwrap();
    std::fs::write(dir.join("x_removed.tsv"), "gene_symbol\tcohen_d\tfdr\nG1\t1\t0.01\nG2\t2\t0.01\nG3\t3\t0.5\nG4\t4\t0.01\nG5\t5\t0.01\n").unwrap();
    std::fs::write(dir.join("y_removed.tsv"), "gene_symbol\tcohen_d\tfdr\nG1\t4\t0.01\nG2\t5\t0.5\nG3\t2\t0.01\nG4\t3\t0.01\nG5\t1\t0.01\n").unwrap();
    std::fs::write(
        dir.join("loading.tsv"),
        "protein\tconsensus_loading\nG1\t0.1\nG2\t0.2\nG3\t0.3\nG4\t0.4\nG5\t0.5\n",
    )
    .unwrap();
    let manifest = dir.join("tables.tsv");
    std::fs::write(
        &manifest,
        format!(
            "label\tpath\teffect_col\tfeature_col\tq_col\tstage\n\
             x\t{}\tcohen_d\t\tfdr\tprimary\n\
             y\t{}\tcohen_d\t\tfdr\tprimary\n\
             x\t{}\tcohen_d\t\tfdr\tremoved\n\
             y\t{}\tcohen_d\t\tfdr\tremoved\n\
             axis1\t{}\tconsensus_loading\tprotein\t\tprimary\n",
            dir.join("x_primary.tsv").display(),
            dir.join("y_primary.tsv").display(),
            dir.join("x_removed.tsv").display(),
            dir.join("y_removed.tsv").display(),
            dir.join("loading.tsv").display()
        ),
    )
    .unwrap();
    manifest
}

#[test]
fn concordance_reports_rho_hits_jaccard_and_delta() {
    let tmp = tempfile::tempdir().unwrap();
    let manifest = write_tables(tmp.path());
    let out = tmp.path().join("conc.tsv");
    let delta = tmp.path().join("delta.tsv");
    let r = run_atman(&[
        "concordance",
        "--manifest",
        manifest.to_str().unwrap(),
        "--top-n",
        "2",
        "--q-threshold",
        "0.05",
        "--n-bootstrap",
        "200",
        "--seed",
        "5",
        "--output",
        out.to_str().unwrap(),
        "--output-delta",
        delta.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&out);
    assert_eq!(rows.len(), 4);
    let xy = rows
        .iter()
        .find(|r| r["a"] == "x" && r["b"] == "y" && r["stage"] == "primary")
        .unwrap();
    assert_eq!(xy["n_features"], "5");
    // scipy.stats.spearmanr([1,2,3,4,5],[1,3,2,-5,4]) -> rho 0.3, p 0.6238376647810728
    assert!((num(xy, "rho") - 0.3).abs() < 1e-12);
    assert!((num(xy, "p") - 0.6238376647810728).abs() < 1e-9);
    assert_eq!(xy["n_hits_both"], "3");
    assert!((num(xy, "sign_concordance") - 2.0 / 3.0).abs() < 1e-12);
    assert!((num(xy, "jaccard_top_n") - 1.0).abs() < 1e-12);
    assert!(num(xy, "ci_lo") <= 0.3 && 0.3 <= num(xy, "ci_hi"));
    let xa = rows
        .iter()
        .find(|r| r["a"] == "x" && r["b"] == "axis1")
        .unwrap();
    assert!((num(xa, "rho") - 1.0).abs() < 1e-12);
    assert!(xa["n_hits_both"].is_empty());
    let xy_r = rows
        .iter()
        .find(|r| r["a"] == "x" && r["b"] == "y" && r["stage"] == "removed")
        .unwrap();
    assert!((num(xy_r, "rho") + 0.8).abs() < 1e-12);
    let d = read_rows(&delta);
    assert_eq!(d.len(), 1);
    assert_eq!(d[0]["stage_from"], "primary");
    assert_eq!(d[0]["stage_to"], "removed");
    // removed: rho -0.8 ; delta = -0.8 - 0.3 = -1.1
    assert!((num(&d[0], "delta_rho") + 1.1).abs() < 1e-12);
    assert!(num(&d[0], "delta_lo") <= -1.1 + 1e-9 && -1.1 - 1e-9 <= num(&d[0], "delta_hi"));
    assert!(tmp.path().join("conc.tsv.run.json").exists());
}

#[test]
fn concordance_pairs_flag_restricts_output_and_is_deterministic() {
    let tmp = tempfile::tempdir().unwrap();
    let manifest = write_tables(tmp.path());
    let run = |p: &std::path::Path| {
        let r = run_atman(&[
            "concordance",
            "--manifest",
            manifest.to_str().unwrap(),
            "--pairs",
            "x:axis1",
            "--n-bootstrap",
            "50",
            "--seed",
            "1",
            "--output",
            p.to_str().unwrap(),
        ]);
        assert!(
            r.status.success(),
            "stderr:\n{}",
            String::from_utf8_lossy(&r.stderr)
        );
        std::fs::read_to_string(p).unwrap()
    };
    let a = run(&tmp.path().join("a.tsv"));
    let b = run(&tmp.path().join("b.tsv"));
    assert_eq!(a, b);
    assert_eq!(a.lines().count(), 2);
}
