use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_de(path: &std::path::Path, rows: &[(&str, f64, f64, f64)]) {
    let mut text = String::from(
        "panel\tassay_id\tgene_symbol\tuniprot\tcomparison\tn_pairs\tmean_a\tmean_b\tmean_diff\tt\tdf\tp_value\tbh_q\tskip_reason\n",
    );
    for (gene, effect, t, p) in rows {
        text.push_str(&format!(
            "p\t{gene}\t{gene}\t\tCase-Control\t5\t\t\t{effect}\t{t}\t4\t{p}\t{p}\t\n"
        ));
    }
    std::fs::write(path, text).unwrap();
}

#[test]
fn meta_combines_shared_genes_and_skips_missing_only_genes() {
    let tmp = tempfile::tempdir().unwrap();
    let c1 = tmp.path().join("cohort1.tsv");
    let c2 = tmp.path().join("cohort2.tsv");
    let out = tmp.path().join("meta.tsv");
    write_de(&c1, &[("G1", 2.0, 4.0, 0.01), ("G2", -1.0, -2.0, 0.05)]);
    write_de(&c2, &[("G1", 1.0, 2.0, 0.05), ("G3", 3.0, 6.0, 0.001)]);

    let output = run_atman(&[
        "meta",
        "--inputs",
        &format!("{},{}", c1.display(), c2.display()),
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("comparison\tpanel\tgene_symbol\tn_cohorts"));
    assert!(text.contains("Case-Control\tp\tG1\t2\t"));
    assert!(!text.contains("\tG2\t"));
    assert!(!text.contains("\tG3\t"));
    assert!(text.contains("\t1\n") || text.contains("\t1."));
}

#[test]
fn meta_reports_sign_inconsistency_and_heterogeneity() {
    let tmp = tempfile::tempdir().unwrap();
    let c1 = tmp.path().join("cohort1.tsv");
    let c2 = tmp.path().join("cohort2.tsv");
    let c3 = tmp.path().join("cohort3.tsv");
    let out = tmp.path().join("meta.tsv");
    write_de(&c1, &[("G1", 2.0, 4.0, 0.01)]);
    write_de(&c2, &[("G1", -1.0, -2.0, 0.05)]);
    write_de(&c3, &[("G1", 1.5, 3.0, 0.02)]);

    let output = run_atman(&[
        "meta",
        "--inputs",
        &format!("{},{},{}", c1.display(), c2.display(), c3.display()),
        "--output",
        out.to_str().unwrap(),
        "--method",
        "all",
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(out).unwrap();
    let row = text.lines().nth(1).expect("meta row");
    let fields: Vec<&str> = row.split('\t').collect();
    let q_heterogeneity: f64 = fields[14].parse().unwrap();
    let sign_consistency: f64 = fields[18].parse().unwrap();
    assert!(q_heterogeneity > 0.0);
    assert!((sign_consistency - (2.0 / 3.0)).abs() < 1e-12);
}
