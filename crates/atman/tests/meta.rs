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

fn write_program_de(path: &std::path::Path, rows: &[(&str, &str, &str, f64, f64, f64)]) {
    let mut text =
        String::from("cohort\tprogram\tcategory\tpoint_delta\tci_lo\tci_hi\tp_sign_stable\n");
    for (cohort, program, category, effect, lo, hi) in rows {
        text.push_str(&format!(
            "{cohort}\t{program}\t{category}\t{effect}\t{lo}\t{hi}\t0.05\n"
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

#[test]
fn meta_module_reports_program_sign_consistency() {
    let tmp = tempfile::tempdir().unwrap();
    let c1 = tmp.path().join("program1.tsv");
    let c2 = tmp.path().join("program2.tsv");
    let c3 = tmp.path().join("program3.tsv");
    let out = tmp.path().join("module_meta.tsv");
    write_program_de(
        &c1,
        &[
            ("c1", "A01", "humoral", 1.0, 0.5, 1.5),
            ("c1", "A02", "synaptic", -1.0, -1.5, -0.5),
        ],
    );
    write_program_de(
        &c2,
        &[
            ("c2", "A01", "humoral", 0.8, 0.3, 1.3),
            ("c2", "A03", "myeloid", 1.0, 0.5, 1.5),
        ],
    );
    write_program_de(&c3, &[("c3", "A01", "humoral", -0.2, -0.7, 0.3)]);

    let output = run_atman(&[
        "meta",
        "--inputs",
        &format!("{},{},{}", c1.display(), c2.display(), c3.display()),
        "--level",
        "module",
        "--report",
        "sign-consistency",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("program\tcategory\tn_cohorts"));
    assert!(text.contains("A01\thumoral\t3\t"));
    assert!(!text.contains("A02\tsynaptic"));
    assert!(!text.contains("A03\tmyeloid"));
    let row = text.lines().nth(1).unwrap();
    let fields: Vec<&str> = row.split('\t').collect();
    assert_eq!(fields[0], "A01");
    assert_eq!(fields[2], "3");
    assert_eq!(fields[16], "2");
    assert_eq!(fields[17], "1");
    assert!((fields[15].parse::<f64>().unwrap() - (2.0 / 3.0)).abs() < 1e-12);
}

#[test]
fn meta_module_sign_binomial_excludes_zero_effect_cohorts() {
    // Four cohorts on program A01: three positive, one with effect == 0.0.
    // Zero-effect cohorts must be excluded from BOTH the success count k and
    // the trial count n. So k = 3, n = 3 (not 4).
    //   corrected: binomial_upper_tail(3, 3) = C(3,3)/2^3       = 1/8  = 0.125
    //   old (bug): binomial_upper_tail(3, 4) = (C(4,3)+C(4,4))/2^4 = 5/16 = 0.3125
    let tmp = tempfile::tempdir().unwrap();
    let c1 = tmp.path().join("program1.tsv");
    let c2 = tmp.path().join("program2.tsv");
    let c3 = tmp.path().join("program3.tsv");
    let c4 = tmp.path().join("program4.tsv");
    let out = tmp.path().join("module_meta.tsv");
    write_program_de(&c1, &[("c1", "A01", "humoral", 1.0, 0.5, 1.5)]);
    write_program_de(&c2, &[("c2", "A01", "humoral", 0.8, 0.3, 1.3)]);
    write_program_de(&c3, &[("c3", "A01", "humoral", 1.2, 0.7, 1.7)]);
    write_program_de(&c4, &[("c4", "A01", "humoral", 0.0, -0.5, 0.5)]);

    let output = run_atman(&[
        "meta",
        "--inputs",
        &format!(
            "{},{},{},{}",
            c1.display(),
            c2.display(),
            c3.display(),
            c4.display()
        ),
        "--level",
        "module",
        "--report",
        "sign-consistency",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(out).unwrap();
    let row = text.lines().nth(1).unwrap();
    let fields: Vec<&str> = row.split('\t').collect();
    assert_eq!(fields[0], "A01");
    assert_eq!(fields[2], "4", "n_cohorts counts all four cohorts");
    assert_eq!(fields[16], "3", "n_positive");
    assert_eq!(fields[17], "0", "n_negative");
    let sign_binomial_p: f64 = fields[18].parse().unwrap();
    // Corrected value: trials = n_positive + n_negative = 3.
    assert!(
        (sign_binomial_p - 0.125).abs() < 1e-12,
        "expected corrected p=0.125, got {sign_binomial_p}"
    );
    // Confirm the fix changed behavior: the old buggy value was 0.3125.
    assert!(
        (sign_binomial_p - 0.3125).abs() > 1e-9,
        "p-value still matches the old buggy behavior"
    );
}
