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

const SAMPLES: &str = "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
S1\tS1\tCtrl\t1\tbio\t1\nS2\tS2\tCtrl\t1\tbio\t2\nS3\tS3\tCase\t0\tbio\t3\nS4\tS4\tCase\t0\tbio\t4\nS5\tS5\tCase\t0\tbio\t5\n";

const MEAS_HEADER: &str = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n";

fn measurement(sid: &str, gene: &str, v: f64, order: usize) -> String {
    format!("maxquant_lfq\t{sid}\t{gene}\t{gene}\tP\t{v}\t{v}\t{v}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}\n")
}

#[test]
fn residuals_emits_variance_tables_and_canonical_dir() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    std::fs::create_dir_all(&input).unwrap();
    std::fs::write(input.join("samples.tsv"), SAMPLES).unwrap();
    std::fs::write(
        input.join("proteins.tsv"),
        "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\nmaxquant_lfq\tGA\tP1\tGA\tP\t\nmaxquant_lfq\tGB\tP2\tGB\tP\t\nmaxquant_lfq\tGC\tP3\tGC\tP\t\n",
    )
    .unwrap();
    let mut m = String::from(MEAS_HEADER);
    let axis = [1.0, 2.0, 3.0, 4.0, 5.0];
    let gb = [1.0, 3.0, 2.0, 3.0, 1.0];
    let mut order = 1;
    for (i, sid) in ["S1", "S2", "S3", "S4", "S5"].iter().enumerate() {
        m.push_str(&measurement(sid, "GA", 2.0 * axis[i], order));
        order += 1;
        m.push_str(&measurement(sid, "GB", gb[i], order));
        order += 1;
        if i < 2 {
            m.push_str(&measurement(sid, "GC", 1.0, order));
            order += 1;
        }
    }
    std::fs::write(input.join("measurements.tsv"), m).unwrap();
    let cov = tmp.path().join("axis.tsv");
    std::fs::write(
        &cov,
        "sample_id\taxis1_raw\nS1\t1\nS2\t2\nS3\t3\nS4\t4\nS5\t5\n",
    )
    .unwrap();
    let long = tmp.path().join("resid.tsv");
    let var = tmp.path().join("var.tsv");
    let summary = tmp.path().join("summary.tsv");
    let canon = tmp.path().join("resid_canonical");
    let r = run_atman(&[
        "residuals",
        "--input-dir",
        input.to_str().unwrap(),
        "--design",
        "~ axis1_raw",
        "--covariates-tsv",
        cov.to_str().unwrap(),
        "--max-missing-fraction",
        "0.5",
        "--min-samples",
        "3",
        "--output",
        long.to_str().unwrap(),
        "--output-variance",
        var.to_str().unwrap(),
        "--output-variance-summary",
        summary.to_str().unwrap(),
        "--output-canonical-dir",
        canon.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let v = read_rows(&var);
    assert_eq!(
        v.iter()
            .map(|r| r["gene_symbol"].as_str())
            .collect::<Vec<_>>(),
        vec!["GA", "GB"]
    );
    let ga = &v[0];
    assert_eq!(ga["n"], "5");
    assert!((ga["r2"].parse::<f64>().unwrap() - 1.0).abs() < 1e-12);
    assert!(v[1]["r2"].parse::<f64>().unwrap().abs() < 1e-12);
    let s = read_rows(&summary);
    assert_eq!(s.len(), 1);
    assert_eq!(s[0]["design"], "~ axis1_raw");
    assert_eq!(s[0]["n_samples"], "5");
    assert_eq!(s[0]["n_proteins"], "2");
    assert!((s[0]["frac_variance"].parse::<f64>().unwrap() - 40.0 / 44.0).abs() < 1e-12);
    assert!((s[0]["median_r2"].parse::<f64>().unwrap() - 0.5).abs() < 1e-12);
    assert!((s[0]["frac_r2_gt_0_25"].parse::<f64>().unwrap() - 0.5).abs() < 1e-12);
    let canon_m = std::fs::read_to_string(canon.join("measurements.tsv")).unwrap();
    assert_eq!(canon_m.lines().count(), 11);
    assert!(!canon_m.contains("\tGC\t"));
    // GA residuals are exactly zero.
    let ga_rows: Vec<&str> = canon_m.lines().filter(|l| l.contains("\tGA\t")).collect();
    assert_eq!(ga_rows.len(), 5);
    assert!(ga_rows.iter().all(|l| l.split('\t').nth(6) == Some("0")));
    assert_eq!(
        std::fs::read_to_string(canon.join("samples.tsv")).unwrap(),
        SAMPLES
    );
    assert!(canon.join("proteins.tsv").exists());
    let de = run_atman(&[
        "de",
        "--input-dir",
        canon.to_str().unwrap(),
        "--output-dir",
        tmp.path().join("de").to_str().unwrap(),
        "--test",
        "welch-t",
        "--groups",
        "Case-Ctrl",
        "--include-controls",
        "--min-pairs",
        "2",
    ]);
    assert!(
        de.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&de.stderr)
    );
    let sidecar = std::fs::read_to_string(tmp.path().join("resid.tsv.run.json")).unwrap();
    assert!(sidecar.contains("axis.tsv"));
    assert!(sidecar.contains("\"n_dropped_missingness\": 1"));
}

#[test]
fn residuals_accepts_expression_terms() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    std::fs::create_dir_all(&input).unwrap();
    std::fs::write(input.join("samples.tsv"), "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tage\nS1\tS1\tC\t1\tbio\t1\t10\nS2\tS2\tC\t1\tbio\t2\t100\nS3\tS3\tC\t1\tbio\t3\t1000\nS4\tS4\tC\t1\tbio\t4\t10000\n").unwrap();
    let mut m = String::from(MEAS_HEADER);
    for (i, sid) in ["S1", "S2", "S3", "S4"].iter().enumerate() {
        m.push_str(&measurement(sid, "GA", (i as f64 + 1.0) * 3.0, i + 1));
    }
    std::fs::write(input.join("measurements.tsv"), m).unwrap();
    let long = tmp.path().join("resid.tsv");
    let r = run_atman(&[
        "residuals",
        "--input-dir",
        input.to_str().unwrap(),
        "--design",
        "~ log10(age)",
        "--output",
        long.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let text = std::fs::read_to_string(&long).unwrap();
    assert!(text.lines().skip(1).all(|l| l.ends_with("\t0")), "{text}");
}
