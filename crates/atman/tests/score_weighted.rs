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

const MEAS_HEADER: &str = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n";

fn measurement(sid: &str, gene: &str, v: f64, order: usize) -> String {
    format!("maxquant_lfq\t{sid}\t{gene}\t{gene}\tP\t{v}\t{v}\t{v}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}\n")
}

fn proteins(genes: &[&str]) -> String {
    let mut s = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for g in genes {
        s.push_str(&format!("maxquant_lfq\t{g}\tU{g}\t{g}\tP\t\n"));
    }
    s
}

#[test]
fn score_weighted_scores_subjects_and_summarizes_contrast() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    std::fs::create_dir_all(&input).unwrap();
    std::fs::write(input.join("samples.tsv"), "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\nS1\tS1\tCase\t0\tbio\t1\nS2\tS2\tCase\t0\tbio\t2\nS3\tS3\tCtrl\t1\tbio\t3\nS4\tS4\tCtrl\t1\tbio\t4\n").unwrap();
    std::fs::write(
        input.join("proteins.tsv"),
        proteins(&["GA", "GB", "GC", "GD"]),
    )
    .unwrap();
    let mut m = String::from(MEAS_HEADER);
    let ga = [1.0, 2.0, 3.0, 4.0];
    let gb = [4.0, 3.0, 2.0, 1.0];
    let mut order = 1;
    for (i, sid) in ["S1", "S2", "S3", "S4"].iter().enumerate() {
        m.push_str(&measurement(sid, "GA", ga[i], order));
        order += 1;
        m.push_str(&measurement(sid, "GB", gb[i], order));
        order += 1;
        m.push_str(&measurement(sid, "GC", 5.0, order));
        order += 1;
        if i == 0 {
            m.push_str(&measurement(sid, "GD", 9.0, order));
            order += 1;
        }
    }
    std::fs::write(input.join("measurements.tsv"), m).unwrap();
    let weights = tmp.path().join("weights.tsv");
    std::fs::write(
        &weights,
        "gene_symbol\tcohen_d\tfdr\nGA\t1\t0.1\nGB\t-1\t0.1\nGC\t2\t0.1\nGD\t5\t0.1\nGZ\t3\t0.1\n",
    )
    .unwrap();
    let out = tmp.path().join("scores.tsv");
    let summary = tmp.path().join("summary.tsv");
    let r = run_atman(&[
        "score",
        "weighted",
        "--input-dir",
        input.to_str().unwrap(),
        "--weights",
        weights.to_str().unwrap(),
        "--weight-col",
        "cohen_d",
        "--signature",
        "sig",
        "--max-missing-fraction",
        "0.5",
        "--min-shared",
        "2",
        "--groups",
        "Case-Ctrl",
        "--output",
        out.to_str().unwrap(),
        "--output-summary",
        summary.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&out);
    assert_eq!(rows.len(), 4);
    for r in &rows {
        assert_eq!(r["signature"], "sig");
        assert_eq!(r["n_shared"], "3");
        assert_eq!(r["n_used"], "2");
    }
    let s1 = rows.iter().find(|r| r["sample_id"] == "S1").unwrap();
    assert!((num(s1, "score") + 1.161895003862225).abs() < 1e-9);
    let s4 = rows.iter().find(|r| r["sample_id"] == "S4").unwrap();
    assert!((num(s4, "score") - 1.161895003862225).abs() < 1e-9);
    assert_eq!(s4["is_control"], "1");
    let sm = read_rows(&summary);
    assert_eq!(sm.len(), 1);
    assert_eq!(sm[0]["comparison"], "Case-Ctrl");
    assert_eq!(sm[0]["n_shared_proteins"], "3");
    assert_eq!(sm[0]["n_case"], "2");
    assert_eq!(sm[0]["n_control"], "2");
    assert!((num(&sm[0], "cohen_d") + 2.8284271247461903).abs() < 1e-9);
    assert!((num(&sm[0], "auc")).abs() < 1e-12);
    assert!(tmp.path().join("scores.tsv.run.json").exists());
}

#[test]
fn score_weighted_min_shared_blanks_score() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    std::fs::create_dir_all(&input).unwrap();
    std::fs::write(input.join("samples.tsv"), "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\nS1\tS1\tA\t0\tbio\t1\nS2\tS2\tA\t0\tbio\t2\nS3\tS3\tB\t1\tbio\t3\n").unwrap();
    std::fs::write(input.join("proteins.tsv"), proteins(&["GA"])).unwrap();
    let mut m = String::from(MEAS_HEADER);
    for (i, sid) in ["S1", "S2", "S3"].iter().enumerate() {
        m.push_str(&measurement(sid, "GA", i as f64, i + 1));
    }
    std::fs::write(input.join("measurements.tsv"), m).unwrap();
    let weights = tmp.path().join("w.tsv");
    std::fs::write(&weights, "gene_symbol\tcohen_d\nGA\t1\n").unwrap();
    let out = tmp.path().join("scores.tsv");
    let r = run_atman(&[
        "score",
        "weighted",
        "--input-dir",
        input.to_str().unwrap(),
        "--weights",
        weights.to_str().unwrap(),
        "--weight-col",
        "cohen_d",
        "--min-shared",
        "2",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&out);
    assert!(rows
        .iter()
        .all(|r| r["score"].is_empty() && r["n_used"] == "1"));
}

#[test]
fn score_weighted_collapse_genes() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    std::fs::create_dir_all(&input).unwrap();
    std::fs::write(input.join("samples.tsv"), "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\nS1\tS1\tA\t0\tbio\t1\nS2\tS2\tA\t0\tbio\t2\nS3\tS3\tB\t1\tbio\t3\nS4\tS4\tB\t1\tbio\t4\n").unwrap();
    std::fs::write(input.join("proteins.tsv"), "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\nmaxquant_lfq\tA1\tA1\tG\tP\t\nmaxquant_lfq\tA2\tA2\tG\tP\t\n").unwrap();
    let mut m = String::from(MEAS_HEADER);
    let mut order = 1;
    // G / A1 = [1,2,3,4]; G / A2 = [4,3,2,1] ⇒ per-sample mean is constant 2.5.
    for (i, sid) in ["S1", "S2", "S3", "S4"].iter().enumerate() {
        m.push_str(&format!("maxquant_lfq\t{sid}\tA1\tG\tP\t{v}\t{v}\t{v}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}\n", v = (i + 1) as f64));
        order += 1;
        m.push_str(&format!("maxquant_lfq\t{sid}\tA2\tG\tP\t{v}\t{v}\t{v}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{order}\n", v = (4 - i) as f64));
        order += 1;
    }
    std::fs::write(input.join("measurements.tsv"), m).unwrap();
    let weights = tmp.path().join("w.tsv");
    std::fs::write(&weights, "gene_symbol\tcohen_d\nG\t1\n").unwrap();
    let run = |rule: &str| {
        let out = tmp.path().join(format!("scores_{rule}.tsv"));
        let r = run_atman(&[
            "score",
            "weighted",
            "--input-dir",
            input.to_str().unwrap(),
            "--weights",
            weights.to_str().unwrap(),
            "--weight-col",
            "cohen_d",
            "--min-shared",
            "1",
            "--collapse-genes",
            rule,
            "--output",
            out.to_str().unwrap(),
        ]);
        assert!(
            r.status.success(),
            "stderr:\n{}",
            String::from_utf8_lossy(&r.stderr)
        );
        read_rows(&out)
    };
    let none = run("none");
    assert!((num(&none[0], "score") + 1.161895003862225).abs() < 1e-9);
    let sidecar = std::fs::read_to_string(tmp.path().join("scores_none.tsv.run.json")).unwrap();
    assert!(sidecar.contains("\"n_assays\": 2"), "{sidecar}");
    assert!(
        sidecar.contains("\"n_genes_after_collapse\": 1"),
        "{sidecar}"
    );
    let maxo = run("max-observed");
    assert_eq!(maxo[0]["score"], none[0]["score"]);
    let mean = run("mean");
    assert!(mean
        .iter()
        .all(|r| r["score"].is_empty() && r["n_used"] == "0"));
}
