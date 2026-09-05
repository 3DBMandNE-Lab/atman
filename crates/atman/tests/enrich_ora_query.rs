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
fn enrich_ora_query_tsv_tests_each_query_against_a_custom_universe() {
    let tmp = tempfile::tempdir().unwrap();
    let sets = tmp.path().join("sets.tsv");
    std::fs::write(
        &sets,
        "set_name\tgene_symbol\nS1\tG1\nS1\tG2\nS1\tG3\nS2\tG4\nS2\tG5\nS3\tG99\n",
    )
    .unwrap();
    let universe = tmp.path().join("universe.tsv");
    std::fs::write(
        &universe,
        "gene_symbol\nG1\nG2\nG3\nG4\nG5\nG6\nG7\nG8\nG9\nG10\n",
    )
    .unwrap();
    let queries = tmp.path().join("queries.tsv");
    std::fs::write(
        &queries,
        "query\tgene_symbol\nup\tG1\nup\tG2\ndown\tG4\ndown\tG42\n",
    )
    .unwrap();
    let out = tmp.path().join("ora.tsv");
    let r = run_atman(&[
        "enrich",
        "ora",
        "--query-tsv",
        queries.to_str().unwrap(),
        "--gene-sets",
        sets.to_str().unwrap(),
        "--universe",
        universe.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
    let rows = read_rows(&out);
    // S3 has no gene in the universe ⇒ skipped; 2 sets x 2 queries.
    assert_eq!(rows.len(), 4);
    let up_s1 = rows
        .iter()
        .find(|r| r["query"] == "up" && r["set_name"] == "S1")
        .unwrap();
    assert_eq!(up_s1["hit_count"], "2");
    assert_eq!(up_s1["overlap_size"], "2");
    assert_eq!(up_s1["universe_size"], "10");
    // P(X >= 2 | N = 10, K = 3, n = 2) = C(3,2)/C(10,2) = 3/45
    assert!((up_s1["p_value"].parse::<f64>().unwrap() - 3.0 / 45.0).abs() < 1e-12);
    let down_s2 = rows
        .iter()
        .find(|r| r["query"] == "down" && r["set_name"] == "S2")
        .unwrap();
    // G42 is outside the universe ⇒ one hit.
    assert_eq!(down_s2["hit_count"], "1");
    assert_eq!(down_s2["overlap_size"], "1");
    let sidecar = std::fs::read_to_string(tmp.path().join("ora.tsv.run.json")).unwrap();
    assert!(sidecar.contains("query_tsv"));
    assert!(sidecar.contains("\"universe\""));

    // --query-tsv without --universe is rejected.
    let r = run_atman(&[
        "enrich",
        "ora",
        "--query-tsv",
        queries.to_str().unwrap(),
        "--gene-sets",
        sets.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("--universe"));
}
