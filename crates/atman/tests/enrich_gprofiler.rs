//! Integration tests for `atman enrich gprofiler`. We never hit the live
//! endpoint: tests pre-populate the cache with canned g:Profiler JSON and run
//! with --offline so the command only exercises parsing + TSV writing.
//!
//! The cache key is a SHA-256 of the canonicalized request; tests derive it
//! by first running the command in --offline mode against an empty cache
//! (which fails with a key in the error message), then writing the canned
//! JSON under that key and re-running.

use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

const CANNED_JSON: &str = r#"{
  "result": [
    {
      "query": "program_01",
      "source": "GO:BP",
      "native": "GO:0001234",
      "name": "biological process",
      "p_value": 1.2e-5,
      "intersection_size": 2,
      "query_size": 3,
      "term_size": 40,
      "effective_domain_size": 20000,
      "intersections": [["AQP4","GFAP"]]
    },
    {
      "query": "program_01",
      "source": "KEGG",
      "native": "KEGG:hsa00000",
      "name": "example pathway",
      "p_value": 0.01,
      "intersection_size": 1,
      "query_size": 3,
      "term_size": 10,
      "effective_domain_size": 20000,
      "intersections": [["S100B"]]
    }
  ]
}"#;

#[test]
fn gprofiler_parses_cached_response_in_offline_mode() {
    let tmp = tempfile::tempdir().unwrap();
    let query_tsv = tmp.path().join("queries.tsv");
    let background = tmp.path().join("background.tsv");
    let cache_dir = tmp.path().join("cache");
    let output = tmp.path().join("enrichment.tsv");

    std::fs::write(
        &query_tsv,
        "query\tgene_symbol\nprogram_01\tAQP4\nprogram_01\tGFAP\nprogram_01\tS100B\n",
    )
    .unwrap();
    std::fs::write(&background, "gene_symbol\nAQP4\nGFAP\nS100B\nOTHER\n").unwrap();

    // First offline run: cache miss expected; error message contains the key.
    let miss = run_atman(&[
        "enrich",
        "gprofiler",
        "--query-tsv",
        query_tsv.to_str().unwrap(),
        "--background",
        background.to_str().unwrap(),
        "--organism",
        "hsapiens",
        "--sources",
        "GO:BP,KEGG",
        "--threshold-method",
        "fdr",
        "--user-threshold",
        "0.05",
        "--cache-dir",
        cache_dir.to_str().unwrap(),
        "--ontology-version",
        "test-pinned",
        "--offline",
        "--output",
        output.to_str().unwrap(),
    ]);
    assert!(!miss.status.success());
    let stderr = String::from_utf8_lossy(&miss.stderr).to_string();
    let key = extract_key(&stderr).expect("cache key in error");

    // Write canned response at that key and rerun.
    std::fs::create_dir_all(&cache_dir).unwrap();
    std::fs::write(cache_dir.join(format!("{key}.json")), CANNED_JSON).unwrap();

    let hit = run_atman(&[
        "enrich",
        "gprofiler",
        "--query-tsv",
        query_tsv.to_str().unwrap(),
        "--background",
        background.to_str().unwrap(),
        "--organism",
        "hsapiens",
        "--sources",
        "GO:BP,KEGG",
        "--threshold-method",
        "fdr",
        "--user-threshold",
        "0.05",
        "--cache-dir",
        cache_dir.to_str().unwrap(),
        "--ontology-version",
        "test-pinned",
        "--offline",
        "--output",
        output.to_str().unwrap(),
    ]);
    assert!(
        hit.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&hit.stderr)
    );

    let text = std::fs::read_to_string(&output).unwrap();
    assert!(text.starts_with(
        "query\tsource\tnative\tname\tp_value\tintersection_size\tquery_size\tterm_size\teffective_domain_size\tintersections\n"
    ));
    let rows: Vec<Vec<&str>> = text
        .lines()
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0][0], "program_01");
    assert_eq!(rows[0][1], "GO:BP");
    assert_eq!(rows[0][2], "GO:0001234");
    assert_eq!(rows[0][5], "2");
    assert_eq!(rows[1][1], "KEGG");
    assert_eq!(rows[1][9], "S100B");
}

#[test]
fn gprofiler_cache_miss_fails_in_offline_mode() {
    let tmp = tempfile::tempdir().unwrap();
    let query_tsv = tmp.path().join("queries.tsv");
    let background = tmp.path().join("background.tsv");
    let cache_dir = tmp.path().join("empty_cache");
    let output = tmp.path().join("enrichment.tsv");
    std::fs::write(&query_tsv, "query\tgene_symbol\np\tAQP4\n").unwrap();
    std::fs::write(&background, "AQP4\nGFAP\n").unwrap();

    let res = run_atman(&[
        "enrich",
        "gprofiler",
        "--query-tsv",
        query_tsv.to_str().unwrap(),
        "--background",
        background.to_str().unwrap(),
        "--cache-dir",
        cache_dir.to_str().unwrap(),
        "--offline",
        "--output",
        output.to_str().unwrap(),
    ]);
    assert!(!res.status.success());
    let stderr = String::from_utf8_lossy(&res.stderr);
    assert!(stderr.contains("cache miss"), "stderr: {}", stderr);
}

fn extract_key(stderr: &str) -> Option<String> {
    for token in stderr.split(|c: char| !c.is_ascii_hexdigit()) {
        if token.len() == 64 {
            return Some(token.to_string());
        }
    }
    None
}
