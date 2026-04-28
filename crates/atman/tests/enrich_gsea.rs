use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_de_results(path: &std::path::Path, top_genes: &[&str], rest: &[&str]) {
    let mut s = String::from(
        "panel\tassay_id\tgene_symbol\tuniprot\tcomparison\tn_pairs\tmean_a\tmean_b\tmean_diff\tt\tdf\tp_value\tbh_q\tskip_reason\n",
    );
    // Top hits: high positive mean_diff, very low p.
    for (i, g) in top_genes.iter().enumerate() {
        let mean_diff = 3.0 - 0.05 * (i as f64);
        let p = 1e-8 * (10f64.powf(0.5 * (i as f64)));
        s.push_str(&format!(
            "p\t{g}\t{g}\t\tCase-Control\t10\t\t\t{:.6}\t{:.4}\t8\t{:.3e}\t{:.3e}\t\n",
            mean_diff, mean_diff * 4.0, p, p
        ));
    }
    // Background: alternating-sign small effects, high p — symmetric around zero
    // so any randomly-chosen subset is uninformative under GSEA.
    for (i, g) in rest.iter().enumerate() {
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
        let mean_diff = sign * (0.10 - 0.0003 * (i as f64));
        let p = (0.4 + 0.001 * (i as f64)).min(0.999);
        s.push_str(&format!(
            "p\t{g}\t{g}\t\tCase-Control\t10\t\t\t{:.6}\t{:.4}\t8\t{:.6}\t{:.6}\t\n",
            mean_diff, mean_diff * 0.5, p, p
        ));
    }
    std::fs::write(path, s).unwrap();
}

#[test]
fn enrich_gsea_top_loaded_set_is_significant() {
    let tmp = tempfile::tempdir().unwrap();
    let de = tmp.path().join("de_results.tsv");
    let sets = tmp.path().join("gene_sets.tsv");
    let out = tmp.path().join("gsea.tsv");
    let top: Vec<String> = (0..15).map(|i| format!("HIT{:03}", i)).collect();
    let rest: Vec<String> = (0..200).map(|i| format!("BG{:03}", i)).collect();
    let top_refs: Vec<&str> = top.iter().map(String::as_str).collect();
    let rest_refs: Vec<&str> = rest.iter().map(String::as_str).collect();
    write_de_results(&de, &top_refs, &rest_refs);

    // Two sets: one made entirely of top-loaded hits, one made of background
    // genes spread evenly through the ranked list (every 20th position).
    let mut sets_text = String::from("set_name\tgene_symbol\n");
    for g in &top[0..10] {
        sets_text.push_str(&format!("hit_panel\t{g}\n"));
    }
    for i in (0..200).step_by(20) {
        sets_text.push_str(&format!("background_panel\tBG{:03}\n", i));
    }
    std::fs::write(&sets, sets_text).unwrap();

    let output = run_atman(&[
        "enrich",
        "gsea",
        "--de-results",
        de.to_str().unwrap(),
        "--gene-sets",
        sets.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--n-permutations",
        "500",
        "--seed",
        "13",
        "--min-set-size",
        "5",
    ]);
    assert!(
        output.status.success(),
        "atman enrich gsea failed: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );

    let body = std::fs::read_to_string(&out).unwrap();
    let lines: Vec<&str> = body.lines().collect();
    assert_eq!(
        lines[0], "set_name\tset_size\tes\tnes\tp_value\tbh_q\tleading_edge",
        "header mismatch"
    );
    let rows: Vec<Vec<&str>> = lines[1..].iter().map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows.len(), 2, "expected one row per set");

    // Row order: hit_panel should sort first (smaller q, larger |NES|).
    assert_eq!(rows[0][0], "hit_panel", "hit_panel should rank first");
    let hit_es: f64 = rows[0][2].parse().unwrap();
    let hit_p: f64 = rows[0][4].parse().unwrap();
    assert!(
        hit_es > 0.9,
        "hit_panel ES should be ≈ 1.0 (top-loaded), got {hit_es}"
    );
    assert!(
        hit_p < 0.01,
        "hit_panel p-value should be < 0.01 at 500 perms, got {hit_p}"
    );

    let bg_p: f64 = rows[1][4].parse().unwrap();
    assert!(
        bg_p > 0.05,
        "background_panel p-value should be > 0.05, got {bg_p}"
    );

    // Sidecar must exist with input hashes.
    let sidecar = tmp.path().join("gsea.tsv.run.json");
    assert!(sidecar.exists(), "sidecar missing");
    let sidecar_body = std::fs::read_to_string(&sidecar).unwrap();
    assert!(sidecar_body.contains("\"command\""));
    assert!(sidecar_body.contains("enrich gsea"));
    assert!(sidecar_body.contains("\"de_results\""));
    assert!(sidecar_body.contains("\"gene_sets\""));
}

#[test]
fn enrich_gsea_two_runs_same_seed_match_byte_for_byte() {
    let tmp = tempfile::tempdir().unwrap();
    let de = tmp.path().join("de_results.tsv");
    let sets = tmp.path().join("gene_sets.tsv");
    let a = tmp.path().join("a.tsv");
    let b = tmp.path().join("b.tsv");
    let top: Vec<String> = (0..20).map(|i| format!("HIT{:03}", i)).collect();
    let rest: Vec<String> = (0..150).map(|i| format!("BG{:03}", i)).collect();
    let top_refs: Vec<&str> = top.iter().map(String::as_str).collect();
    let rest_refs: Vec<&str> = rest.iter().map(String::as_str).collect();
    write_de_results(&de, &top_refs, &rest_refs);
    let mut s = String::from("set_name\tgene_symbol\n");
    for g in &top[0..12] {
        s.push_str(&format!("planted\t{g}\n"));
    }
    std::fs::write(&sets, s).unwrap();
    for out in [&a, &b] {
        let r = run_atman(&[
            "enrich",
            "gsea",
            "--de-results",
            de.to_str().unwrap(),
            "--gene-sets",
            sets.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--n-permutations",
            "300",
            "--seed",
            "7",
        ]);
        assert!(r.status.success());
    }
    assert_eq!(
        std::fs::read(&a).unwrap(),
        std::fs::read(&b).unwrap(),
        "same-seed runs must produce byte-identical output"
    );
}

/// R-parity test: atman `enrich gsea` ES values must match the
/// fgsea::fgseaSimple reference within tight numerical tolerance on
/// the planted fixture in `tests/fixtures/gsea_reference.tsv`.
/// Reference is regenerated by `gsea_reference.R`; the TSV is checked
/// in so this test runs without R available.
#[test]
fn enrich_gsea_es_matches_fgsea_simple_reference() {
    let fixture_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let de = fixture_dir.join("gsea_de_results.tsv");
    let sets = fixture_dir.join("gsea_gene_sets.tsv");
    let ref_path = fixture_dir.join("gsea_reference.tsv");
    assert!(de.exists(), "fixture missing: {de:?} (run gsea_reference.R)");
    assert!(sets.exists(), "fixture missing: {sets:?}");
    assert!(ref_path.exists(), "fixture missing: {ref_path:?}");

    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("gsea.tsv");
    let r = run_atman(&[
        "enrich",
        "gsea",
        "--de-results",
        de.to_str().unwrap(),
        "--gene-sets",
        sets.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--n-permutations",
        "100",
        "--seed",
        "20260428",
        "--min-set-size",
        "1",
    ]);
    assert!(
        r.status.success(),
        "atman enrich gsea failed: {}",
        String::from_utf8_lossy(&r.stderr)
    );

    let mut atman_es: std::collections::BTreeMap<String, (usize, f64)> = std::collections::BTreeMap::new();
    let body = std::fs::read_to_string(&out).unwrap();
    let mut lines = body.lines();
    let header = lines.next().unwrap();
    let cols: Vec<&str> = header.split('\t').collect();
    let name_i = cols.iter().position(|c| *c == "set_name").unwrap();
    let size_i = cols.iter().position(|c| *c == "set_size").unwrap();
    let es_i = cols.iter().position(|c| *c == "es").unwrap();
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        atman_es.insert(
            f[name_i].to_string(),
            (f[size_i].parse().unwrap(), f[es_i].parse().unwrap()),
        );
    }

    let ref_body = std::fs::read_to_string(&ref_path).unwrap();
    let mut ref_lines = ref_body.lines();
    let _ = ref_lines.next(); // header
    for line in ref_lines {
        let f: Vec<&str> = line.split('\t').collect();
        let name = f[0].to_string();
        let ref_size: usize = f[1].parse().unwrap();
        let ref_es: f64 = f[2].parse().unwrap();
        let (ours_size, ours_es) = atman_es
            .get(&name)
            .unwrap_or_else(|| panic!("atman missing set {name}"));
        assert_eq!(*ours_size, ref_size, "set_size mismatch for {name}");
        assert!(
            (ours_es - ref_es).abs() < 1e-12,
            "ES mismatch for {name}: atman={ours_es} fgsea={ref_es}"
        );
    }
}

#[test]
fn enrich_gsea_requires_comparison_when_multiple_present() {
    let tmp = tempfile::tempdir().unwrap();
    let de = tmp.path().join("de_results.tsv");
    let sets = tmp.path().join("gene_sets.tsv");
    let out = tmp.path().join("gsea.tsv");
    let mut s = String::from(
        "panel\tassay_id\tgene_symbol\tuniprot\tcomparison\tn_pairs\tmean_a\tmean_b\tmean_diff\tt\tdf\tp_value\tbh_q\tskip_reason\n",
    );
    for i in 0..30 {
        let g = format!("G{i:03}");
        s.push_str(&format!(
            "p\t{g}\t{g}\t\tCase-Control\t10\t\t\t1.0\t4.0\t8\t1e-4\t1e-3\t\n"
        ));
        s.push_str(&format!(
            "p\t{g}\t{g}\t\tCase-Other\t10\t\t\t0.5\t2.0\t8\t0.05\t0.1\t\n"
        ));
    }
    std::fs::write(&de, s).unwrap();
    let mut sets_text = String::from("set_name\tgene_symbol\n");
    for i in 0..6 {
        sets_text.push_str(&format!("planted\tG{i:03}\n"));
    }
    std::fs::write(&sets, sets_text).unwrap();

    let r = run_atman(&[
        "enrich",
        "gsea",
        "--de-results",
        de.to_str().unwrap(),
        "--gene-sets",
        sets.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--n-permutations",
        "100",
    ]);
    assert!(
        !r.status.success(),
        "should error when multiple comparisons present without --comparison"
    );
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(
        stderr.contains("--comparison") || stderr.contains("comparison"),
        "stderr should mention --comparison: {stderr}"
    );
}
