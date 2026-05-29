use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

/// Build a 2-sample × 10-gene fixture where:
/// - Sample S1 has a clear "up" signature: signature genes hold the top 3 ranks.
/// - Sample S2 has the same signature placed in the middle (uninformative).
fn write_fixture(dir: &std::path::Path) -> std::path::PathBuf {
    std::fs::create_dir_all(dir).unwrap();
    std::fs::write(
        dir.join("samples.tsv"),
        "\
sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n\
S1\tS1\tCase\t0\tcsf\t1\n\
S2\tS2\tControl\t1\tcsf\t2\n",
    )
    .unwrap();

    let n_genes = 10;
    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for i in 0..n_genes {
        proteins.push_str(&format!(
            "spectronaut_report\tP{:05}\tP{:05}\tG{:02}\tspectronaut\t\n",
            i, i, i
        ));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    let mut measurements = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut idx = 0;
    // S1: G07, G08, G09 have the highest abundance; rest ascending.
    for i in 0..n_genes {
        let abundance = i as f64;
        idx += 1;
        measurements.push_str(&format!(
            "spectronaut_report\tS1\tP{:05}\tG{:02}\tspectronaut\t{a:.3}\t{a:.3}\t{a:.3}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{idx}\n",
            i, i, a = abundance
        ));
    }
    // S2: signature genes G07/G08/G09 land at positions 4, 5, 6 (middle).
    let s2_values: Vec<f64> = vec![
        0.0, // G00 -> rank 1
        1.0, // G01 -> rank 2
        2.0, // G02 -> rank 3
        8.0, // G03 -> rank 9
        9.0, // G04 -> rank 10
        7.0, // G05 -> rank 8
        6.0, // G06 -> rank 7
        3.0, // G07 -> rank 4   (signature)
        4.0, // G08 -> rank 5   (signature)
        5.0, // G09 -> rank 6   (signature)
    ];
    for (i, v) in s2_values.iter().enumerate() {
        idx += 1;
        measurements.push_str(&format!(
            "spectronaut_report\tS2\tP{:05}\tG{:02}\tspectronaut\t{v:.3}\t{v:.3}\t{v:.3}\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{idx}\n",
            i, i
        ));
    }
    std::fs::write(dir.join("measurements.tsv"), measurements).unwrap();

    let sets = dir.join("gene_sets.tsv");
    std::fs::write(
        &sets,
        "\
set_name\tgene_symbol\n\
top_three\tG07\n\
top_three\tG08\n\
top_three\tG09\n\
absent_from_data\tQ01\n\
absent_from_data\tQ02\n\
absent_from_data\tQ03\n",
    )
    .unwrap();
    sets
}

#[test]
fn score_signatures_singscore_separates_top_loaded_from_middle() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let out = tmp.path().join("signature_scores.tsv");
    let sets = write_fixture(&input);

    let r = run_atman(&[
        "score",
        "signatures",
        "--input-dir",
        input.to_str().unwrap(),
        "--gene-sets",
        sets.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--min-set-size",
        "1",
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );

    let body = std::fs::read_to_string(&out).unwrap();
    let mut lines = body.lines();
    let header = lines.next().unwrap();
    assert_eq!(
        header,
        "sample_id\tsubject_id\tcondition\tset_name\tmethod\tscore\tn_genes_declared\tn_genes_observed_in_sample\tn_proteins_in_sample"
    );
    let rows: Vec<Vec<&str>> = lines.map(|l| l.split('\t').collect()).collect();

    // Locate rows by (set_name, sample_id).
    let by_key = |set: &str, sid: &str| -> &Vec<&str> {
        rows.iter()
            .find(|r| r[3] == set && r[0] == sid)
            .unwrap_or_else(|| panic!("row missing for set={set} sample={sid}"))
    };

    // S1 / top_three: G07/G08/G09 are top 3 → score = 0.5 (perfect top).
    let s1_top = by_key("top_three", "S1");
    let s1_score: f64 = s1_top[5].parse().unwrap();
    assert!(
        (s1_score - 0.5).abs() < 1e-12,
        "S1 top_three score should be 0.5, got {s1_score}"
    );
    assert_eq!(s1_top[6], "3", "n_genes_declared");
    assert_eq!(s1_top[7], "3", "n_genes_observed_in_sample");
    assert_eq!(s1_top[8], "10", "n_proteins_in_sample");

    // S2 / top_three: signature lies at ranks 4,5,6 (mean rank 5 of 10).
    // score = (2*5 - 11) / (2*7) = -1/14 ≈ -0.0714
    let s2_top = by_key("top_three", "S2");
    let s2_score: f64 = s2_top[5].parse().unwrap();
    assert!(
        (s2_score - (-1.0 / 14.0)).abs() < 1e-12,
        "S2 top_three score should be -1/14 ≈ -0.0714, got {s2_score}"
    );

    // absent_from_data: 0 of 3 signature genes observed → NaN rows.
    for sid in ["S1", "S2"] {
        let r = by_key("absent_from_data", sid);
        let v = r[5];
        assert!(
            v == "NaN" || v == "nan",
            "absent_from_data score should be NaN for {sid}, got {v}"
        );
        assert_eq!(r[7], "0", "n_genes_observed_in_sample for absent set");
    }

    // Sidecar must exist.
    let sidecar = tmp.path().join("signature_scores.tsv.run.json");
    assert!(sidecar.exists(), "sidecar missing");
    let sidecar_body = std::fs::read_to_string(&sidecar).unwrap();
    assert!(sidecar_body.contains("score signatures"));
    assert!(sidecar_body.contains("\"singscore\""));
    assert!(sidecar_body.contains("\"gene_sets\""));
}

/// R-parity test: atman `score signatures --method singscore` TotalScores
/// must match the singscore::simpleScore reference within tight numerical
/// tolerance on the planted fixture in
/// `tests/fixtures/singscore_reference.tsv`. Reference is regenerated by
/// `singscore_reference.R`; the TSV is checked in so this test runs
/// without R available.
#[test]
fn score_signatures_singscore_matches_reference() {
    let fixture_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let samples = fixture_dir.join("singscore_samples.tsv");
    let proteins = fixture_dir.join("singscore_proteins.tsv");
    let measurements = fixture_dir.join("singscore_measurements.tsv");
    let sets = fixture_dir.join("singscore_gene_sets.tsv");
    let reference = fixture_dir.join("singscore_reference.tsv");
    for f in [&samples, &proteins, &measurements, &sets, &reference] {
        assert!(
            f.exists(),
            "fixture missing: {f:?} (run singscore_reference.R)"
        );
    }

    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    std::fs::create_dir_all(&input).unwrap();
    std::fs::copy(&samples, input.join("samples.tsv")).unwrap();
    std::fs::copy(&proteins, input.join("proteins.tsv")).unwrap();
    std::fs::copy(&measurements, input.join("measurements.tsv")).unwrap();
    let out = tmp.path().join("scores.tsv");

    let r = run_atman(&[
        "score",
        "signatures",
        "--input-dir",
        input.to_str().unwrap(),
        "--gene-sets",
        sets.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--min-set-size",
        "1",
    ]);
    assert!(
        r.status.success(),
        "atman score signatures failed: {}",
        String::from_utf8_lossy(&r.stderr)
    );

    let mut atman: std::collections::BTreeMap<(String, String), f64> =
        std::collections::BTreeMap::new();
    let body = std::fs::read_to_string(&out).unwrap();
    let mut lines = body.lines();
    let header = lines.next().unwrap();
    let cols: Vec<&str> = header.split('\t').collect();
    let set_i = cols.iter().position(|c| *c == "set_name").unwrap();
    let sid_i = cols.iter().position(|c| *c == "sample_id").unwrap();
    let score_i = cols.iter().position(|c| *c == "score").unwrap();
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        atman.insert(
            (f[set_i].to_string(), f[sid_i].to_string()),
            f[score_i].trim().parse().unwrap(),
        );
    }

    let ref_body = std::fs::read_to_string(&reference).unwrap();
    let mut ref_lines = ref_body.lines();
    let _ = ref_lines.next();
    let mut compared = 0usize;
    for line in ref_lines {
        let f: Vec<&str> = line.split('\t').collect();
        let key = (f[0].to_string(), f[1].to_string());
        let ref_score: f64 = f[2].trim().parse().unwrap();
        let ours = atman
            .get(&key)
            .unwrap_or_else(|| panic!("atman missing {key:?}"));
        assert!(
            (ours - ref_score).abs() < 1e-12,
            "score mismatch for {key:?}: atman={ours} singscore={ref_score}"
        );
        compared += 1;
    }
    assert_eq!(
        compared, 12,
        "expected 12 reference rows (3 sets × 4 samples), saw {compared}"
    );
}

#[test]
fn score_signatures_two_runs_byte_identical() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let out_a = tmp.path().join("a.tsv");
    let out_b = tmp.path().join("b.tsv");
    let sets = write_fixture(&input);
    for out in [&out_a, &out_b] {
        let r = run_atman(&[
            "score",
            "signatures",
            "--input-dir",
            input.to_str().unwrap(),
            "--gene-sets",
            sets.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--min-set-size",
            "1",
        ]);
        assert!(r.status.success());
    }
    assert_eq!(
        std::fs::read(&out_a).unwrap(),
        std::fs::read(&out_b).unwrap(),
        "two runs must produce byte-identical output (singscore is deterministic)"
    );
}
