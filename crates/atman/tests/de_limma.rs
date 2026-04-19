use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_minimal_canonical(dir: &std::path::Path) {
    std::fs::create_dir_all(dir).unwrap();
    // Two groups, five subjects each, two features.
    let samples = "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n".to_string()
        + &(1..=10)
            .map(|i| {
                let cond = if i <= 5 { "A" } else { "B" };
                format!("S{i:02}\tS{i:02}\t{cond}\t0\tplasma\t{i}\n")
            })
            .collect::<String>();
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();
    let proteins = "platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n\
                    olink_explore_ngs\tA001\tQ00001\tGENE1\tP1\t\n\
                    olink_explore_ngs\tA002\tQ00002\tGENE2\tP1\t\n";
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();
    let mut meas = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order = 0;
    for i in 1..=10 {
        for a in ["A001", "A002"] {
            order += 1;
            let effect = if a == "A001" { 1.0 } else { 0.0 };
            let shift = if i <= 5 { 0.0 } else { effect };
            let noise = ((i * 7 + order) as f64).sin() * 0.1;
            let v = shift + noise;
            meas.push_str(&format!(
                "olink_explore_ngs\tS{i:02}\t{a}\tGENE{a_suffix}\tP1\t{v:.6}\t\
                 {v:.6}\t{v:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n",
                a_suffix = &a[1..],
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &meas).unwrap();
    std::fs::write(dir.join("qc_measurements.tsv"), &meas).unwrap();
}

#[test]
fn atman_de_test_limma_smoke() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("de_out");
    write_minimal_canonical(&input);
    std::fs::create_dir_all(&output).unwrap();

    let out = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "limma",
        "--groups",
        "A-B",
        "--min-pairs",
        "3",
        "--trend",
        "false",
        "--robust",
        "false",
    ]);
    assert!(
        out.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let results = std::fs::read_to_string(output.join("de_results.tsv")).unwrap();
    assert!(results.lines().next().unwrap().contains("s2_posterior"));
    assert!(results.lines().next().unwrap().contains("lfc_threshold"));
    // Body rows should have a populated s2_posterior (non-empty column).
    let body = results.lines().nth(1).unwrap();
    let cols: Vec<&str> = body.split('\t').collect();
    let header: Vec<&str> = results.lines().next().unwrap().split('\t').collect();
    let idx = header.iter().position(|c| *c == "s2_posterior").unwrap();
    assert!(!cols[idx].is_empty());

    // Sidecar.
    let sidecar = output.join("de_results.tsv.run.json");
    assert!(sidecar.exists(), "sidecar missing");
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();
    assert_eq!(j["args"]["test"], "limma");
    assert_eq!(j["args"]["trend"], false);
    assert_eq!(j["args"]["robust"], false);
    assert_eq!(j["args"]["lfc-threshold"], 0.0);
}

fn read_fixture_tsv(path: &std::path::Path) -> Vec<std::collections::HashMap<String, String>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let mut rows = Vec::new();
    for line in lines {
        let cells: Vec<&str> = line.split('\t').collect();
        let mut m = std::collections::HashMap::new();
        for (h, c) in header.iter().zip(cells.iter()) {
            m.insert(h.to_string(), c.to_string());
        }
        rows.push(m);
    }
    rows
}

fn copy_fixture_canonical(dst: &std::path::Path) {
    let src = std::path::Path::new("tests/fixtures");
    std::fs::create_dir_all(dst).unwrap();
    for name in [
        "limma_fixture_samples.tsv",
        "limma_fixture_proteins.tsv",
        "limma_fixture_measurements.tsv",
    ] {
        let s = std::fs::read_to_string(src.join(name)).unwrap();
        let dst_name = match name {
            "limma_fixture_samples.tsv" => "samples.tsv",
            "limma_fixture_proteins.tsv" => "proteins.tsv",
            "limma_fixture_measurements.tsv" => "measurements.tsv",
            _ => unreachable!(),
        };
        std::fs::write(dst.join(dst_name), &s).unwrap();
        std::fs::write(dst.join("qc_measurements.tsv"), &s).unwrap(); // mirror
    }
}

fn run_limma_and_load(trend: bool) -> Vec<std::collections::HashMap<String, String>> {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("canonical");
    let output = tmp.path().join("de_out");
    copy_fixture_canonical(&input);
    std::fs::create_dir_all(&output).unwrap();
    let trend_arg = if trend { "true" } else { "false" };
    let out = run_atman(&[
        "de",
        "--input-dir",
        input.to_str().unwrap(),
        "--output-dir",
        output.to_str().unwrap(),
        "--test",
        "limma",
        "--groups",
        "A-B",
        "--min-pairs",
        "3",
        "--trend",
        trend_arg,
        "--robust",
        "false",
    ]);
    assert!(
        out.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    read_fixture_tsv(&output.join("de_results.tsv"))
}

fn abs_diff(a: &str, b: &str) -> f64 {
    let a: f64 = a.parse().unwrap_or(f64::NAN);
    let b: f64 = b.parse().unwrap_or(f64::NAN);
    (a - b).abs()
}

fn index_by_assay(
    rows: &[std::collections::HashMap<String, String>],
    key_col: &str,
) -> std::collections::HashMap<String, std::collections::HashMap<String, String>> {
    rows.iter()
        .map(|row| (row.get(key_col).unwrap().clone(), row.clone()))
        .collect()
}

#[test]
fn matches_r_limma_notrend_within_tolerance() {
    let rust_rows = run_limma_and_load(false);
    let ref_rows = read_fixture_tsv(std::path::Path::new(
        "tests/fixtures/limma_reference_notrend.tsv",
    ));
    // Rust sorts rows by (comparison, bh_q, |mean_diff|, gene_symbol);
    // the R fixture is in F001-F100 input order. Match by assay_id, not zip.
    let ref_by_assay = index_by_assay(&ref_rows, "feature");
    assert_eq!(rust_rows.len(), ref_rows.len());
    for rust in &rust_rows {
        let assay = rust.get("assay_id").unwrap();
        let r = ref_by_assay
            .get(assay)
            .unwrap_or_else(|| panic!("no R ref for {assay}"));
        // Rust reports mean_a - mean_b (contrast [0,-1]); R records (mean_b - mean_a)/se. Compare by magnitude.
        let rust_t_abs: f64 = rust.get("t").unwrap().parse::<f64>().unwrap().abs();
        let r_t_abs: f64 = r.get("t").unwrap().parse::<f64>().unwrap().abs();
        assert!(
            (rust_t_abs - r_t_abs).abs() < 1e-4,
            "assay={assay} |t| rust={rust_t_abs} r={r_t_abs}"
        );
        assert!(
            abs_diff(rust.get("p_value").unwrap(), r.get("P_Value").unwrap()) < 1e-4,
            "assay={assay} p_value diverges"
        );
        assert!(
            abs_diff(rust.get("s2_posterior").unwrap(), r.get("s2_posterior").unwrap())
                < 1e-4,
            "assay={assay} s2_posterior diverges"
        );
        assert!(
            abs_diff(rust.get("df_total").unwrap(), r.get("df_total").unwrap()) < 1e-4,
            "assay={assay} df_total diverges"
        );
    }
}

#[test]
fn trend_mode_populates_expected_columns_with_finite_stats() {
    // Sanity check, not a numerical-parity test. Trend-mode parity against
    // limma::eBayes(trend=TRUE) is deferred per spec commit aac8fd6: atman's
    // parametric-quadratic trend + non-trend-eBayes-on-ratios is algorithmically
    // different from limma's covariate-spline path in fitFDist(covariate=...),
    // and the two are expected to diverge by up to ~10% on |t|/s2_posterior.
    // Tight parity would require porting splines + effects-based df estimator
    // + per-feature s20.
    let rust_rows = run_limma_and_load(true);
    let ref_rows = read_fixture_tsv(std::path::Path::new(
        "tests/fixtures/limma_reference_trend.tsv",
    ));
    assert_eq!(rust_rows.len(), ref_rows.len());
    assert_eq!(rust_rows.len(), 100);
    for rust in &rust_rows {
        let assay = rust.get("assay_id").unwrap();

        let s2_trend_str = rust.get("s2_trend").expect("s2_trend column present");
        assert!(
            !s2_trend_str.is_empty(),
            "assay={assay} s2_trend unpopulated"
        );
        let s2_trend: f64 = s2_trend_str
            .parse()
            .unwrap_or_else(|_| panic!("assay={assay} s2_trend not parseable: {s2_trend_str}"));
        assert!(
            s2_trend.is_finite() && s2_trend > 0.0,
            "assay={assay} s2_trend not finite/positive: {s2_trend}"
        );

        let s2_post_str = rust
            .get("s2_posterior")
            .expect("s2_posterior column present");
        assert!(
            !s2_post_str.is_empty(),
            "assay={assay} s2_posterior unpopulated"
        );
        let s2_post: f64 = s2_post_str
            .parse()
            .unwrap_or_else(|_| panic!("assay={assay} s2_posterior not parseable: {s2_post_str}"));
        assert!(
            s2_post.is_finite() && s2_post > 0.0,
            "assay={assay} s2_posterior not finite/positive: {s2_post}"
        );

        let t_str = rust.get("t").expect("t column present");
        let t: f64 = t_str
            .parse()
            .unwrap_or_else(|_| panic!("assay={assay} t not parseable: {t_str}"));
        assert!(t.is_finite(), "assay={assay} t not finite: {t}");

        let skip_reason = rust.get("skip_reason").map(String::as_str).unwrap_or("");
        let p_str = rust.get("p_value").expect("p_value column present");
        if skip_reason.is_empty() {
            let p: f64 = p_str
                .parse()
                .unwrap_or_else(|_| panic!("assay={assay} p_value not parseable: {p_str}"));
            assert!(
                p.is_finite() && (0.0..=1.0).contains(&p),
                "assay={assay} p_value out of [0,1]: {p}"
            );
        }

        let df_str = rust.get("df_total").expect("df_total column present");
        let df: f64 = df_str
            .parse()
            .unwrap_or_else(|_| panic!("assay={assay} df_total not parseable: {df_str}"));
        assert!(
            df.is_finite() && df > 0.0,
            "assay={assay} df_total not finite/positive: {df}"
        );
    }
}
