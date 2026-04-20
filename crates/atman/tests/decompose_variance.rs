//! Integration test for `atman decompose variance`.
//!
//! Builds a synthetic activations table with three archetypes:
//! - `cohort_shared` — activations driven by a per-sample random
//!   value (no cohort / condition effect).
//! - `cohort_confounded` — activations driven entirely by cohort
//!   label.
//! - `condition_specific` — activations driven by the condition
//!   label.
//!
//! Runs `atman decompose variance --factors "cohort + condition"`
//! and asserts that each archetype's dominant variance source
//! matches its planted role.

use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn parse_tsv(
    path: &Path,
) -> (Vec<String>, Vec<std::collections::HashMap<String, String>>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines.next().unwrap().split('\t').map(String::from).collect();
    let rows = lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(h, c)| (h.clone(), c.to_string()))
                .collect()
        })
        .collect();
    (header, rows)
}

/// Simple LCG independent of atman's Xoshiro.
struct Lcg(u64);
impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as f64 / (u32::MAX as f64)
    }
}

fn write_variance_fixture(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    // 24 samples: 12 per cohort (A, B), 6 per condition (X, Y).
    // Canonical samples.tsv columns (read_samples parses by position)
    // followed by an extra `cohort` column that the formula parser
    // picks up via the samples-extras reader.
    let mut samples = String::from(
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tcohort\n",
    );
    let mut pairs = Vec::new();
    let mut order = 0usize;
    for cohort in ["A", "B"] {
        for condition in ["X", "Y"] {
            for i in 0..6 {
                order += 1;
                let sid = format!("{cohort}{condition}_{i:02}");
                pairs.push((sid.clone(), cohort.to_string(), condition.to_string()));
                samples.push_str(&format!(
                    "{sid}\t{sid}\t{condition}\t0\tplasma\t{order}\t{cohort}\n"
                ));
            }
        }
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    // activations.tsv per atman decompose ica contract:
    // cohort\tsubject_id\tsample_id\tprogram\tactivation
    let mut rng = Lcg::new(20260418);
    let mut act = String::from("cohort\tsubject_id\tsample_id\tprogram\tactivation\n");
    for (sid, cohort, condition) in &pairs {
        for program in ["cohort_shared", "cohort_confounded", "condition_specific"] {
            let value = match program {
                "cohort_shared" => (rng.next() - 0.5) * 2.0, // pure noise
                "cohort_confounded" => {
                    let base = if cohort == "A" { 2.0 } else { -2.0 };
                    base + (rng.next() - 0.5) * 0.2
                }
                "condition_specific" => {
                    let base = if condition == "X" { -1.5 } else { 1.5 };
                    base + (rng.next() - 0.5) * 0.2
                }
                _ => 0.0,
            };
            act.push_str(&format!(
                "bench\t{sid}\t{sid}\t{program}\t{value:.6}\n"
            ));
        }
    }
    std::fs::write(dir.join("activations.tsv"), act).unwrap();
}

#[test]
fn decompose_variance_attributes_variance_to_planted_factor() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = tmp.path().join("fixture");
    write_variance_fixture(&fixture);
    let out = tmp.path().join("variance.tsv");

    let status = run_atman(&[
        "decompose",
        "variance",
        "--activations",
        fixture.join("activations.tsv").to_str().unwrap(),
        "--samples",
        fixture.join("samples.tsv").to_str().unwrap(),
        "--factors",
        "cohort + condition",
        "--min-samples",
        "4",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "decompose variance failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&status.stdout),
        String::from_utf8_lossy(&status.stderr)
    );

    let (header, rows) = parse_tsv(&out);
    for col in [
        "archetype_id",
        "total_var",
        "var_residual",
        "var_cohort",
        "var_condition",
        "max_abs_t_cohort",
        "max_abs_t_condition",
    ] {
        assert!(header.iter().any(|h| h == col), "missing column {col}");
    }
    let by_id: std::collections::HashMap<&str, &std::collections::HashMap<String, String>> = rows
        .iter()
        .map(|r| (r["archetype_id"].as_str(), r))
        .collect();

    // cohort_confounded: most variance should be in cohort, not condition.
    let cc = by_id["cohort_confounded"];
    let var_cohort: f64 = cc["var_cohort"].parse().unwrap();
    let var_condition: f64 = cc["var_condition"].parse().unwrap();
    let var_residual: f64 = cc["var_residual"].parse().unwrap();
    assert!(
        var_cohort > var_condition,
        "cohort_confounded: var_cohort={var_cohort}, var_condition={var_condition}"
    );
    assert!(
        var_cohort > var_residual,
        "cohort_confounded: var_cohort={var_cohort}, var_residual={var_residual}"
    );

    // condition_specific: most variance should be in condition, not cohort.
    let cs = by_id["condition_specific"];
    let var_cohort: f64 = cs["var_cohort"].parse().unwrap();
    let var_condition: f64 = cs["var_condition"].parse().unwrap();
    let var_residual: f64 = cs["var_residual"].parse().unwrap();
    assert!(
        var_condition > var_cohort,
        "condition_specific: var_condition={var_condition}, var_cohort={var_cohort}"
    );
    assert!(
        var_condition > var_residual,
        "condition_specific: var_condition={var_condition}, var_residual={var_residual}"
    );

    // cohort_shared (pure noise): neither factor dominates the
    // residual at the level we see in the planted archetypes.
    let sh = by_id["cohort_shared"];
    let var_cohort: f64 = sh["var_cohort"].parse().unwrap();
    let var_condition: f64 = sh["var_condition"].parse().unwrap();
    let var_residual: f64 = sh["var_residual"].parse().unwrap();
    assert!(
        var_cohort < var_residual,
        "cohort_shared: var_cohort={var_cohort}, var_residual={var_residual}"
    );
    assert!(
        var_condition < var_residual,
        "cohort_shared: var_condition={var_condition}, var_residual={var_residual}"
    );

    // Sidecar.
    let sidecar = out.parent().unwrap().join("variance.tsv.run.json");
    assert!(sidecar.exists());
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();
    assert_eq!(j["args"]["factors"], "cohort + condition");
}

#[test]
fn decompose_variance_supports_random_intercept_group() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = tmp.path().join("fixture");
    write_variance_fixture(&fixture);
    let out = tmp.path().join("variance_mixed.tsv");

    let status = run_atman(&[
        "decompose",
        "variance",
        "--activations",
        fixture.join("activations.tsv").to_str().unwrap(),
        "--samples",
        fixture.join("samples.tsv").to_str().unwrap(),
        "--factors",
        "condition + (1|cohort)",
        "--min-samples",
        "4",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&status.stderr)
    );
    let (header, rows) = parse_tsv(&out);
    assert!(header.iter().any(|h| h == "var_random_cohort"));
    assert!(header.iter().any(|h| h == "icc_random_cohort"));
    // cohort_confounded under cohort as random intercept should show
    // non-zero ICC.
    let cc = rows
        .iter()
        .find(|r| r["archetype_id"] == "cohort_confounded")
        .unwrap();
    let icc: f64 = cc["icc_random_cohort"].parse().unwrap();
    assert!(
        icc > 0.3,
        "cohort_confounded ICC should be substantial under cohort grouping; got {icc}"
    );
}

#[test]
fn decompose_variance_refuses_single_level_factor() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = tmp.path().join("fixture");
    write_variance_fixture(&fixture);
    let out = tmp.path().join("variance_bad.tsv");

    // sample_type is "plasma" for every sample — one level, no variance.
    let status = run_atman(&[
        "decompose",
        "variance",
        "--activations",
        fixture.join("activations.tsv").to_str().unwrap(),
        "--samples",
        fixture.join("samples.tsv").to_str().unwrap(),
        "--factors",
        "sample_type",
        "--min-samples",
        "4",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(!status.status.success());
    let stderr = String::from_utf8_lossy(&status.stderr);
    assert!(
        stderr.contains("only one level"),
        "expected single-level error; got: {stderr}"
    );
}

#[test]
fn decompose_variance_is_deterministic() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = tmp.path().join("fixture");
    write_variance_fixture(&fixture);
    let out1 = tmp.path().join("r1.tsv");
    let out2 = tmp.path().join("r2.tsv");
    for out in [&out1, &out2] {
        let status = run_atman(&[
            "decompose",
            "variance",
            "--activations",
            fixture.join("activations.tsv").to_str().unwrap(),
            "--samples",
            fixture.join("samples.tsv").to_str().unwrap(),
            "--factors",
            "cohort + condition",
            "--output",
            out.to_str().unwrap(),
        ]);
        assert!(status.status.success());
    }
    let a = std::fs::read_to_string(&out1).unwrap();
    let b = std::fs::read_to_string(&out2).unwrap();
    assert_eq!(a, b);
}
