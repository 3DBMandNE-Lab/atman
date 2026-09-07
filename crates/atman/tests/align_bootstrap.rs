//! Integration test for `atman align bootstrap`.
//!
//! Two synthetic cohorts with a planted universal archetype (same
//! loading pattern in both cohorts) and cohort-specific archetypes
//! (A has its own planted signal on proteins 4..6; B has its own on
//! proteins 7..9). Runs subject-level bootstrap and asserts:
//!
//! 1. Output schema is as specified.
//! 2. At least one point-estimate archetype is detected.
//! 3. Every archetype row has finite, bounded statistics.
//! 4. Determinism across runs under fixed --seed.

use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

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
    fn heavy_tail(&mut self) -> f64 {
        let u = self.next() * 2.0 - 1.0;
        u * u * u
    }
}

fn write_cohort(
    dir: &Path,
    cohort_label: &str,
    planted_a_specific: bool,
    planted_b_specific: bool,
    seed: u64,
) {
    std::fs::create_dir_all(dir).unwrap();
    let n_samples = 30usize;
    let n_proteins = 10usize;

    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=n_samples {
        samples.push_str(&format!(
            "{cohort_label}_S{i:03}\t{cohort_label}_S{i:03}\tN/A\t0\tplasma\t{i}\n"
        ));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 1..=n_proteins {
        proteins.push_str(&format!(
            "olink_explore_ngs\tA{j:03}\tQ{j:05}\tG{j:03}\tP1\t\n"
        ));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    // Universal archetype loadings on proteins 1..4; cohort-specific
    // loadings on 5..7 (A) or 8..10 (B).
    let mut rng = Lcg::new(seed);
    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order: u64 = 0;
    for i in 1..=n_samples {
        // Per-sample heavy-tailed source drives each archetype.
        let universal_source = rng.heavy_tail() * 2.0;
        let a_source = rng.heavy_tail() * 2.0;
        let b_source = rng.heavy_tail() * 2.0;
        for j in 1..=n_proteins {
            order += 1;
            let universal_weight = if (1..=4).contains(&j) { 1.0 } else { 0.0 };
            let a_weight = if planted_a_specific && (5..=7).contains(&j) {
                1.0
            } else {
                0.0
            };
            let b_weight = if planted_b_specific && (8..=10).contains(&j) {
                1.0
            } else {
                0.0
            };
            let noise = (rng.next() - 0.5) * 0.3;
            let value = universal_weight * universal_source
                + a_weight * a_source
                + b_weight * b_source
                + noise;
            qc.push_str(&format!(
                "olink_explore_ngs\t{cohort_label}_S{i:03}\tA{j:03}\tG{j:03}\tP1\t{value:.6}\t\
                 {value:.6}\t{value:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
}

/// A cohort with no planted structure at all.
///
/// The planted fixture shares a strong universal archetype across both
/// cohorts, so every resample recovers it and the match rate sits near
/// 1.0. That is the wrong regime for testing the low-match-rate
/// diagnostic: the warning never fires and an iff assertion against it
/// passes trivially. Noise gives a match rate in the range the warning
/// exists for.
fn write_noise_cohort(dir: &Path, cohort_label: &str, seed: u64) {
    std::fs::create_dir_all(dir).unwrap();
    let n_samples = 30usize;
    let n_proteins = 12usize;

    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=n_samples {
        samples.push_str(&format!(
            "{cohort_label}_S{i:03}\t{cohort_label}_S{i:03}\tN/A\t0\tplasma\t{i}\n"
        ));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 1..=n_proteins {
        proteins.push_str(&format!(
            "olink_explore_ngs\tA{j:03}\tQ{j:05}\tG{j:03}\tP1\t\n"
        ));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    let mut rng = Lcg::new(seed);
    let mut m = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order: u64 = 0;
    for i in 1..=n_samples {
        for j in 1..=n_proteins {
            order += 1;
            let a = 10.0 + rng.heavy_tail() * 2.0;
            m.push_str(&format!(
                "olink_explore_ngs\t{cohort_label}_S{i:03}\tA{j:03}\tG{j:03}\tP1\t\
                 {a:.4}\t{a:.4}\t{a:.4}\tnpx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), m).unwrap();
}

fn parse_tsv(path: &Path) -> (Vec<String>, Vec<std::collections::HashMap<String, String>>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect();
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

#[test]
fn align_bootstrap_runs_on_two_cohort_synthetic_fixture() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_cohort(&a, "A", true, false, 11);
    write_cohort(&b, "B", false, true, 22);
    let out = tmp.path().join("bootstrap").join("summary.tsv");

    let status = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        &format!("{},{}", a.display(), b.display()),
        "--labels",
        "A,B",
        "--k",
        "2",
        "--n-boot",
        "10",
        "--seed",
        "20260418",
        "--cosine-tau",
        "0.1",
        "--match-tau",
        "0.1",
        "--min-subjects",
        "10",
        "--max-iter",
        "80",
        "--tol",
        "1e-3",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "align bootstrap failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&status.stdout),
        String::from_utf8_lossy(&status.stderr)
    );

    let (header, rows) = parse_tsv(&out);
    for col in [
        "archetype_id",
        "observed_n_cohorts",
        "observed_cohorts",
        "bootstrap_mean_n_cohorts",
        "bootstrap_prob_universal",
        "bootstrap_prob_multi",
        "ci_lower_n_cohorts",
        "ci_upper_n_cohorts",
        "bootstrap_match_rate",
        "alignment_entropy",
        "bca_lower_n_cohorts",
        "bca_upper_n_cohorts",
        "bca_fallback_to_percentile",
    ] {
        assert!(header.iter().any(|h| h == col), "missing column {col}");
    }
    for r in &rows {
        let n: usize = r["observed_n_cohorts"].parse().unwrap();
        assert!(
            (2..=2).contains(&n),
            "expected observed_n_cohorts in [2,2] for a universal archetype; got {n}"
        );
        let mean: f64 = r["bootstrap_mean_n_cohorts"].parse().unwrap();
        assert!(mean >= 0.0);
        let prob_u: f64 = r["bootstrap_prob_universal"].parse().unwrap();
        assert!((0.0..=1.0 + 1e-9).contains(&prob_u));
        let prob_m: f64 = r["bootstrap_prob_multi"].parse().unwrap();
        assert!((0.0..=1.0 + 1e-9).contains(&prob_m));
        let mr: f64 = r["bootstrap_match_rate"].parse().unwrap();
        assert!((0.0..=1.0 + 1e-9).contains(&mr));
        let entropy: f64 = r["alignment_entropy"].parse().unwrap();
        assert!(entropy.is_finite() && entropy >= 0.0);
        let bca_lo: f64 = r["bca_lower_n_cohorts"].parse().unwrap();
        let bca_hi: f64 = r["bca_upper_n_cohorts"].parse().unwrap();
        assert!(bca_lo.is_finite() && bca_hi.is_finite());
        assert!(bca_lo <= bca_hi, "BCa CI inverted: [{bca_lo}, {bca_hi}]");
    }

    // Sidecar shape.
    let sidecar_path = out.parent().unwrap().join("summary.tsv.run.json");
    assert!(sidecar_path.exists());
    let j: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&sidecar_path).unwrap()).unwrap();
    assert_eq!(j["args"]["k"], 2);
    assert_eq!(j["args"]["n-boot"], 10);
    assert_eq!(j["args"]["labels"], "A,B");
}

#[test]
fn align_bootstrap_refuses_single_cohort() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    write_cohort(&a, "A", true, false, 1);
    let out = tmp.path().join("summary.tsv");
    let status = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        a.to_str().unwrap(),
        "--k",
        "2",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(!status.status.success());
    let stderr = String::from_utf8_lossy(&status.stderr);
    assert!(stderr.contains("at least 2"), "unexpected: {stderr}");
}

#[test]
fn align_bootstrap_refuses_small_cohorts() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_cohort(&a, "A", true, false, 1);
    write_cohort(&b, "B", false, true, 2);
    let out = tmp.path().join("summary.tsv");
    let status = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        &format!("{},{}", a.display(), b.display()),
        "--k",
        "2",
        "--n-boot",
        "2",
        "--min-subjects",
        "1000",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(!status.status.success());
    let stderr = String::from_utf8_lossy(&status.stderr);
    assert!(stderr.contains("min-subjects"), "unexpected: {stderr}");
}

#[test]
fn align_bootstrap_magnitude_universal_vs_specific_on_strong_fixture() {
    // Stronger-SNR version of the planted fixture: 40 subjects per
    // cohort, n_boot=30. The *universal* archetype (loaded on the
    // shared proteins in both cohorts) should have high
    // bootstrap_prob_multi and low alignment_entropy; its BCa CI
    // should include 2. This is the magnitude contract that DEBT-4
    // asserts — "the bootstrap has power to distinguish planted
    // universal from noise," not just "the numbers are finite."
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_strong_cohort(&a, "A", true, false, 101);
    write_strong_cohort(&b, "B", false, true, 202);
    let out = tmp.path().join("summary.tsv");
    let status = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        &format!("{},{}", a.display(), b.display()),
        "--labels",
        "A,B",
        "--k",
        "2",
        "--n-boot",
        "30",
        "--seed",
        "20260420",
        "--cosine-tau",
        "0.1",
        "--match-tau",
        "0.1",
        "--min-subjects",
        "10",
        "--max-iter",
        "80",
        "--tol",
        "1e-3",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "align bootstrap failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&status.stdout),
        String::from_utf8_lossy(&status.stderr)
    );
    let (_, rows) = parse_tsv(&out);
    assert!(
        !rows.is_empty(),
        "no archetypes recovered; fixture is broken"
    );
    // At least one archetype should behave like a universal archetype:
    // prob_multi >= 0.7, entropy reasonably low, BCa upper ≥ 2.
    let has_universal = rows.iter().any(|r| {
        let prob_m: f64 = r["bootstrap_prob_multi"].parse().unwrap();
        let entropy: f64 = r["alignment_entropy"].parse().unwrap();
        let bca_hi: f64 = r["bca_upper_n_cohorts"].parse().unwrap();
        prob_m >= 0.7 && entropy < 1.0 && bca_hi >= 2.0
    });
    assert!(
        has_universal,
        "expected at least one universal archetype with prob_multi>=0.7 \
         and entropy<1.0 and bca_upper>=2; got rows:\n{}",
        rows.iter()
            .map(|r| format!(
                "  id={} prob_multi={} entropy={} bca=[{}, {}]",
                r["archetype_id"],
                r["bootstrap_prob_multi"],
                r["alignment_entropy"],
                r["bca_lower_n_cohorts"],
                r["bca_upper_n_cohorts"]
            ))
            .collect::<Vec<_>>()
            .join("\n")
    );
}

fn write_strong_cohort(
    dir: &Path,
    cohort_label: &str,
    planted_a_specific: bool,
    planted_b_specific: bool,
    seed: u64,
) {
    // 40 subjects × 15 proteins. Proteins 1-6 are universal; 7-10 are
    // A-specific; 11-15 are B-specific. Noise is low relative to the
    // planted sources so the universal archetype is unambiguous.
    std::fs::create_dir_all(dir).unwrap();
    let n_samples = 40usize;
    let n_proteins = 15usize;
    let mut samples =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for i in 1..=n_samples {
        samples.push_str(&format!(
            "{cohort_label}_S{i:03}\t{cohort_label}_S{i:03}\tN/A\t0\tplasma\t{i}\n"
        ));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();
    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 1..=n_proteins {
        proteins.push_str(&format!(
            "olink_explore_ngs\tA{j:03}\tQ{j:05}\tG{j:03}\tP1\t\n"
        ));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();
    let mut rng = Lcg::new(seed);
    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order: u64 = 0;
    for i in 1..=n_samples {
        let universal_source = rng.heavy_tail() * 3.0;
        let a_source = rng.heavy_tail() * 3.0;
        let b_source = rng.heavy_tail() * 3.0;
        for j in 1..=n_proteins {
            order += 1;
            let universal_weight = if (1..=6).contains(&j) { 1.0 } else { 0.0 };
            let a_weight = if planted_a_specific && (7..=10).contains(&j) {
                1.0
            } else {
                0.0
            };
            let b_weight = if planted_b_specific && (11..=15).contains(&j) {
                1.0
            } else {
                0.0
            };
            let noise = (rng.next() - 0.5) * 0.1;
            let value = universal_weight * universal_source
                + a_weight * a_source
                + b_weight * b_source
                + noise;
            qc.push_str(&format!(
                "olink_explore_ngs\t{cohort_label}_S{i:03}\tA{j:03}\tG{j:03}\tP1\t{value:.6}\t\
                 {value:.6}\t{value:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
}

#[test]
fn align_bootstrap_is_deterministic_across_runs() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_cohort(&a, "A", true, false, 11);
    write_cohort(&b, "B", false, true, 22);
    let out1 = tmp.path().join("r1.tsv");
    let out2 = tmp.path().join("r2.tsv");
    for out in [&out1, &out2] {
        let status = run_atman(&[
            "align",
            "bootstrap",
            "--cohorts",
            &format!("{},{}", a.display(), b.display()),
            "--labels",
            "A,B",
            "--k",
            "2",
            "--n-boot",
            "4",
            "--seed",
            "20260418",
            "--cosine-tau",
            "0.1",
            "--match-tau",
            "0.1",
            "--min-subjects",
            "10",
            "--max-iter",
            "40",
            "--tol",
            "1e-3",
            "--output",
            out.to_str().unwrap(),
        ]);
        assert!(
            status.status.success(),
            "{}",
            String::from_utf8_lossy(&status.stderr)
        );
    }
    let a = std::fs::read_to_string(&out1).unwrap();
    let b = std::fs::read_to_string(&out2).unwrap();
    assert_eq!(
        a, b,
        "align bootstrap must be deterministic under fixed seed"
    );
}

/// The match-rate diagnostic must fire exactly when the summary
/// contains a conditional row that would mislead.
///
/// The warning is the entire user-facing mechanism for the
/// conditional-column defect: the percentile CI and
/// `bootstrap_mean_n_cohorts` are computed over matched replicates only,
/// so at a low match rate they read far stronger than the archetype is.
/// Nothing else at the terminal says so. If it silently stopped firing,
/// the TSV would look identical and no other test would notice.
///
/// The first version of this test used the planted fixture, where every
/// resample recovers the shared archetype and the match rate sits near
/// 1.0. The iff held, the test passed, and it still passed with the
/// warning compiled out — because the branch it was written to check
/// never ran. Hence `assert!(n_low > 0, ...)`: the precondition that
/// makes the assertion mean anything is itself asserted.
#[test]
fn the_match_rate_warning_fires_exactly_when_a_conditional_row_is_low() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("noise_a");
    let b = tmp.path().join("noise_b");
    write_noise_cohort(&a, "A", 101);
    write_noise_cohort(&b, "B", 202);
    let out = tmp.path().join("bootstrap").join("summary.tsv");

    let res = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        &format!("{},{}", a.display(), b.display()),
        "--labels",
        "A,B",
        "--k",
        "3",
        "--n-boot",
        "20",
        "--seed",
        "20260418",
        "--match-tau",
        "0.7",
        "--min-subjects",
        "10",
        "--max-iter",
        "80",
        "--tol",
        "1e-3",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        res.status.success(),
        "align bootstrap failed:\nstderr:\n{}",
        String::from_utf8_lossy(&res.stderr)
    );
    let err = String::from_utf8_lossy(&res.stderr).to_string();

    let (_, rows) = parse_tsv(&out);
    let rates: Vec<f64> = rows
        .iter()
        .map(|r| {
            r["bootstrap_match_rate"]
                .parse::<f64>()
                .expect("match rate parses")
        })
        .collect();
    assert!(!rates.is_empty(), "fixture produced no archetypes");

    let n_low = rates.iter().filter(|&&r| r < 0.50).count();
    assert!(
        n_low > 0,
        "this fixture must produce at least one archetype below a 0.50 \
         match rate, or the assertions below never exercise the warning. \
         Rates: {rates:?}"
    );

    // The always-on summary line, with the median it reports checked
    // against the file rather than merely present.
    let mut sorted = rates.clone();
    sorted.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let median = sorted[sorted.len() / 2];
    assert!(
        err.contains(&format!("median {median:.3}")),
        "summary line must report the median match rate {median:.3}; got:\n{err}"
    );

    assert!(
        err.contains("matched in fewer than"),
        "{n_low} of {} archetypes are below 0.50, so the warning must \
         fire. stderr:\n{err}",
        rates.len()
    );

    // It must name the columns it is about, or it does not help.
    for col in [
        "bootstrap_mean_n_cohorts",
        "ci_lower/ci_upper_n_cohorts",
        "bootstrap_prob_universal",
    ] {
        assert!(
            err.contains(col),
            "the warning must name {col:?} so the reader knows which \
             columns are affected; got:\n{err}"
        );
    }
}

/// The complement: on a fixture where every archetype recurs, the
/// warning must stay silent. A diagnostic that always fires is noise
/// and gets ignored, which is the same outcome as not having it.
#[test]
fn the_match_rate_warning_stays_silent_when_every_archetype_recurs() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_cohort(&a, "A", true, false, 11);
    write_cohort(&b, "B", false, true, 22);
    let out = tmp.path().join("bootstrap").join("summary.tsv");

    let res = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        &format!("{},{}", a.display(), b.display()),
        "--labels",
        "A,B",
        "--k",
        "2",
        "--n-boot",
        "20",
        "--seed",
        "20260418",
        "--cosine-tau",
        "0.1",
        "--match-tau",
        "0.1",
        "--min-subjects",
        "10",
        "--max-iter",
        "80",
        "--tol",
        "1e-3",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(res.status.success());
    let err = String::from_utf8_lossy(&res.stderr).to_string();

    let (_, rows) = parse_tsv(&out);
    let n_low = rows
        .iter()
        .filter(|r| r["bootstrap_match_rate"].parse::<f64>().unwrap() < 0.50)
        .count();
    assert_eq!(
        n_low, 0,
        "this fixture is meant to recover every archetype; if it no \
         longer does, the silence asserted below is not meaningful"
    );
    assert!(
        !err.contains("matched in fewer than"),
        "no archetype is below 0.50, so the warning must not fire:\n{err}"
    );
}

/// The sidecar's `column_populations` block must classify the columns
/// actually written to the TSV.
///
/// A unit test checks the block against the header constant. This checks
/// it against the file on disk, which is what a downstream figure script
/// reads. The block exists because docstrings and stderr do not reach
/// that reader, so a block that drifted from the real header would be an
/// authoritative-looking record that had stopped describing reality.
#[test]
fn the_sidecar_classifies_every_column_actually_written() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("cohort_a");
    let b = tmp.path().join("cohort_b");
    write_cohort(&a, "A", true, false, 11);
    write_cohort(&b, "B", false, true, 22);
    let out = tmp.path().join("bootstrap").join("summary.tsv");

    let res = run_atman(&[
        "align",
        "bootstrap",
        "--cohorts",
        &format!("{},{}", a.display(), b.display()),
        "--labels",
        "A,B",
        "--k",
        "2",
        "--n-boot",
        "10",
        "--seed",
        "20260418",
        "--cosine-tau",
        "0.1",
        "--match-tau",
        "0.1",
        "--min-subjects",
        "10",
        "--max-iter",
        "80",
        "--tol",
        "1e-3",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(res.status.success());

    let sidecar = out.with_extension("tsv.run.json");
    let text = std::fs::read_to_string(&sidecar)
        .unwrap_or_else(|e| panic!("reading sidecar {sidecar:?}: {e}"));
    let doc: serde_json::Value = serde_json::from_str(&text).unwrap();
    let cp = &doc["column_populations"];
    assert!(
        cp.is_object(),
        "the sidecar must carry column_populations; got:\n{text}"
    );

    let mut classified: Vec<String> = Vec::new();
    for group in ["conditional_on_recovery", "over_all_iterations"] {
        for c in cp[group]["columns"].as_array().expect("columns array") {
            classified.push(c.as_str().unwrap().to_string());
        }
    }

    // Identity and provenance columns are not over a resampling
    // population. Everything else must be classified.
    const NOT_A_STATISTIC: &[&str] = &[
        "archetype_id",
        "observed_n_cohorts",
        "observed_cohorts",
        "bca_fallback_to_percentile",
        "bootstrap_match_rate",
    ];
    let (header, _) = parse_tsv(&out);
    for col in &header {
        if NOT_A_STATISTIC.contains(&col.as_str()) {
            continue;
        }
        assert!(
            classified.contains(col),
            "column {col:?} is written to the TSV but not classified in \
             the sidecar's column_populations"
        );
    }

    // The conditional group must actually name the conditional columns.
    let conditional: Vec<&str> = cp["conditional_on_recovery"]["columns"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    for col in [
        "ci_lower_n_cohorts",
        "ci_upper_n_cohorts",
        "bootstrap_mean_n_cohorts",
    ] {
        assert!(
            conditional.contains(&col),
            "{col:?} is conditional on recovery and must be listed as such"
        );
    }
}
