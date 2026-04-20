//! Priority 10 integration test: `atman decompose unmix`.
//!
//! Builds a synthetic canonical Atman directory with 60 subjects ×
//! 30 proteins drawn from 3 planted endmembers (disjoint feature
//! blocks) plus Dirichlet-mixed abundances with small Gaussian noise,
//! then runs `atman decompose unmix --k 3 --method vca
//! --abundance fcls` and verifies:
//!
//! 1. Output schema: endmembers.tsv, abundances.tsv,
//!    unmix_diagnostics.tsv with the expected columns.
//! 2. Compositional constraints on FCLS abundances: each row sums
//!    to 1 ± 1e-5 and every value is ≥ -1e-9.
//! 3. The recovered endmembers match the planted ones at cosine ≥
//!    0.85 under best-permutation matching.
//! 4. Refuses cleanly when `k > n_samples / 2`.
//! 5. Refuses cleanly when `--transform clr` is combined with
//!    `--abundance fcls` without `--allow-unconstrained-simplex`.

use std::collections::HashMap;
use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn parse_tsv(
    path: &Path,
) -> (Vec<String>, Vec<HashMap<String, String>>) {
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

/// Deterministic LCG + Box-Muller Gaussian stream.
struct Lcg(std::num::Wrapping<u64>);
impl Lcg {
    fn new(seed: u64) -> Self {
        Self(std::num::Wrapping(seed))
    }
    fn next_u(&mut self) -> f64 {
        self.0 = self.0 * std::num::Wrapping(6364136223846793005_u64)
            + std::num::Wrapping(1442695040888963407_u64);
        ((self.0 .0 >> 33) as f64) / (u32::MAX as f64)
    }
    fn next_z(&mut self) -> f64 {
        let a = self.next_u().max(1e-12);
        let b = self.next_u();
        (-2.0_f64 * a.ln()).sqrt() * (2.0 * std::f64::consts::PI * b).cos()
    }
}

/// Write a 60-subject × 30-protein canonical directory with 3
/// planted endmembers on disjoint feature blocks (G001..G010 for
/// E1, G011..G020 for E2, G021..G030 for E3). Subjects 0..2 are
/// pure endmember anchors; subjects 3..59 are Dirichlet-mixed.
fn write_planted_cohort(dir: &Path, seed: u64, planted: &mut [Vec<f64>]) -> Vec<Vec<f64>> {
    std::fs::create_dir_all(dir).unwrap();
    let n = 60usize;
    let p = 30usize;
    let k = 3usize;
    let block = p / k;
    let mut rng = Lcg::new(seed);
    // Planted endmembers: zero background, ~1.0 on own block.
    for i in 0..k {
        for f in 0..block {
            planted[i][i * block + f] = 1.0 + 0.05 * rng.next_z();
        }
    }
    // Abundances: pure anchors at sample 0..2, Dirichlet for rest.
    let mut abundances = vec![vec![0.0_f64; k]; n];
    for i in 0..k {
        abundances[i][i] = 1.0;
    }
    for si in k..n {
        let mut row = vec![0.0_f64; k];
        for v in row.iter_mut() {
            *v = rng.next_u().max(1e-12);
            *v = -v.ln(); // Exp(1)
        }
        let s: f64 = row.iter().sum();
        for v in row.iter_mut() {
            *v /= s;
        }
        abundances[si] = row;
    }

    let mut samples = String::from(
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n",
    );
    for i in 1..=n {
        samples.push_str(&format!(
            "S{i:03}\tS{i:03}\tN/A\t0\tplasma\t{i}\n"
        ));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    let mut proteins =
        String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for j in 1..=p {
        proteins.push_str(&format!(
            "olink_explore_ngs\tA{j:03}\tQ{j:05}\tG{j:03}\tP1\t\n"
        ));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    let mut qc = String::from(
        "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\t\
         abundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\t\
         detection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n",
    );
    let mut order: u64 = 0;
    for si in 0..n {
        for j in 0..p {
            order += 1;
            let mut val = 0.0;
            for ei in 0..k {
                val += abundances[si][ei] * planted[ei][j];
            }
            val += 0.01 * rng.next_z();
            qc.push_str(&format!(
                "olink_explore_ngs\tS{si:03}\tA{j:03}\tG{j:03}\tP1\t{val:.6}\t\
                 {val:.6}\t{val:.6}\tlog2_npx\tPASS\tPASS\t\t0\t0\t\t\t{order}\n",
                si = si + 1,
                j = j + 1,
            ));
        }
    }
    std::fs::write(dir.join("qc_measurements.tsv"), &qc).unwrap();
    std::fs::write(dir.join("measurements.tsv"), &qc).unwrap();
    abundances
}

fn best_permuted_cosine(
    planted: &[Vec<f64>],
    recovered: &[Vec<f64>],
) -> Vec<f64> {
    planted
        .iter()
        .map(|p| {
            recovered
                .iter()
                .map(|r| {
                    let dot: f64 = p.iter().zip(r.iter()).map(|(a, b)| a * b).sum();
                    let np = p.iter().map(|v| v * v).sum::<f64>().sqrt();
                    let nr = r.iter().map(|v| v * v).sum::<f64>().sqrt();
                    if np > 0.0 && nr > 0.0 {
                        (dot / (np * nr)).abs()
                    } else {
                        0.0
                    }
                })
                .fold(0.0_f64, f64::max)
        })
        .collect()
}

#[test]
fn unmix_recovers_planted_endmembers_and_satisfies_simplex_constraints() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("cohort");
    let output = tmp.path().join("unmix_out");
    let mut planted = vec![vec![0.0_f64; 30]; 3];
    write_planted_cohort(&input, 20260420, &mut planted);

    let status = run_atman(&[
        "decompose", "unmix",
        "--input-dir", input.to_str().unwrap(),
        "--k", "3",
        "--method", "vca",
        "--abundance", "fcls",
        "--transform", "none",
        "--seed", "42",
        "--fcls-max-iter", "1000",
        "--fcls-tol", "1e-9",
        "--output-dir", output.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "decompose unmix failed:\nstderr:\n{}",
        String::from_utf8_lossy(&status.stderr)
    );

    // Endmembers TSV schema.
    let (header_em, rows_em) = parse_tsv(&output.join("endmembers.tsv"));
    for col in [
        "endmember_id", "source_sample_id", "protein", "loading", "rank_in_endmember",
    ] {
        assert!(header_em.iter().any(|h| h == col), "missing {col}");
    }
    // 3 endmembers × 30 proteins = 90 rows.
    assert_eq!(rows_em.len(), 90);

    // Reconstruct recovered endmember matrix [k × p] from TSV.
    let mut recovered = vec![vec![0.0_f64; 30]; 3];
    for r in &rows_em {
        let em: &str = r["endmember_id"].as_str();
        let ei: usize = em
            .trim_start_matches('E')
            .parse::<usize>()
            .unwrap()
            - 1;
        let prot: &str = r["protein"].as_str();
        let pi: usize = prot
            .trim_start_matches('G')
            .parse::<usize>()
            .unwrap()
            - 1;
        recovered[ei][pi] = r["loading"].parse().unwrap();
    }
    let cos = best_permuted_cosine(&planted, &recovered);
    for (i, c) in cos.iter().enumerate() {
        assert!(
            *c >= 0.85,
            "planted endmember {i} not recovered at cosine ≥ 0.85; got {c}"
        );
    }

    // Abundance TSV: sample_id + E{i}/ci_lower/ci_upper per endmember.
    // 60 rows, each sums to ~1 across point-estimate columns, non-
    // negative entries.
    let (header_ab, rows_ab) = parse_tsv(&output.join("abundances.tsv"));
    assert_eq!(header_ab.len(), 1 + 3 * 3);
    assert_eq!(rows_ab.len(), 60);
    for r in &rows_ab {
        let mut s = 0.0_f64;
        for col in &["E001", "E002", "E003"] {
            let v: f64 = r[*col].parse().unwrap();
            assert!(
                v >= -1e-9,
                "negative abundance for {}: {v}",
                r["sample_id"]
            );
            s += v;
        }
        assert!(
            (s - 1.0).abs() < 1e-5,
            "abundance row-sum for {} = {s} (expected 1.0 ± 1e-5)",
            r["sample_id"]
        );
    }

    // Diagnostics TSV carries residual + row-sum sanity.
    let (header_d, rows_d) = parse_tsv(&output.join("unmix_diagnostics.tsv"));
    for col in ["sample_id", "reconstruction_residual_norm", "abundance_sum"] {
        assert!(header_d.iter().any(|h| h == col), "missing {col}");
    }
    assert_eq!(rows_d.len(), 60);
    for r in &rows_d {
        let residual: f64 = r["reconstruction_residual_norm"].parse().unwrap();
        // Residual should be small (noise-scale); bound ≤ 1 is
        // loose but catches catastrophic misfits.
        assert!(
            residual >= 0.0 && residual < 1.0,
            "implausibly large residual for {}: {residual}",
            r["sample_id"]
        );
    }

    // Sidecar exists.
    let sidecar = output.join("endmembers.tsv.run.json");
    assert!(sidecar.exists());
}

#[test]
fn unmix_refuses_k_greater_than_half_samples() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("cohort");
    let output = tmp.path().join("out");
    let mut planted = vec![vec![0.0_f64; 30]; 3];
    write_planted_cohort(&input, 1, &mut planted);

    // k = 40 > n/2 = 30.
    let out = run_atman(&[
        "decompose", "unmix",
        "--input-dir", input.to_str().unwrap(),
        "--k", "40",
        "--output-dir", output.to_str().unwrap(),
    ]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("under-determined") || stderr.contains("n/2"),
        "expected under-determined error, got: {stderr}"
    );
}

#[test]
fn unmix_refuses_clr_plus_fcls_without_escape_hatch() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("cohort");
    let output = tmp.path().join("out");
    let mut planted = vec![vec![0.0_f64; 30]; 3];
    write_planted_cohort(&input, 2, &mut planted);

    let out = run_atman(&[
        "decompose", "unmix",
        "--input-dir", input.to_str().unwrap(),
        "--k", "3",
        "--transform", "clr",
        "--abundance", "fcls",
        "--output-dir", output.to_str().unwrap(),
    ]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("allow-unconstrained-simplex") || stderr.contains("refused"),
        "expected CLR+FCLS refusal, got: {stderr}"
    );
}

#[test]
fn unmix_k_auto_selects_planted_k_on_planted_3_fixture() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("cohort");
    let output = tmp.path().join("k_auto_out");
    let mut planted = vec![vec![0.0_f64; 30]; 3];
    write_planted_cohort(&input, 20260420, &mut planted);

    let status = run_atman(&[
        "decompose", "unmix",
        "--input-dir", input.to_str().unwrap(),
        "--k", "auto",
        "--k-min", "2",
        "--k-max", "6",
        "--k-elbow-threshold", "0.10",
        "--method", "vca",
        "--abundance", "fcls",
        "--transform", "none",
        "--seed", "42",
        "--output-dir", output.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "decompose unmix --k auto failed:\nstderr:\n{}",
        String::from_utf8_lossy(&status.stderr)
    );
    let k_path = output.join("k_selection.tsv");
    assert!(k_path.exists(), "k_selection.tsv must be emitted under --k auto");
    let (header, rows) = parse_tsv(&k_path);
    for col in ["k", "mean_residual_norm", "marginal_improvement", "chosen"] {
        assert!(header.iter().any(|h| h == col), "missing {col}");
    }
    let chosen_rows: Vec<&HashMap<String, String>> = rows
        .iter()
        .filter(|r| r["chosen"] == "1")
        .collect();
    assert_eq!(chosen_rows.len(), 1, "expected exactly one chosen k");
    let chosen_k: usize = chosen_rows[0]["k"].parse().unwrap();
    assert!(
        chosen_k == 3 || chosen_k == 4,
        "--k auto should pick ≈ 3 on planted-3 fixture; got {chosen_k}"
    );
    // Final endmembers/abundances reflect the chosen k.
    let (_, em_rows) = parse_tsv(&output.join("endmembers.tsv"));
    let distinct: std::collections::HashSet<String> =
        em_rows.iter().map(|r| r["endmember_id"].clone()).collect();
    assert_eq!(distinct.len(), chosen_k);
}

#[test]
fn unmix_n_boot_emits_ci_columns_with_finite_ordered_bounds() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("cohort");
    let output = tmp.path().join("ci_out");
    let mut planted = vec![vec![0.0_f64; 30]; 3];
    write_planted_cohort(&input, 20260420, &mut planted);

    let status = run_atman(&[
        "decompose", "unmix",
        "--input-dir", input.to_str().unwrap(),
        "--k", "3",
        "--n-boot", "20",
        "--seed", "42",
        "--output-dir", output.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "decompose unmix --n-boot 20 failed:\nstderr:\n{}",
        String::from_utf8_lossy(&status.stderr)
    );

    let (header_em, rows_em) = parse_tsv(&output.join("endmembers.tsv"));
    for col in ["loading_ci_lower", "loading_ci_upper"] {
        assert!(header_em.iter().any(|h| h == col), "missing {col} on endmembers.tsv");
    }
    for r in &rows_em {
        let lo: f64 = r["loading_ci_lower"].parse().unwrap();
        let hi: f64 = r["loading_ci_upper"].parse().unwrap();
        assert!(lo.is_finite() && hi.is_finite());
        assert!(lo <= hi + 1e-9);
    }
    let (header_ab, rows_ab) = parse_tsv(&output.join("abundances.tsv"));
    for col in ["E001_ci_lower", "E001_ci_upper"] {
        assert!(header_ab.iter().any(|h| h == col), "missing {col} on abundances.tsv");
    }
    for r in &rows_ab {
        let lo: f64 = r["E001_ci_lower"].parse().unwrap();
        let hi: f64 = r["E001_ci_upper"].parse().unwrap();
        assert!(lo.is_finite() && hi.is_finite());
        assert!(lo <= hi + 1e-9);
    }
}

#[test]
fn unmix_is_deterministic_under_fixed_seed() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("cohort");
    let out1 = tmp.path().join("r1");
    let out2 = tmp.path().join("r2");
    let mut planted = vec![vec![0.0_f64; 30]; 3];
    write_planted_cohort(&input, 3, &mut planted);

    for out in [&out1, &out2] {
        let status = run_atman(&[
            "decompose", "unmix",
            "--input-dir", input.to_str().unwrap(),
            "--k", "3",
            "--seed", "42",
            "--output-dir", out.to_str().unwrap(),
        ]);
        assert!(status.status.success());
    }
    let a = std::fs::read_to_string(out1.join("endmembers.tsv")).unwrap();
    let b = std::fs::read_to_string(out2.join("endmembers.tsv")).unwrap();
    assert_eq!(a, b, "endmembers must be deterministic under fixed seed");
    let c = std::fs::read_to_string(out1.join("abundances.tsv")).unwrap();
    let d = std::fs::read_to_string(out2.join("abundances.tsv")).unwrap();
    assert_eq!(c, d, "abundances must be deterministic under fixed seed");
}
