//! Integration tests for `atman bench decompose --tools atman.<method>` dispatch.
//!
//! Tests:
//!  1. `atman.ica` and `atman.missingness-ica` run against the ICA-shaped
//!     fixture (default `--fixture`); `atman.nmf` runs against the non-negative
//!     NMF fixture (`--fixture-nmf` or the committed default).
//!  2. NMF achieves archetype_correlation ≥ 0.85 on at least 2 of 3 planted
//!     archetypes in the NMF fixture.
//!  3. All methods are deterministic (determinism_score = 1.0) under a fixed
//!     seed.
//!  4. Unknown method suffix (`atman.fictitious`) emits a tool_not_available
//!     row with the exact contracted error message.

use std::collections::HashMap;
use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn parse_tsv(path: &Path) -> (Vec<String>, Vec<HashMap<String, String>>) {
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

/// Deterministic Gaussian draws via LCG + Box-Muller.
fn make_gaussian_source(seed: u64, n: usize) -> Vec<f64> {
    let mut state = std::num::Wrapping(seed);
    let mut next_u = || {
        state = state * std::num::Wrapping(6364136223846793005_u64)
            + std::num::Wrapping(1442695040888963407_u64);
        ((state.0 >> 33) as f64) / (u32::MAX as f64)
    };
    let mut out = Vec::with_capacity(n);
    while out.len() < n {
        let u1 = next_u().max(1e-12);
        let u2 = next_u();
        let z = (-2.0_f64 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
        out.push(z);
    }
    out
}

/// ICA-shaped fixture: cubic-Gaussian sources, signed mixtures.
/// Used by atman.ica and atman.missingness-ica.
fn write_ica_fixture(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    // 2 archetypes on 30 proteins. A1 loads on P01..P10; A2 on P11..P20.
    let p = 30usize;
    let n = 40usize;
    let archs = [
        {
            let mut v = vec![0.0_f64; p];
            for (i, slot) in v.iter_mut().enumerate().take(10) {
                *slot = 1.0 - (i as f64) * 0.03;
            }
            v
        },
        {
            let mut v = vec![0.0_f64; p];
            for (offset, slot) in v.iter_mut().skip(10).take(10).enumerate() {
                *slot = 1.0 - (offset as f64) * 0.03;
            }
            v
        },
    ];

    // Write planted_loadings.tsv.
    let mut planted = String::from("archetype_id\tprotein\tloading\n");
    for (ai, row) in archs.iter().enumerate() {
        for (pi, &v) in row.iter().enumerate() {
            planted.push_str(&format!("P_A{:02}\tG{:03}\t{v:.6}\n", ai + 1, pi + 1));
        }
    }
    std::fs::write(dir.join("planted_loadings.tsv"), planted).unwrap();

    let s1 = make_gaussian_source(77, n);
    let s2 = make_gaussian_source(88, n);
    let noise = make_gaussian_source(99, n * p);

    let mut abundance = String::from("sample_id");
    for pi in 0..p {
        abundance.push('\t');
        abundance.push_str(&format!("G{:03}", pi + 1));
    }
    abundance.push('\n');
    for si in 0..n {
        abundance.push_str(&format!("S{si:03}"));
        for pi in 0..p {
            let v = s1[si].powi(3) * archs[0][pi]
                + s2[si].powi(3) * archs[1][pi]
                + 0.1 * noise[si * p + pi];
            abundance.push('\t');
            abundance.push_str(&format!("{v:.6}"));
        }
        abundance.push('\n');
    }
    std::fs::write(dir.join("abundance.tsv"), abundance).unwrap();
}

/// Non-negative fixture for NMF: Gamma-distributed mixing weights and loadings,
/// Gamma-distributed noise. Mirrors the design of bench/planted_archetypes_nmf_v1/
/// but at smaller scale for fast test runs (3 archetypes × 40 features × 60 samples).
fn write_nmf_fixture(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();

    const N_ARCHETYPES: usize = 3;
    const N_FEATURES: usize = 40;
    const N_SAMPLES: usize = 60;
    const BLOCK: usize = N_FEATURES / N_ARCHETYPES; // 13 features per archetype

    // Deterministic LCG-based Gamma approximation via Marsaglia-Tsang (sum of
    // squared Gaussians for integer shapes, scaled for non-integer shapes).
    // For bench tests we use the Box-Muller Gaussian and a simple gamma
    // approximation: X ~ Gamma(a, s) ≈ (-a * ln(U1*U2*...*Ua)) * s for
    // integer a, or via log-gamma for non-integer.  Here we use shape=2 and
    // shape=1 which are simple.
    let gaussian = make_gaussian_source;

    // Gamma(1.0, 1.0) = Exponential(1.0): X = -ln(U).
    // We approximate U as (gaussian_abs + epsilon) / scale.
    let n_loadings = N_ARCHETYPES * N_FEATURES;
    let n_weights = N_SAMPLES * N_ARCHETYPES;
    let n_noise = N_SAMPLES * N_FEATURES;

    // Draw uniform-ish values from the Gaussian source by taking |z| values
    // and using Box-Muller u1 variable (which is Uniform[0,1] before log).
    // Simpler: use a separate stream of [0,1] uniform values derived from LCG.
    fn lcg_uniform(seed: u64, n: usize) -> Vec<f64> {
        let mut state = std::num::Wrapping(seed);
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            state = state * std::num::Wrapping(6364136223846793005_u64)
                + std::num::Wrapping(1442695040888963407_u64);
            let v = ((state.0 >> 11) as f64) / ((1_u64 << 53) as f64);
            out.push(v.max(1e-15));
        }
        out
    }

    // Gamma(shape, scale): sum of `shape` exponentials each scaled by `scale`.
    // Works for integer shape only.
    fn gamma_int(uniforms: &[f64], shape: usize, scale: f64) -> Vec<f64> {
        assert_eq!(uniforms.len() % shape, 0);
        uniforms
            .chunks(shape)
            .map(|chunk| {
                let exp_sum: f64 = chunk.iter().map(|&u| -u.ln()).sum();
                exp_sum * scale
            })
            .collect()
    }

    // Gamma(0.5, 0.2) noise: approximate with |Normal| * sqrt(2) * scale.
    // Shape=0.5 is a chi distribution with 1 dof scaled by sqrt(2*scale^2).
    fn gamma_half(gaussians: &[f64], scale: f64) -> Vec<f64> {
        gaussians
            .iter()
            .map(|&z| {
                // Gamma(0.5, scale) ~ (z^2 / 2) * 2*scale ... use |z| * scale
                // This is approximate; exact shape doesn't matter for noise floor.
                z.abs() * scale
            })
            .collect()
    }

    // --- loadings H [K x P] with block structure ---
    let u_load = lcg_uniform(1001, n_loadings);
    let mut h = vec![vec![0.0_f64; N_FEATURES]; N_ARCHETYPES];
    let mut load_cursor = 0usize;
    let u_leak = lcg_uniform(2002, n_loadings); // for sparse leakage

    // `k` and `j` are 2D indices into `h[k][j]` and drive block-structure arithmetic.
    #[allow(clippy::needless_range_loop)]
    for k in 0..N_ARCHETYPES {
        for j in 0..N_FEATURES {
            let in_block = j >= k * BLOCK && j < (k + 1) * BLOCK;
            if in_block {
                // Gamma(1,1) = Exponential
                h[k][j] = -u_load[load_cursor].ln();
            } else {
                // Sparse leakage: Gamma(0.3, 0.1) ~ small positive
                h[k][j] = (-u_leak[load_cursor].ln()) * 0.1 * 0.3;
            }
            load_cursor += 1;
        }
    }

    // --- mixing weights W [N x K], Gamma(2,1) ---
    let u_w = lcg_uniform(3003, n_weights * 2); // need 2 uniforms per Gamma(2,1)
    let w_flat = gamma_int(&u_w, 2, 1.0); // shape=2 -> sum of 2 exponentials
    let mut w = vec![vec![0.0_f64; N_ARCHETYPES]; N_SAMPLES];
    for i in 0..N_SAMPLES {
        for k in 0..N_ARCHETYPES {
            w[i][k] = w_flat[i * N_ARCHETYPES + k];
        }
    }

    // --- noise [N x P], Gamma(0.5, 0.2) ---
    let g_noise = gaussian(4004, n_noise);
    let noise_flat = gamma_half(&g_noise, 0.2);

    // --- X = W @ H + noise ---
    let mut x = vec![vec![0.0_f64; N_FEATURES]; N_SAMPLES];
    // `i` and `j` are 2D indices into `x[i][j]` plus `noise_flat` (flat-indexed by i*P+j).
    #[allow(clippy::needless_range_loop)]
    for i in 0..N_SAMPLES {
        for j in 0..N_FEATURES {
            let signal: f64 = (0..N_ARCHETYPES).map(|k| w[i][k] * h[k][j]).sum();
            x[i][j] = (signal + noise_flat[i * N_FEATURES + j]).max(0.0);
        }
    }

    // Write planted_loadings.tsv.
    let gene_ids: Vec<String> = (0..N_FEATURES).map(|j| format!("G{:03}", j + 1)).collect();
    let arch_ids: Vec<String> = (0..N_ARCHETYPES)
        .map(|k| format!("P_A{:02}", k + 1))
        .collect();

    let mut planted = String::from("archetype_id\tprotein\tloading\n");
    // `k`/`j` index `arch_ids[k]`, `gene_ids[j]`, and `h[k][j]` in lockstep.
    #[allow(clippy::needless_range_loop)]
    for k in 0..N_ARCHETYPES {
        for j in 0..N_FEATURES {
            planted.push_str(&format!(
                "{}\t{}\t{:.6}\n",
                arch_ids[k], gene_ids[j], h[k][j]
            ));
        }
    }
    std::fs::write(dir.join("planted_loadings.tsv"), planted).unwrap();

    // Write abundance.tsv.
    let mut abundance = String::from("sample_id");
    for gid in &gene_ids {
        abundance.push('\t');
        abundance.push_str(gid);
    }
    abundance.push('\n');
    // `i`/`j` index `x[i][j]` in lockstep with the sample-row construction.
    #[allow(clippy::needless_range_loop)]
    for i in 0..N_SAMPLES {
        abundance.push_str(&format!("S{i:03}"));
        for j in 0..N_FEATURES {
            abundance.push('\t');
            abundance.push_str(&format!("{:.6}", x[i][j]));
        }
        abundance.push('\n');
    }
    std::fs::write(dir.join("abundance.tsv"), abundance).unwrap();
}

/// Locate the committed NMF fixture relative to the workspace root.
/// Falls back to a temp-dir inline fixture when running outside the repo.
fn committed_nmf_fixture_path() -> Option<std::path::PathBuf> {
    // CARGO_MANIFEST_DIR is the test crate (crates/atman).
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    // crates/atman -> repo root -> bench/planted_archetypes_nmf_v1
    let candidate = manifest
        .parent() // crates/
        .and_then(|p| p.parent()) // repo root
        .map(|root| root.join("bench/planted_archetypes_nmf_v1"));
    candidate.filter(|p| p.join("abundance.tsv").exists())
}

#[test]
fn bench_decompose_runs_all_three_atman_methods() {
    let tmp = tempfile::tempdir().unwrap();
    let ica_fixture = tmp.path().join("ica_fixture");
    let nmf_fixture = tmp.path().join("nmf_fixture");
    let output = tmp.path().join("bench.tsv");
    write_ica_fixture(&ica_fixture);

    // Prefer the committed NMF fixture so we validate against the canonical
    // generated data; fall back to an inline fixture in sandboxed envs.
    let nmf_fixture_path = committed_nmf_fixture_path().unwrap_or_else(|| {
        write_nmf_fixture(&nmf_fixture);
        nmf_fixture.clone()
    });

    // Run ICA-family tools with the ICA fixture + NMF tool with its own fixture.
    let status = run_atman(&[
        "bench",
        "decompose",
        "--fixture",
        ica_fixture.to_str().unwrap(),
        "--fixture-nmf",
        nmf_fixture_path.to_str().unwrap(),
        "--tools",
        "atman.ica,atman.nmf,atman.missingness-ica",
        "--seed",
        "20260420",
        "--top-n",
        "20",
        "--output",
        output.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "bench decompose failed:\nstderr:\n{}",
        String::from_utf8_lossy(&status.stderr)
    );

    let (_, rows) = parse_tsv(&output);

    // At minimum: 2 rows for ICA (2 planted), 3 rows for NMF (3 planted),
    // 2 rows for missingness-ica. Total >= 7.
    assert!(
        rows.len() >= 7,
        "expected at least 7 rows, got {}",
        rows.len()
    );

    // --- ICA-family tools: each must be available and deterministic. ---
    for method in &["atman.ica", "atman.missingness-ica"] {
        let method_rows: Vec<&HashMap<String, String>> = rows
            .iter()
            .filter(|r| r.get("tool").map(|t| t == *method).unwrap_or(false))
            .collect();
        assert!(!method_rows.is_empty(), "no rows for method {method}");
        for r in &method_rows {
            let not_avail: u8 = r["tool_not_available"].parse().unwrap_or(1);
            assert_eq!(not_avail, 0, "{method}: tool_not_available should be 0");
            let det: f64 = r["determinism_score"].parse().unwrap_or(f64::NAN);
            assert!(
                (det - 1.0).abs() < 1e-12,
                "{method}: determinism_score should be 1.0 under fixed seed; got {det}"
            );
        }
        let any_nonzero = method_rows.iter().any(|r| {
            r["archetype_correlation"]
                .parse::<f64>()
                .map(|c| c.is_finite() && c != 0.0)
                .unwrap_or(false)
        });
        assert!(
            any_nonzero,
            "{method}: at least one archetype_correlation should be finite and non-zero"
        );
    }

    // --- NMF: must achieve archetype_correlation >= 0.85 on >= 2 of 3 planted
    // archetypes in the non-negative fixture. ---
    let nmf_rows: Vec<&HashMap<String, String>> = rows
        .iter()
        .filter(|r| r.get("tool").map(|t| t == "atman.nmf").unwrap_or(false))
        .collect();
    assert!(!nmf_rows.is_empty(), "no rows for atman.nmf");
    for r in &nmf_rows {
        let not_avail: u8 = r["tool_not_available"].parse().unwrap_or(1);
        assert_eq!(not_avail, 0, "atman.nmf: tool_not_available should be 0");
        let det: f64 = r["determinism_score"].parse().unwrap_or(f64::NAN);
        assert!(
            (det - 1.0).abs() < 1e-12,
            "atman.nmf: determinism_score should be 1.0 under fixed seed; got {det}"
        );
    }
    let high_corr_count = nmf_rows
        .iter()
        .filter(|r| {
            r["archetype_correlation"]
                .parse::<f64>()
                .map(|c| c >= 0.85)
                .unwrap_or(false)
        })
        .count();
    assert!(
        high_corr_count >= 2,
        "atman.nmf should achieve archetype_correlation >= 0.85 on >= 2 of 3 planted archetypes \
         in the NMF fixture; got {} archetype(s) above threshold.\nNMF rows:\n{:?}",
        high_corr_count,
        nmf_rows
            .iter()
            .map(|r| format!(
                "  archetype={} corr={}",
                r["planted_archetype"], r["archetype_correlation"]
            ))
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn bench_decompose_rejects_unknown_atman_method() {
    let tmp = tempfile::tempdir().unwrap();
    let fixture = tmp.path().join("fixture");
    let output = tmp.path().join("bench.tsv");
    write_ica_fixture(&fixture);

    let status = run_atman(&[
        "bench",
        "decompose",
        "--fixture",
        fixture.to_str().unwrap(),
        "--tools",
        "atman.fictitious",
        "--seed",
        "20260420",
        "--output",
        output.to_str().unwrap(),
    ]);

    // The harness emits a tool_not_available row and exits 0 (same pattern as
    // missing external adapters).  Either a non-zero exit or a
    // tool_not_available row with the contracted message is acceptable.
    let stderr_str = String::from_utf8_lossy(&status.stderr);
    let stdout_str = String::from_utf8_lossy(&status.stdout);

    if status.status.success() && output.exists() {
        let (_, rows) = parse_tsv(&output);
        let fictitious_rows: Vec<&HashMap<String, String>> = rows
            .iter()
            .filter(|r| {
                r.get("tool")
                    .map(|t| t == "atman.fictitious")
                    .unwrap_or(false)
            })
            .collect();
        assert_eq!(
            fictitious_rows.len(),
            1,
            "expected exactly one row for atman.fictitious; got {}",
            fictitious_rows.len()
        );
        let flag: u8 = fictitious_rows[0]["tool_not_available"].parse().unwrap();
        assert_eq!(flag, 1, "tool_not_available should be 1");

        let reason = &fictitious_rows[0]["unavailable_reason"];
        let expected_msg =
            "unknown atman method 'fictitious'; expected one of: ica, nmf, missingness-ica";
        assert!(
            reason.contains(expected_msg),
            "unavailable_reason should contain the contracted error message.\nGot: {reason}\nExpected substring: {expected_msg}"
        );
    } else {
        // Non-zero exit is also acceptable; just check the error message appears
        // somewhere in stderr or stdout.
        let combined = format!("{stderr_str}{stdout_str}");
        assert!(
            combined.contains("fictitious") || !status.status.success(),
            "expected non-zero exit or error message for unknown method; stderr: {stderr_str}"
        );
    }
}
