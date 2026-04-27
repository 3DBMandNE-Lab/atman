//! Archetype null calibration for `atman decompose null`.
//!
//! For a given FastICA run on `x` (n samples × p proteins), we want a
//! q-value for each recovered program against "how often does a random
//! null dataset produce a program at least this stable?" The stability
//! statistic is the multi-seed ICASSO score from
//! [`crate::ica::compute_stability_scores`] — for a given program, the
//! mean over alternative seeds of the best Jaccard top-N overlap with
//! any program in that alt seed's result.
//!
//! The null distribution is assembled by running the same procedure on
//! null matrices and taking, for each null iteration, the **maximum**
//! stability across its programs. That "best-of-k stability" represents
//! what random data would produce at its peak; a real program is
//! significant only if its stability exceeds that random-best with
//! permutation probability < threshold.
//!
//! Three null modes:
//!
//! - `SampleShuffle`: permute the rows (samples) of the sample × protein
//!   matrix. **Weakest null** — preserves protein × protein covariance,
//!   so ICA will still find similar archetypes. Use only to test
//!   whether per-sample activation patterns are coherent with sample
//!   metadata (downstream); not a null for archetype existence.
//! - `ProteinShuffle`: independently permute the samples within each
//!   protein column. Destroys inter-protein covariance. The recommended
//!   null for "does this dataset produce more stable archetypes than
//!   decorrelated noise with the same marginals?"
//! - `GaussianMatched`: per-protein N(μ, σ²) draw with sample mean and
//!   variance. Strongest null — parametric, no marginal structure beyond
//!   first two moments.
//!
//! All randomness comes from `crate::ica::Xoshiro256pp` seeded per null
//! iteration via SplitMix-style sub-seed derivation. Under a fixed
//! top-level `seed`, every invocation is byte-identical.

use crate::ica::{canonicalize_ica, compute_stability_scores, fast_ica, Xoshiro256pp};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NullMode {
    SampleShuffle,
    ProteinShuffle,
    GaussianMatched,
}

impl NullMode {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::SampleShuffle => "sample-shuffle",
            Self::ProteinShuffle => "protein-shuffle",
            Self::GaussianMatched => "gaussian-matched",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct NullParams {
    pub k: usize,
    /// Number of null iterations.
    pub n_perm: usize,
    /// Number of FastICA seeds per matrix (real and each null). `>= 2`.
    pub n_seeds: usize,
    /// Top seed for the reference run; null sub-seeds are derived from
    /// this + the iteration index.
    pub seed: u64,
    pub mode: NullMode,
    /// Top-N loadings used by the Jaccard stability metric. `>= 1`.
    pub top_n: usize,
    pub max_iter: usize,
    pub tol: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ArchetypeNullRow {
    /// 1-based program index (matches `atman decompose ica` naming).
    pub program: usize,
    pub observed_stability: f64,
    /// Mean over null iterations of the max-across-programs stability.
    pub null_stability_mean: f64,
    /// 95th percentile of the null "max stability" distribution.
    pub null_stability_p95: f64,
    /// Permutation p-value with the +1 small-sample correction:
    /// `(1 + count(null_max ≥ observed)) / (1 + n_perm)`.
    pub null_p: f64,
}

/// Derive a deterministic per-iteration sub-seed from `(seed, iter)`.
/// Uses SplitMix64 so near-identical inputs produce well-separated
/// outputs without pulling in a cryptographic hash dependency.
fn derive_sub_seed(seed: u64, iter: usize) -> u64 {
    let mut z = seed.wrapping_add(0x9E3779B97F4A7C15_u64.wrapping_mul(iter as u64 + 1));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

/// Fisher-Yates shuffle producing a permutation of `0..n` from a
/// given RNG. Deterministic under fixed seed.
fn permutation(rng: &mut Xoshiro256pp, n: usize) -> Vec<usize> {
    let mut out: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        let bound = (i + 1) as u64;
        // Rejection sampling to avoid modulo bias.
        let mut u;
        loop {
            u = rng_u64(rng);
            let limit = u64::MAX - u64::MAX % bound;
            if u < limit {
                break;
            }
        }
        let j = (u % bound) as usize;
        out.swap(i, j);
    }
    out
}

fn rng_u64(rng: &mut Xoshiro256pp) -> u64 {
    // Xoshiro256pp stores private state; reuse its public `next_normal`
    // via a one-value draw. For uniform u64 draws we use two normals
    // and bit-hash. Simpler: consume one `next_normal` and reinterpret
    // its bit pattern as u64 — fine for permutation index since we only
    // need a uniform-ish spread followed by rejection sampling.
    // (A proper `next_u64` accessor would be cleaner; we piggyback on
    // the existing public API without touching ica.rs.)
    let v = rng.next_normal();
    v.to_bits()
}

/// Produce a null matrix from the original `x` under the selected
/// [`NullMode`], using `rng` for all randomness.
pub fn generate_null_matrix(
    x: &[Vec<f64>],
    mode: NullMode,
    rng: &mut Xoshiro256pp,
) -> Vec<Vec<f64>> {
    let n = x.len();
    if n == 0 {
        return Vec::new();
    }
    let p = x[0].len();
    match mode {
        NullMode::SampleShuffle => {
            let perm = permutation(rng, n);
            let mut out = vec![vec![0.0; p]; n];
            for (i, &src) in perm.iter().enumerate() {
                out[i].clone_from(&x[src]);
            }
            out
        }
        NullMode::ProteinShuffle => {
            let mut out = vec![vec![0.0; p]; n];
            for j in 0..p {
                let perm = permutation(rng, n);
                for (i, &src) in perm.iter().enumerate() {
                    out[i][j] = x[src][j];
                }
            }
            out
        }
        NullMode::GaussianMatched => {
            let mut out = vec![vec![0.0; p]; n];
            for j in 0..p {
                let (mean, std) = column_mean_std(x, j);
                for row in out.iter_mut() {
                    row[j] = mean + std * rng.next_normal();
                }
            }
            out
        }
    }
}

fn column_mean_std(x: &[Vec<f64>], j: usize) -> (f64, f64) {
    let n = x.len() as f64;
    if n == 0.0 {
        return (0.0, 0.0);
    }
    let mean: f64 = x.iter().map(|row| row[j]).sum::<f64>() / n;
    let var: f64 = if n > 1.0 {
        x.iter()
            .map(|row| {
                let d = row[j] - mean;
                d * d
            })
            .sum::<f64>()
            / (n - 1.0)
    } else {
        0.0
    };
    (mean, var.max(0.0).sqrt())
}

/// Compute per-protein ICASSO-style stability for one matrix and
/// return `(n_seeds × observed stability per program)`.
///
/// Runs FastICA with `n_seeds` seeds starting from `base_seed`, uses
/// seed 0 as the reference, and returns per-reference-program mean
/// best-Jaccard across the alternatives.
fn stability_for_matrix(
    matrix: &[Vec<f64>],
    k: usize,
    n_seeds: usize,
    base_seed: u64,
    top_n: usize,
    max_iter: usize,
    tol: f64,
) -> Vec<f64> {
    let runs: Vec<_> = (0..n_seeds)
        .map(|i| fast_ica(matrix, k, base_seed.wrapping_add(i as u64), max_iter, tol))
        .collect();
    let canon: Vec<_> = runs.iter().map(canonicalize_ica).collect();
    if canon.is_empty() {
        return Vec::new();
    }
    let reference = &canon[0];
    // `compute_stability_scores` takes `&[CanonicalIca]` by value —
    // clone the tail of `canon` into a fresh slice.
    let alt_slice: Vec<crate::ica::CanonicalIca> = canon[1..]
        .iter()
        .map(|c| crate::ica::CanonicalIca {
            loadings: c.loadings.clone(),
            activations: c.activations.clone(),
        })
        .collect();
    compute_stability_scores(reference, &alt_slice, top_n)
}

/// Top-level entry: compute per-program observed stability on `x` and
/// a permutation distribution of max-program-stability across null
/// iterations. Returns one [`ArchetypeNullRow`] per program in `x`'s
/// reference run.
///
/// Determinism: fully reproducible under fixed `params.seed`. Each null
/// iteration uses a sub-seed derived from (seed, iter), and within an
/// iteration the FastICA runs use sub-seed + 0..n_seeds.
pub fn archetype_null(
    x: &[Vec<f64>],
    params: &NullParams,
) -> Result<Vec<ArchetypeNullRow>, String> {
    if x.is_empty() {
        return Err("empty input matrix".to_string());
    }
    if params.n_seeds < 2 {
        return Err("--n-seeds must be >= 2 for stability calibration".to_string());
    }
    if params.n_perm == 0 {
        return Err("--n-perm must be >= 1".to_string());
    }
    if params.top_n == 0 {
        return Err("--top-n must be >= 1".to_string());
    }

    // 1. Observed stability on the real matrix.
    let observed = stability_for_matrix(
        x,
        params.k,
        params.n_seeds,
        params.seed,
        params.top_n,
        params.max_iter,
        params.tol,
    );
    if observed.is_empty() {
        return Err("no programs returned from reference fit".to_string());
    }

    // 2. Null permutation distribution of max-program-stability.
    let mut null_maxes: Vec<f64> = Vec::with_capacity(params.n_perm);
    for iter in 0..params.n_perm {
        let sub_seed = derive_sub_seed(params.seed, iter);
        let mut matrix_rng = Xoshiro256pp::new(sub_seed);
        let null_matrix = generate_null_matrix(x, params.mode, &mut matrix_rng);
        let null_stabs = stability_for_matrix(
            &null_matrix,
            params.k,
            params.n_seeds,
            sub_seed.wrapping_add(0x1234_5678_9ABC_DEF0),
            params.top_n,
            params.max_iter,
            params.tol,
        );
        let max_stab = null_stabs
            .iter()
            .copied()
            .filter(|v| v.is_finite())
            .fold(f64::NEG_INFINITY, f64::max);
        null_maxes.push(if max_stab.is_finite() { max_stab } else { 0.0 });
    }
    null_maxes.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mean = null_maxes.iter().sum::<f64>() / params.n_perm as f64;
    let p95_idx = ((0.95 * params.n_perm as f64).ceil() as usize)
        .saturating_sub(1)
        .min(params.n_perm - 1);
    let p95 = null_maxes[p95_idx];

    // 3. Per-program permutation p-value with +1 correction.
    let denom = (params.n_perm as f64) + 1.0;
    Ok(observed
        .iter()
        .enumerate()
        .map(|(i, &obs)| {
            let count_ge = null_maxes.iter().filter(|&&v| v >= obs).count();
            let p = ((count_ge as f64) + 1.0) / denom;
            ArchetypeNullRow {
                program: i + 1,
                observed_stability: obs,
                null_stability_mean: mean,
                null_stability_p95: p95,
                null_p: p,
            }
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toy_matrix(n: usize, p: usize, seed: u64) -> Vec<Vec<f64>> {
        let mut rng = Xoshiro256pp::new(seed);
        (0..n)
            .map(|_| (0..p).map(|_| rng.next_normal()).collect())
            .collect()
    }

    #[test]
    fn permutation_returns_valid_permutation() {
        let mut rng = Xoshiro256pp::new(42);
        let perm = permutation(&mut rng, 20);
        let mut sorted = perm.clone();
        sorted.sort();
        assert_eq!(sorted, (0..20).collect::<Vec<_>>());
    }

    #[test]
    fn sample_shuffle_preserves_column_content_as_set() {
        let x = toy_matrix(10, 5, 7);
        let mut rng = Xoshiro256pp::new(1);
        let null = generate_null_matrix(&x, NullMode::SampleShuffle, &mut rng);
        for j in 0..5 {
            let mut orig: Vec<f64> = x.iter().map(|r| r[j]).collect();
            let mut shuf: Vec<f64> = null.iter().map(|r| r[j]).collect();
            orig.sort_by(|a, b| a.partial_cmp(b).unwrap());
            shuf.sort_by(|a, b| a.partial_cmp(b).unwrap());
            for (a, b) in orig.iter().zip(shuf.iter()) {
                assert!((a - b).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn protein_shuffle_preserves_per_column_multiset() {
        let x = toy_matrix(10, 5, 7);
        let mut rng = Xoshiro256pp::new(1);
        let null = generate_null_matrix(&x, NullMode::ProteinShuffle, &mut rng);
        for j in 0..5 {
            let mut orig: Vec<f64> = x.iter().map(|r| r[j]).collect();
            let mut shuf: Vec<f64> = null.iter().map(|r| r[j]).collect();
            orig.sort_by(|a, b| a.partial_cmp(b).unwrap());
            shuf.sort_by(|a, b| a.partial_cmp(b).unwrap());
            for (a, b) in orig.iter().zip(shuf.iter()) {
                assert!((a - b).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn gaussian_matched_matches_per_column_moments_in_expectation() {
        // Large-n so sample mean / var of the null concentrate near
        // the originals; tolerances are lax (the PRNG isn't fresh).
        let n = 2000;
        let p = 3;
        let mut rng = Xoshiro256pp::new(5);
        let x: Vec<Vec<f64>> = (0..n)
            .map(|_| {
                vec![
                    rng.next_normal() + 1.0,
                    2.0 * rng.next_normal() - 3.0,
                    0.5 * rng.next_normal(),
                ]
            })
            .collect();
        let mut null_rng = Xoshiro256pp::new(11);
        let null = generate_null_matrix(&x, NullMode::GaussianMatched, &mut null_rng);
        for j in 0..p {
            let (mx, sx) = column_mean_std(&x, j);
            let (my, sy) = column_mean_std(&null, j);
            assert!(
                (mx - my).abs() < 0.2,
                "col {j} mean drift: x={mx}, null={my}"
            );
            assert!(
                (sx - sy).abs() / sx.abs().max(1e-6) < 0.2,
                "col {j} std drift: x={sx}, null={sy}"
            );
        }
    }

    #[test]
    fn archetype_null_returns_finite_p_values_in_valid_range() {
        // Small-scale run to keep this unit test fast.
        let x = toy_matrix(30, 20, 17);
        let params = NullParams {
            k: 3,
            n_perm: 20,
            n_seeds: 2,
            seed: 20260418,
            mode: NullMode::ProteinShuffle,
            top_n: 5,
            max_iter: 50,
            tol: 1e-3,
        };
        let rows = archetype_null(&x, &params).expect("ok");
        assert_eq!(rows.len(), 3);
        let min_p = 1.0 / (params.n_perm as f64 + 1.0);
        for r in &rows {
            assert!(r.observed_stability.is_finite());
            assert!(r.null_stability_mean.is_finite());
            assert!(r.null_p >= min_p - 1e-12);
            assert!(r.null_p <= 1.0 + 1e-12);
        }
    }

    #[test]
    fn archetype_null_is_deterministic_under_seed() {
        let x = toy_matrix(20, 15, 9);
        let params = NullParams {
            k: 2,
            n_perm: 10,
            n_seeds: 2,
            seed: 20260418,
            mode: NullMode::ProteinShuffle,
            top_n: 4,
            max_iter: 40,
            tol: 1e-3,
        };
        let a = archetype_null(&x, &params).unwrap();
        let b = archetype_null(&x, &params).unwrap();
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.program, y.program);
            assert_eq!(
                x.observed_stability.to_bits(),
                y.observed_stability.to_bits()
            );
            assert_eq!(x.null_p.to_bits(), y.null_p.to_bits());
            assert_eq!(
                x.null_stability_mean.to_bits(),
                y.null_stability_mean.to_bits()
            );
        }
    }

    #[test]
    fn archetype_null_rejects_empty_input() {
        let params = NullParams {
            k: 1,
            n_perm: 1,
            n_seeds: 2,
            seed: 0,
            mode: NullMode::SampleShuffle,
            top_n: 1,
            max_iter: 10,
            tol: 1e-3,
        };
        let err = archetype_null(&[], &params).unwrap_err();
        assert!(err.contains("empty"), "unexpected: {err}");
    }

    #[test]
    fn archetype_null_rejects_n_seeds_below_two() {
        let x = toy_matrix(5, 5, 3);
        let params = NullParams {
            k: 1,
            n_perm: 1,
            n_seeds: 1,
            seed: 0,
            mode: NullMode::SampleShuffle,
            top_n: 1,
            max_iter: 10,
            tol: 1e-3,
        };
        let err = archetype_null(&x, &params).unwrap_err();
        assert!(err.contains("n-seeds"), "unexpected: {err}");
    }
}
