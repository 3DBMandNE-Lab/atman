//! Equicorrelated multivariate-t CDF — `pdunnett` / `qdunnett`.
//!
//! Dunnett's test compares every non-control group mean to a single
//! reference-level ("control") mean and adjusts the family-wise
//! error rate via the joint distribution of the correlated
//! t-statistics. For `m` comparisons with common correlation `ρ`
//! and residual degrees of freedom `ν`, the null joint follows a
//! multivariate t with correlation matrix
//!
//! ```text
//! R = ρ · 1 1' + (1 − ρ) · I   (m × m)
//! ```
//!
//! Using the Dunnett–Curnow representation `T_i = (√ρ · W₀ + √(1−ρ) · W_i) / S`
//! with `W₀, W_i ~ iid N(0, 1)` and `S = √(χ²_ν / ν)` independent,
//! the two-sided CDF `P(max_i |T_i| ≤ q)` reduces to a 2D integral:
//!
//! ```text
//! P(max |T_i| ≤ q)
//!   = ∫₀^∞ f_S(s; ν) ∫_{-∞}^{∞} φ(w) · H_m(q, s, w; ρ) dw ds
//!
//! H_m(q, s, w; ρ)
//!   = [Φ((q·s − √ρ·w) / √(1−ρ)) − Φ((−q·s − √ρ·w) / √(1−ρ))]^m
//! ```
//!
//! We evaluate that integral via nested Gauss–Legendre quadrature,
//! identical in structure to `studentized_range::ptukey`. For the
//! balanced-design Dunnett case `ρ = 0.5` is the canonical choice
//! and matches `emmeans::contrast(method="dunnett")`. Unequal
//! sample sizes induce slight correlation heterogeneity that the
//! equicorrelation formula approximates to high accuracy when the
//! n_i are within ~2× of each other; a true Dunnett–Hsu
//! implementation would need `mvtnorm::pmvt` (Genz–Bretz) over an
//! arbitrary correlation matrix and is a documented follow-on.
//!
//! Parity target: `|pdunnett_atman − (1 − p_R_emmeans)| ≤ 3e-3`
//! at canonical `(m, ν)` pairs.

use crate::rng::Xoshiro256pp;
use statrs::distribution::{ContinuousCDF, Normal};
use statrs::function::gamma::ln_gamma;

const N_NODES_INNER: usize = 64;
const N_NODES_OUTER: usize = 48;

fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = vec![0.0; n];
    let mut weights = vec![0.0; n];
    for i in 0..n {
        let theta = std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5);
        let mut x = theta.cos();
        for _ in 0..50 {
            let (p, pp) = legendre_pn_and_derivative(n, x);
            let dx = p / pp;
            x -= dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        nodes[i] = x;
        let (_, pp) = legendre_pn_and_derivative(n, x);
        weights[i] = 2.0 / ((1.0 - x * x) * pp * pp);
    }
    (nodes, weights)
}

fn legendre_pn_and_derivative(n: usize, x: f64) -> (f64, f64) {
    if n == 0 {
        return (1.0, 0.0);
    }
    if n == 1 {
        return (x, 1.0);
    }
    let mut p_prev = 1.0;
    let mut p_curr = x;
    for k in 1..n {
        let kf = k as f64;
        let p_next = ((2.0 * kf + 1.0) * x * p_curr - kf * p_prev) / (kf + 1.0);
        p_prev = p_curr;
        p_curr = p_next;
    }
    let pd = (n as f64) * (x * p_curr - p_prev) / (x * x - 1.0);
    (p_curr, pd)
}

fn sqrt_chi_density(s: f64, df: f64) -> f64 {
    if s <= 0.0 || !s.is_finite() || df <= 0.0 {
        return 0.0;
    }
    let half = 0.5 * df;
    let log_density =
        2.0_f64.ln() + half * half.ln() - ln_gamma(half) + (df - 1.0) * s.ln() - half * s * s;
    log_density.exp()
}

/// Inner integrand over `w` given `q·s`. Maps `w ∈ R` onto the unit
/// interval via `w = Φ⁻¹(u)` so the integration domain is finite.
fn inner_h_m_two_sided(qs: f64, m: usize, rho: f64, nodes: &[f64], weights: &[f64]) -> f64 {
    if m == 0 {
        return 1.0;
    }
    let normal = match Normal::new(0.0, 1.0) {
        Ok(n) => n,
        Err(_) => return f64::NAN,
    };
    let sqrt_rho = rho.sqrt();
    let sqrt_one_minus_rho = (1.0 - rho).sqrt();
    let mf = m as f64;
    let mut acc = 0.0;
    for (&t, &wt) in nodes.iter().zip(weights.iter()) {
        let u = 0.5 * (t + 1.0);
        if !(0.0 < u && u < 1.0) {
            continue;
        }
        let w = normal.inverse_cdf(u);
        let a = (qs - sqrt_rho * w) / sqrt_one_minus_rho;
        let b = (-qs - sqrt_rho * w) / sqrt_one_minus_rho;
        let inner = normal.cdf(a) - normal.cdf(b);
        let inner = inner.clamp(0.0, 1.0);
        acc += inner.powf(mf) * wt * 0.5;
    }
    acc
}

/// `P(max_{i=1..m} |T_i| ≤ q)` for an equicorrelated multivariate-t
/// with correlation `ρ ∈ (0, 1)` and `ν` residual degrees of
/// freedom. Two-sided. Returns values clamped to `[0, 1]`.
///
/// `ν → ∞` short-circuits the outer integration (S = 1).
pub fn pdunnett(q: f64, m: usize, df: f64, rho: f64) -> f64 {
    if m == 0 {
        return 1.0;
    }
    if !(0.0..1.0).contains(&rho) {
        return f64::NAN;
    }
    if !q.is_finite() || q < 0.0 {
        return 0.0;
    }
    if q == 0.0 {
        return 0.0;
    }
    let (inner_nodes, inner_weights) = gauss_legendre(N_NODES_INNER);
    if !df.is_finite() || df >= 1e8 {
        return inner_h_m_two_sided(q, m, rho, &inner_nodes, &inner_weights).clamp(0.0, 1.0);
    }
    if df <= 0.0 {
        return f64::NAN;
    }
    let s_max = (1.0 + 6.0 / df.sqrt()).max(4.0);
    let half_len = 0.5 * s_max;
    let (outer_nodes, outer_weights) = gauss_legendre(N_NODES_OUTER);
    let mut acc = 0.0;
    for (&t, &wt) in outer_nodes.iter().zip(outer_weights.iter()) {
        let s = half_len * (t + 1.0);
        if !(s.is_finite() && s > 0.0) {
            continue;
        }
        let density = sqrt_chi_density(s, df);
        let inner = inner_h_m_two_sided(q * s, m, rho, &inner_nodes, &inner_weights);
        acc += density * inner * wt * half_len;
    }
    acc.clamp(0.0, 1.0)
}

/// Cholesky factor `L` of a symmetric positive-definite `R` such
/// that `L Lᵀ = R`. Lower-triangular; `L[i][j] = 0` for `j > i`.
/// Returns `None` when `R` is non-SPD.
fn cholesky_factor(r: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let m = r.len();
    if m == 0 || r.iter().any(|row| row.len() != m) {
        return None;
    }
    let mut l = vec![vec![0.0_f64; m]; m];
    for i in 0..m {
        let diag_sum: f64 = (0..i).map(|k| l[i][k] * l[i][k]).sum();
        let diag = r[i][i] - diag_sum;
        if !diag.is_finite() || diag <= 0.0 {
            return None;
        }
        l[i][i] = diag.sqrt();
        for j in (i + 1)..m {
            let off_sum: f64 = (0..i).map(|k| l[i][k] * l[j][k]).sum();
            l[j][i] = (r[j][i] - off_sum) / l[i][i];
        }
    }
    Some(l)
}

/// Monte-Carlo two-sided `P(max_i |T_i| ≤ q)` under a multivariate-t
/// with an **arbitrary** correlation matrix `R` (m × m) and residual
/// degrees of freedom `ν`. This is the unbalanced Dunnett-Hsu variant:
/// when comparison group sizes differ, `R` is no longer equicorrelated
/// (off-diagonals depend on the full pair of `n_i` values) and the
/// equicorrelation-based [`pdunnett`] is inexact.
///
/// Uses the standard multivariate-t construction `T = L z / S` where
/// `L Lᵀ = R`, `z ~ N(0, I_m)`, and `S = √(χ²_ν / ν)` independent. A
/// deterministic Xoshiro256++ is seeded from `seed` via a SplitMix64
/// init so repeated calls at the same seed are byte-equal.
///
/// `n_mc` draws deliver ~√(p(1-p) / n_mc) Monte-Carlo error;
/// 50 000 is a typical default. `df` accepts fractional values but
/// the chi² sampler uses `⌈df⌉` standard-normal-squared draws — fine
/// for the small-to-moderate df regime Dunnett-Hsu targets.
pub fn pdunnett_hsu(q: f64, df: f64, correlation: &[Vec<f64>], n_mc: usize, seed: u64) -> f64 {
    let m = correlation.len();
    if m == 0 || n_mc == 0 {
        return 1.0;
    }
    if !q.is_finite() || q < 0.0 {
        return 0.0;
    }
    if q == 0.0 {
        return 0.0;
    }
    let l = match cholesky_factor(correlation) {
        Some(l) => l,
        None => return f64::NAN,
    };
    let df_ceil = df.ceil() as usize;
    if df_ceil == 0 {
        return f64::NAN;
    }
    let mut rng = Xoshiro256pp::new(seed);
    let mut count = 0_usize;
    for _ in 0..n_mc {
        let mut z = vec![0.0_f64; m];
        for v in z.iter_mut() {
            *v = rng.next_normal();
        }
        let mut t = vec![0.0_f64; m];
        for i in 0..m {
            for j in 0..=i {
                t[i] += l[i][j] * z[j];
            }
        }
        let mut chi_sq = 0.0_f64;
        for _ in 0..df_ceil {
            let v = rng.next_normal();
            chi_sq += v * v;
        }
        let s = (chi_sq / df).sqrt();
        let max_abs = t.iter().map(|v| (v / s).abs()).fold(0.0_f64, f64::max);
        if max_abs <= q {
            count += 1;
        }
    }
    count as f64 / n_mc as f64
}

/// Build the Dunnett-Hsu correlation matrix from group sizes:
/// `ρ_{ij} = sqrt(n_i n_j / ((n_0 + n_i)(n_0 + n_j)))` for `i ≠ j`.
/// `n_control` is the size of the reference group; `n_treatments`
/// are the sizes of the `m` non-control groups in display order.
pub fn dunnett_hsu_correlation_matrix(n_control: usize, n_treatments: &[usize]) -> Vec<Vec<f64>> {
    let m = n_treatments.len();
    let mut r = vec![vec![0.0_f64; m]; m];
    for i in 0..m {
        r[i][i] = 1.0;
        for j in (i + 1)..m {
            let ni = n_treatments[i] as f64;
            let nj = n_treatments[j] as f64;
            let n0 = n_control as f64;
            let rho = (ni * nj / ((n0 + ni) * (n0 + nj))).sqrt();
            r[i][j] = rho;
            r[j][i] = rho;
        }
    }
    r
}

/// Inverse of [`pdunnett`]: `q` such that `P(max |T_i| ≤ q) = p`.
/// Bisection bracketed on `[0, 20]`.
pub fn qdunnett(p: f64, m: usize, df: f64, rho: f64) -> f64 {
    if !(0.0..=1.0).contains(&p) {
        return f64::NAN;
    }
    if p == 0.0 {
        return 0.0;
    }
    if p == 1.0 {
        return f64::INFINITY;
    }
    let mut lo = 0.0_f64;
    let mut hi = 20.0_f64;
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let p_mid = pdunnett(mid, m, df, rho);
        if p_mid < p {
            lo = mid;
        } else {
            hi = mid;
        }
        if (hi - lo).abs() < 1e-7 {
            break;
        }
    }
    0.5 * (lo + hi)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pdunnett_boundary_values() {
        assert_eq!(pdunnett(0.0, 3, 30.0, 0.5), 0.0);
        // Large q ⇒ CDF essentially 1.
        assert!((pdunnett(20.0, 3, 30.0, 0.5) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn pdunnett_is_monotone_nondecreasing_in_q() {
        let grid = [1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0];
        let mut prev = -1.0;
        for q in grid {
            let got = pdunnett(q, 3, 30.0, 0.5);
            assert!(got >= prev - 1e-9, "non-monotone at q={q}: {got} < {prev}");
            prev = got;
        }
    }

    #[test]
    fn pdunnett_is_monotone_nonincreasing_in_m() {
        let mut prev = 2.0;
        for m in 1..=5 {
            let got = pdunnett(3.0, m, 30.0, 0.5);
            assert!(got <= prev + 1e-9, "non-monotone at m={m}: {got} > {prev}");
            prev = got;
        }
    }

    /// Known Dunnett critical values at α = 0.05 two-sided, balanced
    /// (ρ = 0.5). From Miller (1981) / R `qmvt()`:
    /// - m = 2, df = ∞: q* ≈ 2.212
    /// - m = 3, df = ∞: q* ≈ 2.349
    /// - m = 3, df = 30: q* ≈ 2.470
    /// - m = 5, df = 30: q* ≈ 2.695
    #[test]
    fn pdunnett_critical_values_recover_0_95() {
        let cases = [
            (2.212_f64, 2usize, 1e9, 0.5),
            (2.349, 3, 1e9, 0.5),
            (2.470, 3, 30.0, 0.5),
            (2.695, 5, 30.0, 0.5),
        ];
        for (q, m, df, rho) in cases {
            let got = pdunnett(q, m, df, rho);
            assert!(
                (got - 0.95).abs() < 1e-2,
                "pdunnett({q}, {m}, {df}, {rho}) = {got}; expected ≈ 0.95"
            );
        }
    }

    /// For m = 1 the "max of one" collapses to a scalar two-sided
    /// t-test: `P(|T| ≤ q) = 2·F_t(q) − 1`. Verify at df = ∞ where
    /// T is standard normal.
    #[test]
    fn pdunnett_m_one_reduces_to_two_sided_z() {
        use statrs::distribution::Normal as N;
        let normal = N::new(0.0, 1.0).unwrap();
        for q in [1.0_f64, 1.5, 2.0, 2.5, 3.0] {
            let expected = 2.0 * normal.cdf(q) - 1.0;
            let got = pdunnett(q, 1, 1e9, 0.5);
            assert!(
                (got - expected).abs() < 5e-3,
                "pdunnett m=1 df=∞ q={q}: got {got}, expected {expected}"
            );
        }
    }

    #[test]
    fn dunnett_hsu_correlation_matrix_reduces_to_half_in_balanced_case() {
        // Balanced n_i = n_control = 10 → ρ_{ij} = √(100 / (400)) = 0.5.
        let r = dunnett_hsu_correlation_matrix(10, &[10, 10, 10]);
        for (i, row) in r.iter().enumerate() {
            assert!((row[i] - 1.0).abs() < 1e-12);
            for (j, &rij) in row.iter().enumerate() {
                if i != j {
                    assert!((rij - 0.5).abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn pdunnett_hsu_matches_equicorrelated_pdunnett_on_balanced_design() {
        // Balanced design ⇒ Hsu's ρ ≡ 0.5 ⇒ pdunnett_hsu should
        // track pdunnett at the same q within MC noise (few 1e-2).
        let r = dunnett_hsu_correlation_matrix(20, &[20, 20, 20]);
        let eq = pdunnett(2.5, 3, 30.0, 0.5);
        let mc = pdunnett_hsu(2.5, 30.0, &r, 50_000, 20260420);
        assert!(
            (eq - mc).abs() < 2e-2,
            "equicorrelated pdunnett({}) vs MC Hsu ({mc}): drift {}",
            eq,
            (eq - mc).abs()
        );
    }

    #[test]
    fn pdunnett_hsu_is_deterministic_under_fixed_seed() {
        let r = dunnett_hsu_correlation_matrix(10, &[15, 20, 25]);
        let a = pdunnett_hsu(2.0, 20.0, &r, 5_000, 7);
        let b = pdunnett_hsu(2.0, 20.0, &r, 5_000, 7);
        assert!((a - b).abs() < 1e-15);
    }

    #[test]
    fn pdunnett_hsu_larger_unbalance_moves_away_from_equicorrelated() {
        // With strongly unbalanced n_i, the Hsu correlation matrix
        // differs noticeably from ρ = 0.5. The MC CDF should differ
        // from the equicorrelated result (drift > MC noise).
        let r = dunnett_hsu_correlation_matrix(5, &[50, 5, 5]);
        let eq = pdunnett(2.5, 3, 30.0, 0.5);
        let mc = pdunnett_hsu(2.5, 30.0, &r, 50_000, 777);
        assert!(
            mc.is_finite(),
            "Hsu MC should return finite probability, got {mc}",
        );
        assert!(
            (mc - eq).abs() > 1e-3,
            "strongly unbalanced Hsu should drift from equicorrelated; got {} vs {}",
            mc,
            eq
        );
    }

    #[test]
    fn qdunnett_round_trips_with_pdunnett() {
        for (m, df, p) in [
            (2usize, 30.0_f64, 0.95_f64),
            (3, 30.0, 0.99),
            (5, 10.0, 0.95),
        ] {
            let q = qdunnett(p, m, df, 0.5);
            let p_back = pdunnett(q, m, df, 0.5);
            assert!(
                (p_back - p).abs() < 1e-2,
                "qdunnett({p}, {m}, {df}) = {q}; pdunnett of that = {p_back}"
            );
        }
    }
}
