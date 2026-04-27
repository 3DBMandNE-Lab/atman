//! Studentized range distribution — CDF and inverse CDF.
//!
//! The studentized range
//!
//! ```text
//! Q = (max(X_1..X_n) - min(X_1..X_n)) / S
//! ```
//!
//! where `X_1..X_n` are independent N(0, σ²) and `S` is an
//! independent estimator of `σ` with `ν` degrees of freedom
//! (`S = σ · √(χ²_ν / ν)`), governs Tukey's HSD pairwise comparisons
//! and related simultaneous-inference procedures. Its CDF has no
//! closed form; R's `ptukey` / `qtukey` evaluate it numerically via
//! Copenhaver–Holland's AS 190 algorithm.
//!
//! We take the same approach but with a pure-Rust implementation:
//! two nested Gauss–Legendre quadratures, one over the chi-density
//! in the `S` variable (outer) and one over the normal range
//! integrand in the reduced variable (inner). Node/weight
//! computation uses Newton iteration on Legendre polynomial zeros
//! via Bonnet's three-term recurrence, so this module has no
//! external data tables.
//!
//! Parity target: `|ptukey_atman(q, n, ν) − ptukey_R(q, n, ν)| ≤ 1e-3`
//! for `q ∈ [0.5, 8]`, `n ∈ {2..10}`, `ν ∈ {2, 5, 10, 30, 100, ∞}`.
//!
//! References:
//! - Copenhaver, M. D. & Holland, B. S. (1988). *Computing the
//!   exact studentized range statistic*. Journal of Statistical
//!   Computation and Simulation 30(1), 1–15.
//! - R source `src/nmath/ptukey.c`.

use statrs::distribution::{ContinuousCDF, Normal};
use statrs::function::gamma::ln_gamma;

/// Number of Gauss–Legendre nodes used for the inner (normal-range)
/// integration. The inner integrand `[Φ(Φ⁻¹(u) + w) − u]^{n−1}` has
/// a sharp peak near `u ∈ (0, 0.1)` for practical `q ≥ 3`, so we
/// need a generous node count to resolve it at ≤ 1e-3 accuracy.
const N_NODES_INNER: usize = 128;

/// Number of nodes for the outer (chi-density) integration. The
/// outer integrand concentrates near `s = 1` with scale `1/√(2ν)`;
/// 48 nodes give ≥ 5 decimals for `ν ≥ 1`.
const N_NODES_OUTER: usize = 48;

/// Compute Gauss–Legendre nodes and weights on the interval
/// [-1, 1] via Newton iteration on Legendre polynomial zeros.
/// Returns `(nodes, weights)` both of length `n`.
fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = vec![0.0; n];
    let mut weights = vec![0.0; n];
    for i in 0..n {
        // Initial guess from the Chebyshev approximation.
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

/// Legendre polynomial P_n(x) and its derivative P'_n(x) via
/// Bonnet's recursion. For `n = 0` returns `(1, 0)`; for `n = 1`
/// returns `(x, 1)`.
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

/// Inner integrand: `P(range(n i.i.d. standard normals) ≤ w)`.
/// Uses the substitution `u = Φ(x)` to fold the infinite integral
/// onto the unit interval:
///
/// ```text
/// W(w, n) = n · ∫_0^1 [Φ(Φ⁻¹(u) + w) − u]^(n−1) du
/// ```
fn range_cdf_infinite_df(w: f64, n: usize, nodes: &[f64], weights: &[f64]) -> f64 {
    if w <= 0.0 {
        return 0.0;
    }
    if !w.is_finite() {
        return 1.0;
    }
    let normal = match Normal::new(0.0, 1.0) {
        Ok(n) => n,
        Err(_) => return f64::NAN,
    };
    let nf = n as f64;
    let mut acc = 0.0;
    for (&t, &wt) in nodes.iter().zip(weights.iter()) {
        // Map t ∈ (-1, 1) → u ∈ (0, 1).
        let u = 0.5 * (t + 1.0);
        if !(0.0 < u && u < 1.0) {
            continue;
        }
        let x = normal.inverse_cdf(u);
        let inner = normal.cdf(x + w) - u;
        // Clamp tiny negative values coming from cancellation.
        let inner = inner.max(0.0);
        let integrand = inner.powf(nf - 1.0);
        acc += integrand * wt * 0.5;
    }
    nf * acc
}

/// Density of `S = √(χ²_ν / ν)` at `s`. For `ν → ∞` this
/// concentrates at `s = 1`; the caller handles that case separately.
fn sqrt_chi_density(s: f64, df: f64) -> f64 {
    if s <= 0.0 || !s.is_finite() || df <= 0.0 {
        return 0.0;
    }
    // log f_S(s) = log 2 + (df/2) log(df/2) − lnΓ(df/2) + (df−1) log s
    //             − (df · s² / 2)
    let half = 0.5 * df;
    let log_density =
        2.0_f64.ln() + half * half.ln() - ln_gamma(half) + (df - 1.0) * s.ln() - half * s * s;
    log_density.exp()
}

/// `P(Q ≤ q | n, ν)` — studentized range CDF with `n` groups and
/// `ν` residual degrees of freedom.
///
/// - `q < 0` returns 0.
/// - `ν = ∞` (or `ν` above 1e8) short-circuits the outer integral.
pub fn ptukey(q: f64, n: usize, df: f64) -> f64 {
    if n < 2 {
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
        return range_cdf_infinite_df(q, n, &inner_nodes, &inner_weights).clamp(0.0, 1.0);
    }
    if df <= 0.0 {
        return f64::NAN;
    }
    // Outer integration domain: S concentrates near 1 with scale 1/√(2ν).
    // s_max = max(4, 1 + 6/√ν) covers ≥ 1 − 10⁻¹² of the chi mass at
    // any practical ν (ν = 1 ⇒ s_max = 7; ν = 100 ⇒ s_max = 4).
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
        let inner = range_cdf_infinite_df(q * s, n, &inner_nodes, &inner_weights);
        acc += density * inner * wt * half_len;
    }
    acc.clamp(0.0, 1.0)
}

/// Inverse of [`ptukey`]: find `q` such that `P(Q ≤ q | n, ν) = p`.
/// Uses bisection bracketed on `[0, 30]`, which covers any practical
/// critical value up to ~1e-15 in the upper tail.
pub fn qtukey(p: f64, n: usize, df: f64) -> f64 {
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
    let mut hi = 30.0_f64;
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let p_mid = ptukey(mid, n, df);
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
    fn ptukey_boundary_values() {
        assert_eq!(ptukey(0.0, 3, 10.0), 0.0);
        assert!(ptukey(f64::INFINITY, 3, 10.0).is_finite());
        // For a huge q the CDF must be essentially 1.
        assert!((ptukey(50.0, 3, 10.0) - 1.0).abs() < 1e-6);
        // Negative q returns 0.
        assert_eq!(ptukey(-1.0, 3, 10.0), 0.0);
    }

    #[test]
    fn gauss_legendre_integrates_polynomial_exactly() {
        // Gauss-Legendre of order n exactly integrates polynomials
        // of degree up to 2n-1. Check: ∫_{-1}^{1} x^4 dx = 2/5.
        let (nodes, weights) = gauss_legendre(8);
        let integral: f64 = nodes
            .iter()
            .zip(weights.iter())
            .map(|(x, w)| x.powi(4) * w)
            .sum();
        assert!((integral - 2.0 / 5.0).abs() < 1e-12);
    }

    /// Canonical α=0.05 critical values from standard studentized
    /// range tables: `qtukey(0.95, k, df)`. For each, verify that
    /// `ptukey(q_crit, k, df) ≈ 0.95`.
    #[test]
    fn ptukey_critical_values_recover_0_95() {
        let cases = [
            (3.4864_f64, 3usize, 30.0_f64),
            (3.9567, 4, 20.0),
            (4.6543, 5, 10.0),
            (3.3145, 3, 1e9), // df = ∞
        ];
        for (q, n, df) in cases {
            let got = ptukey(q, n, df);
            assert!(
                (got - 0.95).abs() < 5e-3,
                "ptukey({q}, {n}, {df}) = {got}; expected ≈ 0.95"
            );
        }
    }

    /// ptukey is strictly increasing in q. Verify monotonicity on
    /// a small grid — a basic sanity check that the CDF doesn't
    /// oscillate.
    #[test]
    fn ptukey_is_monotone_nondecreasing_in_q() {
        let grid = [1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0];
        let mut prev = -1.0;
        for q in grid {
            let got = ptukey(q, 3, 30.0);
            assert!(got >= prev - 1e-9, "non-monotone at q={q}: {got} < {prev}");
            assert!((0.0..=1.0 + 1e-9).contains(&got));
            prev = got;
        }
    }

    /// At fixed q and df, ptukey is strictly decreasing in n_means
    /// (more groups ⇒ wider expected range ⇒ smaller CDF at any q).
    #[test]
    fn ptukey_is_monotone_nonincreasing_in_nmeans() {
        let mut prev = 2.0;
        for n in 2..=6 {
            let got = ptukey(4.0, n, 30.0);
            assert!(got <= prev + 1e-9, "non-monotone at n={n}: {got} > {prev}");
            prev = got;
        }
    }

    #[test]
    fn qtukey_inverts_ptukey_to_tight_tolerance() {
        for (n, df, p) in [
            (3usize, 30.0_f64, 0.95_f64),
            (4, 20.0, 0.95),
            (5, 10.0, 0.95),
            (3, 5.0, 0.99),
        ] {
            let q = qtukey(p, n, df);
            let p_back = ptukey(q, n, df);
            assert!(
                (p_back - p).abs() < 5e-3,
                "qtukey({p}, {n}, {df}) = {q}; ptukey of that = {p_back}"
            );
        }
    }
}
