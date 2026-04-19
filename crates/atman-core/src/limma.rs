//! limma 3.x–grade eBayes primitives for `atman de --test limma`.
//!
//! Pure-Rust port of the minimal subset of limma needed to reproduce
//! `eBayes(lmFit(y, X), trend, robust)` + TREAT + moderated F. Entry
//! point: `limma_fit`. Everything else is internal algorithm plumbing
//! exposed for unit testing.
//!
//! Reference: Smyth 2004 (eBayes), Phipson et al. 2013 (TREAT),
//! Phipson et al. 2016 (robust), Ritchie et al. 2015 (limma v3).

/// Digamma ψ(x) = d/dx ln Γ(x). Asymptotic expansion for x ≥ 10,
/// recurrence `ψ(x) = ψ(x+1) - 1/x` for smaller x. The x ≥ 10 cutoff
/// keeps the truncated series (through 1/x⁸) under ~1e-10 absolute
/// error across tested values.
pub fn digamma(x: f64) -> f64 {
    if !x.is_finite() || x <= 0.0 {
        return f64::NAN;
    }
    if x < 10.0 {
        return digamma(x + 1.0) - 1.0 / x;
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    // ψ(x) = ln x - 1/(2x) - 1/(12x²) + 1/(120x⁴) - 1/(252x⁶) + 1/(240x⁸) - …
    x.ln() - 0.5 * inv
        - inv2 * (1.0 / 12.0 - inv2 * (1.0 / 120.0 - inv2 * (1.0 / 252.0 - inv2 / 240.0)))
}

/// Trigamma ψ'(x). Asymptotic for x ≥ 10, recurrence `ψ'(x) = ψ'(x+1) + 1/x²`.
pub fn trigamma(x: f64) -> f64 {
    if !x.is_finite() || x <= 0.0 {
        return f64::NAN;
    }
    if x < 10.0 {
        return trigamma(x + 1.0) + 1.0 / (x * x);
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    let inv3 = inv2 * inv;
    // ψ'(x) = 1/x + 1/(2x²) + 1/(6x³) - 1/(30x⁵) + 1/(42x⁷) - 1/(30x⁹) + …
    inv + 0.5 * inv2
        + inv3 * (1.0 / 6.0 - inv2 * (1.0 / 30.0 - inv2 * (1.0 / 42.0 - inv2 / 30.0)))
}

/// Tetragamma ψ''(x) = d/dx ψ'(x). Used by `trigamma_inverse` Newton step.
pub fn tetragamma(x: f64) -> f64 {
    if !x.is_finite() || x <= 0.0 {
        return f64::NAN;
    }
    if x < 10.0 {
        return tetragamma(x + 1.0) - 2.0 / (x * x * x);
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    let inv3 = inv2 * inv;
    let inv4 = inv2 * inv2;
    // ψ''(x) = -1/x² - 1/x³ - 1/(2x⁴) + 1/(6x⁶) - 1/(6x⁸) + …
    -inv2 - inv3 - 0.5 * inv4 + inv4 * inv2 * (1.0 / 6.0 - inv2 / 6.0)
}

/// Inverse trigamma: given `y > 0`, find `x > 0` such that `trigamma(x) = y`.
/// Newton iteration on `ψ'(x) - y = 0` using `ψ''(x)` as the derivative.
/// Initial guess from the asymptotic `ψ'(x) ≈ 1/x + 1/(2x²)`.
pub fn trigamma_inverse(y: f64) -> f64 {
    if !y.is_finite() || y <= 0.0 {
        return f64::NAN;
    }
    // Asymptotic inverse: for small y (large x), x ≈ 1/y + 1/2; for large y
    // (small x), x ≈ 1/sqrt(y). Blend: use 0.5 + 1/y.
    let mut x = if y > 1e7 {
        1.0 / y.sqrt()
    } else if y < 1e-6 {
        1.0 / y
    } else {
        0.5 + 1.0 / y
    };
    for _ in 0..50 {
        let tri = trigamma(x);
        let tet = tetragamma(x);
        if tet == 0.0 || !tet.is_finite() {
            break;
        }
        let delta = (tri - y) / tet;
        let new_x = x - delta;
        if new_x <= 0.0 {
            x *= 0.5;
            continue;
        }
        x = new_x;
        if delta.abs() < 1e-12 * x.max(1.0) {
            break;
        }
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reference values computed by Wolfram Alpha / R:digamma / R:trigamma / R:psigamma(.,2):
    //   digamma(1)        = -Euler-Mascheroni = -0.5772156649015329
    //   digamma(2)        = 1 - γ            =  0.4227843350984671
    //   digamma(10)       ≈  2.251752589066721
    //   trigamma(1)       = π²/6             ≈  1.6449340668482264
    //   trigamma(10)      ≈  0.10516633568168574
    //   tetragamma(1)     = -2 ζ(3)          ≈ -2.4041138063191886
    //   tetragamma(10)    ≈ -0.011049834970801975  (= -2ζ(3) + 2·Σ_{k=1..9} 1/k³)
    //   trigamma_inverse(π²/6)    → 1
    //   trigamma_inverse(0.10517) → ~10

    #[test]
    fn digamma_matches_tabulated_values() {
        assert!((digamma(1.0) - (-0.5772156649015329)).abs() < 1e-10);
        assert!((digamma(2.0) - 0.4227843350984671).abs() < 1e-10);
        assert!((digamma(10.0) - 2.251752589066721).abs() < 1e-10);
    }

    #[test]
    fn trigamma_matches_tabulated_values() {
        assert!((trigamma(1.0) - 1.6449340668482264).abs() < 1e-10);
        assert!((trigamma(10.0) - 0.10516633568168574).abs() < 1e-10);
    }

    #[test]
    fn tetragamma_matches_tabulated_values() {
        assert!((tetragamma(1.0) - (-2.4041138063191886)).abs() < 1e-10);
        assert!((tetragamma(10.0) - (-0.011049834970801975)).abs() < 1e-8);
    }

    #[test]
    fn trigamma_inverse_round_trips() {
        for x in [0.5_f64, 1.0, 2.0, 5.0, 10.0, 50.0] {
            let y = trigamma(x);
            let x_back = trigamma_inverse(y);
            assert!(
                (x_back - x).abs() < 1e-6,
                "trigamma_inverse({y}) = {x_back}, want {x}"
            );
        }
    }
}
