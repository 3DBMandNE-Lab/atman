//! Shared statistics primitives used across commands.
//!
//! Every function here is pure, deterministic, and takes `&[f64]` slices.
//! Correlation helpers return `Option<f64>` (None when the computation is
//! undefined — constant vector, length mismatch, or empty input). Callers
//! adapt the `Option` to their local error style: `anyhow::bail` when the
//! operation has user-visible prerequisites, or a sentinel like `0.0` when
//! the result is a matrix cell that must stay numeric.

use std::cmp::Ordering;
use std::collections::BTreeSet;

/// Arithmetic mean; returns 0.0 for an empty slice.
pub fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        0.0
    } else {
        values.iter().sum::<f64>() / values.len() as f64
    }
}

/// Average-tie fractional ranks (1-based). Equivalent to `scipy.stats.rankdata`
/// with `method='average'`.
pub fn ranks(values: &[f64]) -> Vec<f64> {
    let mut indexed: Vec<(usize, f64)> = values.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));
    let mut out = vec![0.0; values.len()];
    let mut i = 0;
    while i < indexed.len() {
        let mut j = i + 1;
        while j < indexed.len() && indexed[j].1 == indexed[i].1 {
            j += 1;
        }
        let rank = (i + 1 + j) as f64 / 2.0;
        for k in i..j {
            out[indexed[k].0] = rank;
        }
        i = j;
    }
    out
}

/// Pearson correlation on paired slices. Returns `None` when lengths differ
/// or either input has zero variance. Result is clamped to `[-1, 1]`.
pub fn pearson(xs: &[f64], ys: &[f64]) -> Option<f64> {
    if xs.len() != ys.len() || xs.is_empty() {
        return None;
    }
    let mx = mean(xs);
    let my = mean(ys);
    let mut sxx = 0.0_f64;
    let mut syy = 0.0_f64;
    let mut sxy = 0.0_f64;
    for (x, y) in xs.iter().zip(ys.iter()) {
        let dx = x - mx;
        let dy = y - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    if sxx == 0.0 || syy == 0.0 {
        return None;
    }
    Some((sxy / (sxx.sqrt() * syy.sqrt())).clamp(-1.0, 1.0))
}

/// Spearman rank correlation; `None` when lengths differ or ranks are constant.
pub fn spearman(xs: &[f64], ys: &[f64]) -> Option<f64> {
    if xs.len() != ys.len() {
        return None;
    }
    pearson(&ranks(xs), &ranks(ys))
}

/// Cosine similarity of two equal-length vectors. Returns `None` when lengths
/// differ or either norm is zero. Caller chooses whether to take `abs()`.
pub fn cosine(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.len() != b.len() || a.is_empty() {
        return None;
    }
    let mut dot = 0.0_f64;
    let mut na = 0.0_f64;
    let mut nb = 0.0_f64;
    for (x, y) in a.iter().zip(b.iter()) {
        dot += x * y;
        na += x * x;
        nb += y * y;
    }
    if na == 0.0 || nb == 0.0 {
        return None;
    }
    Some(dot / (na.sqrt() * nb.sqrt()))
}

/// Indices of the top-`n` entries by absolute value, as a sorted set.
pub fn top_abs_indices(values: &[f64], top_n: usize) -> BTreeSet<usize> {
    let mut indexed: Vec<(usize, f64)> = values
        .iter()
        .enumerate()
        .map(|(i, v)| (i, v.abs()))
        .collect();
    indexed.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
    let keep = indexed.len().min(top_n);
    indexed.into_iter().take(keep).map(|(i, _)| i).collect()
}

/// Jaccard similarity of top-`n` absolute-value index sets.
pub fn jaccard_top_n(a: &[f64], b: &[f64], top_n: usize) -> f64 {
    let set_a = top_abs_indices(a, top_n);
    let set_b = top_abs_indices(b, top_n);
    if set_a.is_empty() && set_b.is_empty() {
        return 0.0;
    }
    let inter = set_a.iter().filter(|i| set_b.contains(i)).count();
    let union = set_a.len() + set_b.len() - inter;
    if union == 0 {
        0.0
    } else {
        inter as f64 / union as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mean_handles_empty() {
        assert_eq!(mean(&[]), 0.0);
        assert_eq!(mean(&[1.0, 2.0, 3.0]), 2.0);
    }

    #[test]
    fn ranks_handle_ties() {
        let r = ranks(&[1.0, 2.0, 2.0, 4.0]);
        assert_eq!(r, vec![1.0, 2.5, 2.5, 4.0]);
    }

    #[test]
    fn pearson_is_perfect_for_linear_input() {
        assert!((pearson(&[1.0, 2.0, 3.0], &[2.0, 4.0, 6.0]).unwrap() - 1.0).abs() < 1e-12);
        assert!(pearson(&[1.0, 1.0, 1.0], &[2.0, 3.0, 4.0]).is_none());
    }

    #[test]
    fn spearman_handles_monotonic() {
        assert!((spearman(&[1.0, 2.0, 3.0], &[10.0, 20.0, 30.0]).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn cosine_returns_signed_value() {
        assert!((cosine(&[1.0, 2.0, 3.0], &[-1.0, -2.0, -3.0]).unwrap() + 1.0).abs() < 1e-12);
    }

    #[test]
    fn jaccard_top_n_matches_hand_computed() {
        let a = vec![0.9, 0.1, -0.8, 0.2];
        let b = vec![0.85, 0.15, 0.0, -0.78];
        assert!((jaccard_top_n(&a, &b, 2) - 1.0 / 3.0).abs() < 1e-12);
    }
}
