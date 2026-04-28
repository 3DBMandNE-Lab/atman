//! Pre-ranked gene-set enrichment analysis (GSEA).
//!
//! Implements the weighted Kolmogorov-Smirnov-style enrichment statistic of
//! Subramanian et al. (2005) with gene-permutation null calibration. Walks
//! down a ranked gene list, incrementing a running sum by
//! `|stat[i]|^p / N_R` at positions whose gene is in the query set and
//! decrementing by `1 / (N - N_H)` at positions whose gene is not, where
//! `N_R = sum |stat[i]|^p` over set positions and `N_H = |set ∩ universe|`.
//!
//! The enrichment score (ES) is the signed maximum deviation from zero.
//! P-values are estimated by re-rolling membership masks under
//! gene-label-equivalent permutations (Korotkevich et al. 2019,
//! `fgseaSimple` formulation): for each permutation, choose `N_H` positions
//! uniformly at random and compute ES under that mask. The two formulations
//! are equivalent under exchangeability and the position-mask form avoids
//! re-ranking the full list per permutation.
//!
//! Matches `fgsea::fgseaSimple` ES on identical inputs within floating-point
//! precision. NES depends on the permutation draws and is reproducible
//! across atman re-runs at the same seed but is not byte-equal to fgsea's
//! NES (the two implementations use different PRNGs).
//!
//! All randomness comes from [`crate::ica::Xoshiro256pp`] seeded by the
//! caller; under a fixed seed, repeated runs produce identical output.

use crate::ica::Xoshiro256pp;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone)]
pub struct GseaConfig {
    pub n_permutations: usize,
    pub seed: u64,
    pub weight_p: f64,
    pub min_set_size: usize,
    pub max_set_size: usize,
}

impl Default for GseaConfig {
    fn default() -> Self {
        Self {
            n_permutations: 1000,
            seed: 42,
            weight_p: 1.0,
            min_set_size: 5,
            max_set_size: 500,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GseaResult {
    pub set_name: String,
    pub set_size: usize,
    pub es: f64,
    pub nes: f64,
    pub p_value: f64,
    pub leading_edge: Vec<String>,
}

/// Run GSEA on a pre-ranked gene list.
///
/// `ranked` is a list of `(gene, stat)` pairs sorted descending by `stat`
/// (ties broken however the caller chose). Genes appearing more than once
/// are de-duplicated by the caller; this function treats duplicates as
/// independent positions.
///
/// Returns one [`GseaResult`] per gene set whose `set ∩ ranked-universe`
/// size lies within `[config.min_set_size, config.max_set_size]`. Sets
/// outside the range are silently skipped (not paper-grade signal).
pub fn gsea(
    ranked: &[(String, f64)],
    gene_sets: &BTreeMap<String, BTreeSet<String>>,
    config: &GseaConfig,
) -> Vec<GseaResult> {
    if ranked.is_empty() || gene_sets.is_empty() {
        return Vec::new();
    }
    let n = ranked.len();
    let abs_scores: Vec<f64> = ranked
        .iter()
        .map(|(_, s)| s.abs().powf(config.weight_p))
        .collect();
    let gene_to_index: BTreeMap<&str, usize> = ranked
        .iter()
        .enumerate()
        .map(|(i, (g, _))| (g.as_str(), i))
        .collect();

    let mut results = Vec::new();
    let mut rng = Xoshiro256pp::new(config.seed);

    for (set_name, set_genes) in gene_sets {
        let positions: Vec<usize> = set_genes
            .iter()
            .filter_map(|g| gene_to_index.get(g.as_str()).copied())
            .collect();
        let n_h = positions.len();
        if n_h < config.min_set_size || n_h > config.max_set_size {
            continue;
        }
        let mut in_set = vec![false; n];
        for &p in &positions {
            in_set[p] = true;
        }
        let (es, peak_pos) = enrichment_score(&in_set, &abs_scores);

        // Leading edge: genes in the set whose position is at or before the
        // peak (ES > 0) or at or after the peak (ES < 0).
        let mut leading: Vec<&str> = if es >= 0.0 {
            positions
                .iter()
                .filter(|&&p| p <= peak_pos)
                .map(|&p| ranked[p].0.as_str())
                .collect()
        } else {
            positions
                .iter()
                .filter(|&&p| p >= peak_pos)
                .map(|&p| ranked[p].0.as_str())
                .collect()
        };
        leading.sort();

        // Permutation null: re-roll membership mask of size n_h, recompute ES.
        let mut perm_es: Vec<f64> = Vec::with_capacity(config.n_permutations);
        let mut perm_mask = vec![false; n];
        for _ in 0..config.n_permutations {
            sample_positions(&mut rng, n, n_h, &mut perm_mask);
            let (perm_es_val, _) = enrichment_score(&perm_mask, &abs_scores);
            perm_es.push(perm_es_val);
        }

        let (p_value, nes) = pvalue_and_nes(es, &perm_es);
        results.push(GseaResult {
            set_name: set_name.clone(),
            set_size: n_h,
            es,
            nes,
            p_value,
            leading_edge: leading.into_iter().map(str::to_string).collect(),
        });
    }
    results.sort_by(|a, b| {
        b.nes
            .abs()
            .partial_cmp(&a.nes.abs())
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.set_name.cmp(&b.set_name))
    });
    results
}

/// Weighted KS-style enrichment score and the position of the peak deviation.
///
/// At each position, the running sum increments by `abs[i] / N_R` if
/// `in_set[i]` else decrements by `1 / N_M`, where `N_R = sum abs[i] over
/// in_set positions` and `N_M = N - N_H`. Returns the signed maximum
/// deviation from zero and the (zero-based) position at which it occurs.
fn enrichment_score(in_set: &[bool], abs_scores: &[f64]) -> (f64, usize) {
    let n = in_set.len();
    let n_h = in_set.iter().filter(|&&b| b).count();
    if n_h == 0 || n_h == n {
        return (0.0, 0);
    }
    let n_r: f64 = abs_scores
        .iter()
        .zip(in_set)
        .filter_map(|(&s, &b)| if b { Some(s) } else { None })
        .sum();
    if n_r <= 0.0 {
        return (0.0, 0);
    }
    let miss_decrement = 1.0 / ((n - n_h) as f64);
    let mut running = 0.0_f64;
    let mut max_dev = 0.0_f64;
    let mut max_pos = 0usize;
    for i in 0..n {
        if in_set[i] {
            running += abs_scores[i] / n_r;
        } else {
            running -= miss_decrement;
        }
        if running.abs() > max_dev.abs() {
            max_dev = running;
            max_pos = i;
        }
    }
    (max_dev, max_pos)
}

/// Empirical two-sided p-value and normalized enrichment score from a
/// permutation null. NES = ES / mean(|ES_perm| with same sign as ES).
/// If no same-sign permutation exists, NES is `f64::NAN` and p-value is
/// `1 / (n_permutations + 1)` (one-sided lower bound).
fn pvalue_and_nes(es: f64, perm_es: &[f64]) -> (f64, f64) {
    if es == 0.0 || perm_es.is_empty() {
        return (1.0, 0.0);
    }
    let same_sign: Vec<f64> = perm_es
        .iter()
        .copied()
        .filter(|&v| v.signum() == es.signum())
        .collect();
    let p_value = if same_sign.is_empty() {
        1.0 / ((perm_es.len() + 1) as f64)
    } else {
        let extreme = same_sign
            .iter()
            .filter(|&&v| v.abs() >= es.abs())
            .count();
        ((extreme + 1) as f64) / ((same_sign.len() + 1) as f64)
    };
    let nes = if same_sign.is_empty() {
        f64::NAN
    } else {
        let mean_abs = same_sign.iter().map(|v| v.abs()).sum::<f64>() / (same_sign.len() as f64);
        if mean_abs > 0.0 {
            es / mean_abs
        } else {
            f64::NAN
        }
    };
    (p_value, nes)
}

/// Sample `k` distinct positions from `0..n` uniformly at random, writing a
/// boolean membership mask into `out`. Uses Floyd's algorithm: O(k) draws.
fn sample_positions(rng: &mut Xoshiro256pp, n: usize, k: usize, out: &mut [bool]) {
    debug_assert_eq!(out.len(), n);
    debug_assert!(k <= n);
    for slot in out.iter_mut() {
        *slot = false;
    }
    if k == 0 {
        return;
    }
    if k >= n {
        for slot in out.iter_mut() {
            *slot = true;
        }
        return;
    }
    // Floyd's combination algorithm: for j = n-k..n, draw r in 0..=j; if
    // already chosen, choose j instead. Yields a uniform k-subset.
    for j in (n - k)..n {
        let r = rng_index(rng, j + 1);
        if out[r] {
            out[j] = true;
        } else {
            out[r] = true;
        }
    }
}

/// Uniform integer in `0..bound` with rejection sampling to remove modulo
/// bias. Mirrors the convention used by `crate::ica_null::permutation` —
/// extracts entropy from one `next_normal()` draw via bit reinterpretation.
fn rng_index(rng: &mut Xoshiro256pp, bound: usize) -> usize {
    debug_assert!(bound > 0);
    let bound_u = bound as u64;
    loop {
        let v = rng.next_normal().to_bits();
        let limit = u64::MAX - u64::MAX % bound_u;
        if v < limit {
            return (v % bound_u) as usize;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ranked_with_planted_top(n: usize, k: usize) -> Vec<(String, f64)> {
        // `k` planted "responder" genes get the top-ranked positions.
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let stat = (n - i) as f64;
            let name = if i < k {
                format!("R{:04}", i)
            } else {
                format!("N{:04}", i)
            };
            out.push((name, stat));
        }
        out
    }

    fn responder_set(k: usize) -> BTreeMap<String, BTreeSet<String>> {
        let mut m = BTreeMap::new();
        let mut s = BTreeSet::new();
        for i in 0..k {
            s.insert(format!("R{:04}", i));
        }
        m.insert("responders".to_string(), s);
        m
    }

    #[test]
    fn empty_inputs_return_empty() {
        let cfg = GseaConfig::default();
        let r = gsea(&[], &BTreeMap::new(), &cfg);
        assert!(r.is_empty());
    }

    #[test]
    fn perfectly_top_ranked_set_yields_es_one() {
        // Set sits at the very top of the ranked list → ES → 1.0.
        let ranked = ranked_with_planted_top(200, 10);
        let sets = responder_set(10);
        let cfg = GseaConfig {
            n_permutations: 200,
            ..GseaConfig::default()
        };
        let out = gsea(&ranked, &sets, &cfg);
        assert_eq!(out.len(), 1);
        let r = &out[0];
        assert!(
            (r.es - 1.0).abs() < 1e-9,
            "ES should equal 1.0 for top-ranked set, got {}",
            r.es
        );
        assert_eq!(r.set_size, 10);
        assert_eq!(r.leading_edge.len(), 10);
    }

    #[test]
    fn perfectly_bottom_ranked_set_yields_negative_es() {
        // Place the set at the very bottom: ES should be negative ≈ -1.
        let mut ranked = Vec::new();
        for i in 0..200 {
            let stat = (200 - i) as f64;
            let name = if i < 190 {
                format!("N{:04}", i)
            } else {
                format!("R{:04}", i - 190)
            };
            ranked.push((name, stat));
        }
        let sets = responder_set(10);
        let cfg = GseaConfig {
            n_permutations: 200,
            ..GseaConfig::default()
        };
        let out = gsea(&ranked, &sets, &cfg);
        let r = &out[0];
        assert!(
            r.es < -0.99,
            "ES should be ≈ -1.0 for bottom-ranked set, got {}",
            r.es
        );
    }

    #[test]
    fn random_set_has_es_near_zero_and_nonsignificant_p() {
        // Randomly placed set members → ES near zero on average and a
        // non-significant p-value at 200 permutations.
        let n = 200;
        let mut ranked = Vec::with_capacity(n);
        for i in 0..n {
            ranked.push((format!("G{:04}", i), (n - i) as f64));
        }
        // Pick a deterministic non-clustered subset: every 20th gene.
        let mut s = BTreeSet::new();
        for i in (0..n).step_by(20) {
            s.insert(format!("G{:04}", i));
        }
        let mut sets = BTreeMap::new();
        sets.insert("scattered".to_string(), s);
        let cfg = GseaConfig {
            n_permutations: 500,
            min_set_size: 1,
            ..GseaConfig::default()
        };
        let out = gsea(&ranked, &sets, &cfg);
        let r = &out[0];
        assert!(
            r.es.abs() < 0.5,
            "ES should be near zero for scattered set, got {}",
            r.es
        );
        assert!(
            r.p_value > 0.05,
            "p-value should be non-significant for scattered set, got {}",
            r.p_value
        );
    }

    #[test]
    fn determinism_under_fixed_seed() {
        let ranked = ranked_with_planted_top(150, 8);
        let sets = responder_set(8);
        let cfg = GseaConfig {
            n_permutations: 300,
            seed: 17,
            ..GseaConfig::default()
        };
        let a = gsea(&ranked, &sets, &cfg);
        let b = gsea(&ranked, &sets, &cfg);
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.set_name, y.set_name);
            assert_eq!(x.es.to_bits(), y.es.to_bits());
            assert_eq!(x.nes.to_bits(), y.nes.to_bits());
            assert_eq!(x.p_value.to_bits(), y.p_value.to_bits());
            assert_eq!(x.leading_edge, y.leading_edge);
        }
    }

    #[test]
    fn set_size_filtering_skips_oversized_sets() {
        let ranked = ranked_with_planted_top(100, 50);
        let mut sets = BTreeMap::new();
        let mut big = BTreeSet::new();
        for i in 0..50 {
            big.insert(format!("R{:04}", i));
        }
        sets.insert("oversized".to_string(), big);
        let cfg = GseaConfig {
            n_permutations: 50,
            max_set_size: 30,
            min_set_size: 5,
            ..GseaConfig::default()
        };
        let out = gsea(&ranked, &sets, &cfg);
        assert!(out.is_empty(), "oversized set should be filtered");
    }

    #[test]
    fn enrichment_score_basic_walk() {
        // Mask: T T F F F F F F F F  with abs_scores all 1.0
        // n=10, n_h=2, miss = 1/8.
        // Walk: +0.5, +1.0 (peak), -0.125, ... -1.0 at end.
        let mask = vec![
            true, true, false, false, false, false, false, false, false, false,
        ];
        let scores = vec![1.0; 10];
        let (es, pos) = enrichment_score(&mask, &scores);
        assert!((es - 1.0).abs() < 1e-12);
        assert_eq!(pos, 1);
    }
}
