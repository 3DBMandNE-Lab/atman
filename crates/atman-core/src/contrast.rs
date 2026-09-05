//! Subject-level contrast statistics shared by `axes contrast`, `axes groups`,
//! `axes anchor`, `score weighted`, `concordance`, and the `de` OLS effect size.
//!
//! Every function is pure and deterministic; all dispersion estimates use
//! ddof = 1.
//!
//! - Cohen d: `(mean_a − mean_b) / s_pooled` with
//!   `s_pooled² = ((n_a−1)s_a² + (n_b−1)s_b²) / (n_a+n_b−2)`.
//! - d CI: `d ± 1.96·sqrt((n_a+n_b)/(n_a·n_b) + d²/(2(n_a+n_b)))`
//!   (Hedges & Olkin large-sample standard error).
//! - AUC: `U_a / (n_a·n_b)` with `U_a = R_a − n_a(n_a+1)/2` on midranks of the
//!   pooled sample (`scipy.stats.mannwhitneyu(a, b).statistic / (n_a n_b)`).
//! - Kruskal–Wallis: H with the tie correction, chi-square on k−1 df
//!   (`scipy.stats.kruskal`).
//! - Spearman p: `t = rho·sqrt((n−2)/(1−rho²))`, two-sided on n−2 df
//!   (`scipy.stats.spearmanr`).
//! - Percentile: linear interpolation on the sorted sample (numpy default,
//!   R type 7).

use crate::stats::{mean, ranks, spearman};
use statrs::distribution::{ChiSquared, ContinuousCDF, StudentsT};

/// Unbiased sample variance (ddof = 1). `None` when fewer than two values.
pub fn sample_variance(values: &[f64]) -> Option<f64> {
    if values.len() < 2 {
        return None;
    }
    let m = mean(values);
    Some(values.iter().map(|v| (v - m).powi(2)).sum::<f64>() / (values.len() as f64 - 1.0))
}

/// Sample standard deviation (ddof = 1).
pub fn sample_sd(values: &[f64]) -> Option<f64> {
    sample_variance(values).map(f64::sqrt)
}

/// Pooled-SD Cohen d, `a` minus `b`. `None` when either group has fewer than
/// two values or the pooled variance is not strictly positive.
pub fn cohen_d_pooled(a: &[f64], b: &[f64]) -> Option<f64> {
    let (na, nb) = (a.len(), b.len());
    if na < 2 || nb < 2 {
        return None;
    }
    let va = sample_variance(a)?;
    let vb = sample_variance(b)?;
    let pooled = ((na as f64 - 1.0) * va + (nb as f64 - 1.0) * vb) / (na + nb - 2) as f64;
    if !pooled.is_finite() || pooled <= 0.0 {
        return None;
    }
    Some((mean(a) - mean(b)) / pooled.sqrt())
}

/// Large-sample 95% CI for Cohen d.
pub fn cohen_d_ci(d: f64, n_a: usize, n_b: usize) -> Option<(f64, f64)> {
    if n_a == 0 || n_b == 0 || !d.is_finite() {
        return None;
    }
    let n = (n_a + n_b) as f64;
    let se = (n / (n_a as f64 * n_b as f64) + d * d / (2.0 * n)).sqrt();
    Some((d - 1.96 * se, d + 1.96 * se))
}

/// Area under the ROC curve for `a` versus `b` via the Mann–Whitney U of `a`.
pub fn auc_mann_whitney(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.is_empty() || b.is_empty() {
        return None;
    }
    let mut pooled = Vec::with_capacity(a.len() + b.len());
    pooled.extend_from_slice(a);
    pooled.extend_from_slice(b);
    if pooled.iter().any(|v| !v.is_finite()) {
        return None;
    }
    let r = ranks(&pooled);
    let r_a: f64 = r[..a.len()].iter().sum();
    let u_a = r_a - (a.len() as f64) * (a.len() as f64 + 1.0) / 2.0;
    Some(u_a / (a.len() as f64 * b.len() as f64))
}

/// Kruskal–Wallis H test result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KruskalWallis {
    pub h: f64,
    pub p_value: f64,
    pub df: usize,
}

/// Kruskal–Wallis across `groups` (empty groups are ignored). `None` when
/// fewer than two non-empty groups, any non-finite value, or all values tie.
pub fn kruskal_wallis(groups: &[Vec<f64>]) -> Option<KruskalWallis> {
    let mut pooled: Vec<f64> = Vec::new();
    let mut sizes: Vec<usize> = Vec::new();
    for g in groups {
        if g.is_empty() {
            continue;
        }
        if g.iter().any(|v| !v.is_finite()) {
            return None;
        }
        pooled.extend_from_slice(g);
        sizes.push(g.len());
    }
    let k = sizes.len();
    if k < 2 {
        return None;
    }
    let n = pooled.len() as f64;
    let r = ranks(&pooled);
    let mut offset = 0usize;
    let mut sum = 0.0_f64;
    for &sz in &sizes {
        let rs: f64 = r[offset..offset + sz].iter().sum();
        sum += rs * rs / sz as f64;
        offset += sz;
    }
    let mut h = 12.0 / (n * (n + 1.0)) * sum - 3.0 * (n + 1.0);
    let mut sorted = pooled.clone();
    sorted.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    let mut tie_sum = 0.0_f64;
    let mut i = 0usize;
    while i < sorted.len() {
        let mut j = i + 1;
        while j < sorted.len() && sorted[j] == sorted[i] {
            j += 1;
        }
        let t = (j - i) as f64;
        tie_sum += t * t * t - t;
        i = j;
    }
    let correction = 1.0 - tie_sum / (n * n * n - n);
    if !correction.is_finite() || correction <= 0.0 {
        return None;
    }
    h /= correction;
    let df = k - 1;
    let chi = ChiSquared::new(df as f64).ok()?;
    let p_value = chi.sf(h).clamp(0.0, 1.0);
    Some(KruskalWallis { h, p_value, df })
}

/// Spearman rho with the two-sided t-approximation p-value. `p` is `NaN`
/// when n < 3 and `0.0` when |rho| = 1.
pub fn spearman_with_p(xs: &[f64], ys: &[f64]) -> Option<(f64, f64)> {
    let rho = spearman(xs, ys)?;
    let n = xs.len();
    if n < 3 {
        return Some((rho, f64::NAN));
    }
    if rho.abs() >= 1.0 - 1e-12 {
        return Some((rho.clamp(-1.0, 1.0), 0.0));
    }
    let df = (n - 2) as f64;
    let t = rho * (df / (1.0 - rho * rho)).sqrt();
    let dist = StudentsT::new(0.0, 1.0, df).ok()?;
    Some((rho, (2.0 * dist.sf(t.abs())).clamp(0.0, 1.0)))
}

/// t-based confidence interval of the mean at `level` (e.g. 0.95).
pub fn t_ci_mean(values: &[f64], level: f64) -> Option<(f64, f64)> {
    let n = values.len();
    if n < 2 || !(0.0..1.0).contains(&level) {
        return None;
    }
    let sd = sample_sd(values)?;
    let m = mean(values);
    let dist = StudentsT::new(0.0, 1.0, (n - 1) as f64).ok()?;
    let h = dist.inverse_cdf(0.5 + level / 2.0) * sd / (n as f64).sqrt();
    Some((m - h, m + h))
}

/// Linear-interpolation percentile of an ascending-sorted slice, `q` in [0, 1].
pub fn percentile(sorted: &[f64], q: f64) -> Option<f64> {
    if sorted.is_empty() || !(0.0..=1.0).contains(&q) {
        return None;
    }
    let idx = q * (sorted.len() - 1) as f64;
    let lo = idx.floor() as usize;
    let hi = idx.ceil() as usize;
    let frac = idx - lo as f64;
    Some(sorted[lo] + frac * (sorted[hi] - sorted[lo]))
}

/// Standardize over the finite entries (ddof = 1); `None` entries stay
/// `None`; a constant column yields all `None`.
pub fn zscore(values: &[Option<f64>]) -> Vec<Option<f64>> {
    let finite: Vec<f64> = values
        .iter()
        .filter_map(|v| v.filter(|x| x.is_finite()))
        .collect();
    let sd = match sample_sd(&finite) {
        Some(s) if s > 0.0 && s.is_finite() => s,
        _ => return vec![None; values.len()],
    };
    let m = mean(&finite);
    values
        .iter()
        .map(|v| v.filter(|x| x.is_finite()).map(|x| (x - m) / sd))
        .collect()
}

/// Median; even lengths average the middle pair.
pub fn median(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    let mut v = values.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = v.len() / 2;
    Some(if v.len() % 2 == 1 {
        v[mid]
    } else {
        (v[mid - 1] + v[mid]) / 2.0
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn variance_and_sd_use_ddof_one() {
        assert!(close(
            sample_variance(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]).unwrap(),
            4.571428571428571,
            1e-12
        ));
        assert!(sample_variance(&[1.0]).is_none());
        assert!(close(
            sample_sd(&[1.0, 2.0, 3.0, 4.0]).unwrap(),
            1.2909944487358056,
            1e-12
        ));
    }

    #[test]
    fn cohen_d_pooled_matches_hand_computation() {
        let a = [2.0, 3.0, 4.0];
        let b = [0.0, 1.0, 2.0];
        assert!(close(cohen_d_pooled(&a, &b).unwrap(), 2.0, 1e-12));
        assert!(close(
            cohen_d_pooled(&[1.0, 2.0, 3.0, 4.0], &[4.0, 6.0]).unwrap(),
            -1.8898223650461363,
            1e-12
        ));
        assert!(cohen_d_pooled(&[1.0], &[2.0, 3.0]).is_none());
        assert!(cohen_d_pooled(&[1.0, 1.0], &[2.0, 2.0]).is_none());
    }

    #[test]
    fn cohen_d_ci_uses_large_sample_se() {
        let (lo, hi) = cohen_d_ci(0.5, 20, 30).unwrap();
        assert!(close(lo, 0.5 - 1.96 * 0.29297326, 1e-7));
        assert!(close(hi, 0.5 + 1.96 * 0.29297326, 1e-7));
        assert!(cohen_d_ci(0.5, 0, 3).is_none());
    }

    #[test]
    fn auc_handles_ties_with_midranks() {
        assert!(close(
            auc_mann_whitney(&[1.0, 2.0, 3.0], &[2.0, 3.0, 4.0]).unwrap(),
            2.0 / 9.0,
            1e-12
        ));
        assert!(close(
            auc_mann_whitney(&[5.0, 6.0], &[1.0, 2.0, 3.0]).unwrap(),
            1.0,
            1e-12
        ));
        assert!(auc_mann_whitney(&[], &[1.0]).is_none());
    }

    #[test]
    fn kruskal_wallis_matches_scipy() {
        let kw = kruskal_wallis(&[
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ])
        .unwrap();
        assert!(close(kw.h, 7.2, 1e-12));
        assert!(close(kw.p_value, 0.02732372244729256, 1e-10));
        assert_eq!(kw.df, 2);
        // scipy.stats.kruskal([1,1,2],[2,3,3]) -> H = 3.3333333333333295, p = 0.06788915486182918
        let kw = kruskal_wallis(&[vec![1.0, 1.0, 2.0], vec![2.0, 3.0, 3.0]]).unwrap();
        assert!(close(kw.h, 3.3333333333333295, 1e-10));
        assert!(close(kw.p_value, 0.06788915486182918, 1e-10));
        assert!(kruskal_wallis(&[vec![1.0, 2.0]]).is_none());
    }

    #[test]
    fn spearman_p_uses_t_approximation() {
        let (rho, p) =
            spearman_with_p(&[1.0, 2.0, 3.0, 4.0, 5.0], &[1.0, 3.0, 2.0, 5.0, 4.0]).unwrap();
        assert!(close(rho, 0.8, 1e-12));
        assert!(close(p, 0.10408803866182788, 1e-10));
        let (rho, p) = spearman_with_p(&[1.0, 2.0, 3.0], &[10.0, 20.0, 30.0]).unwrap();
        assert!(close(rho, 1.0, 1e-12));
        assert_eq!(p, 0.0);
    }

    #[test]
    fn t_ci_of_mean_matches_scipy() {
        let (lo, hi) = t_ci_mean(&[1.0, 2.0, 3.0, 4.0, 5.0], 0.95).unwrap();
        assert!(close(lo, 1.036756838522439, 1e-9), "lo = {lo}");
        assert!(close(hi, 4.9632431614775605, 1e-9), "hi = {hi}");
        assert!(t_ci_mean(&[1.0], 0.95).is_none());
    }

    #[test]
    fn percentile_interpolates_linearly() {
        let s = [1.0, 2.0, 3.0, 4.0];
        assert!(close(percentile(&s, 0.025).unwrap(), 1.075, 1e-12));
        assert!(close(percentile(&s, 0.5).unwrap(), 2.5, 1e-12));
        assert!(close(percentile(&s, 1.0).unwrap(), 4.0, 1e-12));
        assert!(percentile(&[], 0.5).is_none());
    }

    #[test]
    fn zscore_skips_missing_and_uses_ddof_one() {
        let z = zscore(&[Some(1.0), None, Some(2.0), Some(3.0)]);
        assert!(close(z[0].unwrap(), -1.0, 1e-12));
        assert!(z[1].is_none());
        assert!(close(z[2].unwrap(), 0.0, 1e-12));
        assert!(close(z[3].unwrap(), 1.0, 1e-12));
        assert!(zscore(&[Some(1.0), Some(1.0)]).iter().all(|v| v.is_none()));
    }

    #[test]
    fn median_averages_middle_pair() {
        assert_eq!(median(&[3.0, 1.0, 2.0]).unwrap(), 2.0);
        assert_eq!(median(&[4.0, 1.0, 3.0, 2.0]).unwrap(), 2.5);
        assert!(median(&[]).is_none());
    }
}
