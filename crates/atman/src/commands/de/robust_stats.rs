//! Distribution-agnostic companion statistics emitted alongside the
//! parametric DE result on every paired and unpaired path. Covers
//! effect size (Cohen's d_z / Hedges' g), mean-difference CI, Wilcoxon
//! signed-rank / rank-sum p-values, median difference, and 20%
//! trimmed-mean difference.

use atman_core::stats::mean;
use statrs::distribution::{ContinuousCDF, Normal, StudentsT};

#[derive(Debug, Clone, Default)]
pub(super) struct RobustStats {
    pub effect_size: Option<f64>,
    pub effect_size_method: String,
    pub ci_low: Option<f64>,
    pub ci_high: Option<f64>,
    pub wilcoxon_p: Option<f64>,
    pub wilcoxon_method: String,
    pub median_diff: Option<f64>,
    pub trimmed_mean_diff: Option<f64>,
}

pub(super) fn robust_paired(pairs: &[(f64, f64)], min_pairs: usize) -> RobustStats {
    let diffs: Vec<f64> = pairs.iter().map(|(a, b)| a - b).collect();
    if diffs.len() < min_pairs || diffs.len() < 2 || diffs.iter().any(|v| !v.is_finite()) {
        return RobustStats::default();
    }
    let mean_diff = mean(&diffs);
    let sd = sample_sd(&diffs, mean_diff);
    let (ci_low, ci_high) = mean_ci(mean_diff, sd, diffs.len());
    RobustStats {
        effect_size: if sd > 0.0 { Some(mean_diff / sd) } else { None },
        effect_size_method: "cohen_dz".to_string(),
        ci_low,
        ci_high,
        wilcoxon_p: wilcoxon_signed_rank_p(&diffs),
        wilcoxon_method: "signed_rank".to_string(),
        median_diff: median(diffs.clone()),
        trimmed_mean_diff: trimmed_mean(diffs, 0.2),
    }
}

pub(super) fn robust_unpaired(a: &[f64], b: &[f64], min_pairs: usize) -> RobustStats {
    if a.len() < min_pairs
        || b.len() < min_pairs
        || a.len() < 2
        || b.len() < 2
        || a.iter().chain(b.iter()).any(|v| !v.is_finite())
    {
        return RobustStats::default();
    }
    let mean_a = mean(a);
    let mean_b = mean(b);
    let var_a = sample_var(a, mean_a);
    let var_b = sample_var(b, mean_b);
    let mean_diff = mean_a - mean_b;
    let df_pooled = (a.len() + b.len() - 2) as f64;
    let pooled_var = ((a.len() - 1) as f64 * var_a + (b.len() - 1) as f64 * var_b) / df_pooled;
    let hedges_j = 1.0 - 3.0 / (4.0 * df_pooled - 1.0);
    let effect_size = if pooled_var > 0.0 {
        Some(hedges_j * mean_diff / pooled_var.sqrt())
    } else {
        None
    };
    let se = (var_a / a.len() as f64 + var_b / b.len() as f64).sqrt();
    let df = welch_df(var_a, var_b, a.len(), b.len());
    let (ci_low, ci_high) = if se > 0.0 && df.is_finite() && df > 0.0 {
        let dist = StudentsT::new(0.0, 1.0, df).ok();
        let crit = dist.map(|d| d.inverse_cdf(0.975));
        (
            crit.map(|c| mean_diff - c * se),
            crit.map(|c| mean_diff + c * se),
        )
    } else {
        (None, None)
    };
    RobustStats {
        effect_size,
        effect_size_method: "hedges_g".to_string(),
        ci_low,
        ci_high,
        wilcoxon_p: wilcoxon_rank_sum_p(a, b),
        wilcoxon_method: "rank_sum".to_string(),
        median_diff: median(a.to_vec())
            .zip(median(b.to_vec()))
            .map(|(ma, mb)| ma - mb),
        trimmed_mean_diff: trimmed_mean(a.to_vec(), 0.2)
            .zip(trimmed_mean(b.to_vec(), 0.2))
            .map(|(ma, mb)| ma - mb),
    }
}

pub(super) fn median(mut values: Vec<f64>) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = values.len() / 2;
    if values.len() % 2 == 1 {
        Some(values[mid])
    } else {
        Some((values[mid - 1] + values[mid]) / 2.0)
    }
}

fn sample_var(values: &[f64], mean: f64) -> f64 {
    values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (values.len() - 1) as f64
}

fn sample_sd(values: &[f64], mean: f64) -> f64 {
    sample_var(values, mean).sqrt()
}

fn mean_ci(mean: f64, sd: f64, n: usize) -> (Option<f64>, Option<f64>) {
    if n < 2 || sd == 0.0 {
        return (None, None);
    }
    let df = (n - 1) as f64;
    let Some(crit) = StudentsT::new(0.0, 1.0, df)
        .ok()
        .map(|d| d.inverse_cdf(0.975))
    else {
        return (None, None);
    };
    let se = sd / (n as f64).sqrt();
    (Some(mean - crit * se), Some(mean + crit * se))
}

fn welch_df(var_a: f64, var_b: f64, n_a: usize, n_b: usize) -> f64 {
    let va = var_a / n_a as f64;
    let vb = var_b / n_b as f64;
    let den = va.powi(2) / (n_a - 1) as f64 + vb.powi(2) / (n_b - 1) as f64;
    if den == 0.0 {
        f64::NAN
    } else {
        (va + vb).powi(2) / den
    }
}

fn trimmed_mean(mut values: Vec<f64>, proportion: f64) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let trim = ((values.len() as f64) * proportion).floor() as usize;
    if trim * 2 >= values.len() {
        return None;
    }
    Some(mean(&values[trim..(values.len() - trim)]))
}

fn wilcoxon_signed_rank_p(diffs: &[f64]) -> Option<f64> {
    let nonzero: Vec<f64> = diffs.iter().copied().filter(|d| *d != 0.0).collect();
    let n = nonzero.len();
    if n < 2 {
        return None;
    }
    let abs_values: Vec<f64> = nonzero.iter().map(|d| d.abs()).collect();
    let ranks = average_ranks(&abs_values);
    let w_plus: f64 = nonzero
        .iter()
        .zip(ranks.iter())
        .filter_map(|(d, r)| if *d > 0.0 { Some(*r) } else { None })
        .sum();
    let n_f = n as f64;
    let mean_w = n_f * (n_f + 1.0) / 4.0;
    let var_w = n_f * (n_f + 1.0) * (2.0 * n_f + 1.0) / 24.0;
    normal_two_sided_p(w_plus, mean_w, var_w)
}

fn wilcoxon_rank_sum_p(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.len() < 2 || b.len() < 2 {
        return None;
    }
    let mut values = Vec::with_capacity(a.len() + b.len());
    values.extend_from_slice(a);
    values.extend_from_slice(b);
    let ranks = average_ranks(&values);
    let rank_sum_a: f64 = ranks.iter().take(a.len()).sum();
    let n_a = a.len() as f64;
    let n_b = b.len() as f64;
    let u_a = rank_sum_a - n_a * (n_a + 1.0) / 2.0;
    let mean_u = n_a * n_b / 2.0;
    let var_u = n_a * n_b * (n_a + n_b + 1.0) / 12.0;
    normal_two_sided_p(u_a, mean_u, var_u)
}

fn normal_two_sided_p(stat: f64, mean: f64, var: f64) -> Option<f64> {
    if !(var.is_finite() && var > 0.0) {
        return None;
    }
    let z = (stat - mean).abs() / var.sqrt();
    let dist = Normal::new(0.0, 1.0).ok()?;
    Some(2.0 * (1.0 - dist.cdf(z)))
}

fn average_ranks(values: &[f64]) -> Vec<f64> {
    let mut indexed: Vec<(usize, f64)> = values.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let mut ranks = vec![0.0; values.len()];
    let mut i = 0;
    while i < indexed.len() {
        let mut j = i + 1;
        while j < indexed.len() && indexed[j].1 == indexed[i].1 {
            j += 1;
        }
        let avg_rank = (i + 1 + j) as f64 / 2.0;
        for k in i..j {
            ranks[indexed[k].0] = avg_rank;
        }
        i = j;
    }
    ranks
}
