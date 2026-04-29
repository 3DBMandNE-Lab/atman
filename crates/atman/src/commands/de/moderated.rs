//! Moderated (eBayes-style) variance shrinkage for the paired-t /
//! Welch-t paths. Shrinks each protein's observed variance toward
//! the median across proteins with a fixed prior degrees-of-freedom,
//! then recomputes the t-statistic and p-value.

use anyhow::{Context, Result};
use statrs::distribution::StudentsT;

use atman_core::de::two_sided_t_p_value;

use crate::io::DeResultRow;

pub(super) fn apply_moderated_shrinkage(rows: &mut [DeResultRow], prior_df: f64) -> Result<()> {
    let mut vars: Vec<f64> = rows
        .iter()
        .filter_map(paired_variance_from_row)
        .filter(|v| v.is_finite() && *v > 0.0)
        .collect();
    if vars.is_empty() {
        return Ok(());
    }
    vars.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let prior_var = if vars.len() % 2 == 1 {
        vars[vars.len() / 2]
    } else {
        (vars[vars.len() / 2 - 1] + vars[vars.len() / 2]) / 2.0
    };

    for row in rows.iter_mut() {
        let mean_diff = match row.mean_diff {
            Some(v) => v,
            None => continue,
        };
        let df_i = match row.df {
            Some(v) if v > 0.0 => v,
            _ => continue,
        };
        let n = row.n_pairs as f64;
        if n < 2.0 {
            continue;
        }
        let var_i = match paired_variance_from_row(row) {
            Some(v) if v.is_finite() && v >= 0.0 => v,
            _ => continue,
        };
        let post_var = (prior_df * prior_var + df_i * var_i) / (prior_df + df_i);
        if !(post_var.is_finite() && post_var > 0.0) {
            continue;
        }
        let t_mod = mean_diff / (post_var / n).sqrt();
        let df_mod = prior_df + df_i;
        StudentsT::new(0.0, 1.0, df_mod)
            .with_context(|| format!("building t distribution with df={}", df_mod))?;
        let p_mod = two_sided_t_p_value(t_mod, df_mod).unwrap_or(f64::NAN);
        if p_mod.is_finite() {
            row.t = Some(t_mod);
            row.df = Some(df_mod);
            row.p_value = Some(p_mod);
        }
    }
    Ok(())
}

fn paired_variance_from_row(row: &DeResultRow) -> Option<f64> {
    let t = row.t?;
    let mean = row.mean_diff?;
    let n = row.n_pairs as f64;
    if n < 2.0 {
        return None;
    }
    if t == 0.0 {
        return Some(0.0);
    }
    let se = mean / t;
    Some((se * n.sqrt()).powi(2))
}
