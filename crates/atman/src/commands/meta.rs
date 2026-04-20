use anyhow::{bail, Context, Result};
use atman_core::de::bh_fdr;
use clap::{Args as ClapArgs, ValueEnum};
use csv::ReaderBuilder;
use serde_json::json;
use statrs::distribution::{ContinuousCDF, Normal};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{
    atomic_write, hash_labeled_inputs, need_col, sidecar_path_for, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Comma-separated de_results.tsv paths.
    #[arg(long)]
    inputs: String,

    /// Output TSV path.
    #[arg(long)]
    output: PathBuf,

    /// Meta-analysis method. Initial implementation reports all summaries;
    /// this flag is retained for CLI stability.
    #[arg(long, default_value = "fixed-effect")]
    method: String,

    /// Input level to meta-analyze.
    #[arg(long, value_enum, default_value_t = Level::Protein)]
    level: Level,

    /// Optional report mode for module-level inputs.
    #[arg(long)]
    report: Option<String>,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum Level {
    Protein,
    Module,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Key {
    comparison: String,
    panel: String,
    gene_symbol: String,
}

#[derive(Debug, Clone)]
struct CohortEffect {
    effect: f64,
    se: f64,
    p_value: f64,
}

#[derive(Debug)]
struct MetaRow {
    key: Key,
    n_cohorts: usize,
    fixed_effect: f64,
    fixed_se: f64,
    fixed_z: f64,
    fixed_p: f64,
    random_effect: f64,
    random_se: f64,
    random_z: f64,
    random_p: f64,
    tau2: f64,
    q_heterogeneity: f64,
    i2: f64,
    stouffer_z: f64,
    stouffer_p: f64,
    sign_consistency: f64,
    bh_q: Option<f64>,
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    if args.method != "fixed-effect" && args.method != "random-effects" && args.method != "all" {
        bail!("meta supports --method fixed-effect, random-effects, or all");
    }
    if let Some(report) = &args.report {
        if report != "sign-consistency" {
            bail!("meta supports --report sign-consistency");
        }
    }
    let inputs: Vec<PathBuf> = args
        .inputs
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(PathBuf::from)
        .collect();
    if inputs.len() < 2 {
        bail!("meta requires at least two input de_results.tsv files");
    }
    if args.level == Level::Module {
        run_module_meta(&inputs, &args.output, &args.method)?;
        let finished_at = SystemTime::now();
        write_meta_sidecar(&args, &inputs, started_at, finished_at)?;
        return Ok(());
    }
    let mut effects: BTreeMap<Key, Vec<CohortEffect>> = BTreeMap::new();
    for path in &inputs {
        for (key, effect) in read_effects(path)? {
            effects.entry(key).or_default().push(effect);
        }
    }
    let mut rows = Vec::new();
    let mut p_values = Vec::new();
    for (key, cohort_effects) in effects {
        if cohort_effects.len() < 2 {
            continue;
        }
        let row = summarize(key, &cohort_effects)?;
        p_values.push(Some(row.fixed_p));
        rows.push(row);
    }
    let qs = bh_fdr(&p_values);
    for (row, q) in rows.iter_mut().zip(qs) {
        row.bh_q = q;
    }
    rows.sort_by(|a, b| {
        match (a.bh_q, b.bh_q) {
            (Some(qa), Some(qb)) => qa.partial_cmp(&qb).unwrap_or(std::cmp::Ordering::Equal),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => std::cmp::Ordering::Equal,
        }
        .then_with(|| a.key.comparison.cmp(&b.key.comparison))
        .then_with(|| a.key.panel.cmp(&b.key.panel))
        .then_with(|| a.key.gene_symbol.cmp(&b.key.gene_symbol))
    });
    write_rows(&args.output, &rows)?;
    eprintln!(
        "meta: inputs={} rows={} method={}",
        inputs.len(),
        rows.len(),
        args.method
    );

    let finished_at = SystemTime::now();
    write_meta_sidecar(&args, &inputs, started_at, finished_at)?;
    Ok(())
}

fn write_meta_sidecar(
    args: &Args,
    inputs: &[PathBuf],
    started_at: SystemTime,
    finished_at: SystemTime,
) -> Result<()> {
    let labeled: Vec<(String, PathBuf)> = inputs
        .iter()
        .enumerate()
        .map(|(i, p)| (format!("input_{:02}", i + 1), p.clone()))
        .collect();
    let refs: Vec<(&str, &Path)> = labeled
        .iter()
        .map(|(l, p)| (l.as_str(), p.as_path()))
        .collect();
    let inputs_sha256 = hash_labeled_inputs(&refs)?;
    let sidecar = sidecar_path_for(&args.output);
    let level = match args.level {
        Level::Protein => "protein",
        Level::Module => "module",
    };
    write_run_sidecar(
        &sidecar,
        "meta",
        json!({
            "inputs": args.inputs,
            "output": args.output.display().to_string(),
            "method": args.method,
            "level": level,
            "report": args.report,
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
)?;
    eprintln!("meta: sidecar={}", sidecar.display());
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct ModuleKey {
    program: String,
    category: String,
}

#[derive(Debug)]
struct ModuleRow {
    key: ModuleKey,
    n_cohorts: usize,
    fixed_effect: f64,
    fixed_se: f64,
    fixed_z: f64,
    fixed_p: f64,
    random_effect: f64,
    random_se: f64,
    random_z: f64,
    random_p: f64,
    tau2: f64,
    q_heterogeneity: f64,
    i2: f64,
    sign_consistency: f64,
    n_positive: usize,
    n_negative: usize,
    sign_binomial_p: f64,
    bh_q: Option<f64>,
}

fn run_module_meta(inputs: &[PathBuf], output: &Path, method: &str) -> Result<()> {
    let mut effects: BTreeMap<ModuleKey, Vec<CohortEffect>> = BTreeMap::new();
    for path in inputs {
        for (key, effect) in read_module_effects(path)? {
            effects.entry(key).or_default().push(effect);
        }
    }
    let mut rows = Vec::new();
    let mut p_values = Vec::new();
    for (key, cohort_effects) in effects {
        if cohort_effects.len() < 2 {
            continue;
        }
        let row = summarize_module(key, &cohort_effects)?;
        p_values.push(Some(row.fixed_p));
        rows.push(row);
    }
    let qs = bh_fdr(&p_values);
    for (row, q) in rows.iter_mut().zip(qs) {
        row.bh_q = q;
    }
    rows.sort_by(|a, b| {
        match (a.bh_q, b.bh_q) {
            (Some(qa), Some(qb)) => qa.partial_cmp(&qb).unwrap_or(std::cmp::Ordering::Equal),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => std::cmp::Ordering::Equal,
        }
        .then_with(|| a.key.program.cmp(&b.key.program))
        .then_with(|| a.key.category.cmp(&b.key.category))
    });
    write_module_rows(output, &rows)?;
    eprintln!(
        "meta: inputs={} rows={} method={} level=module",
        inputs.len(),
        rows.len(),
        method
    );
    Ok(())
}

fn read_effects(path: &Path) -> Result<Vec<(Key, CohortEffect)>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let comparison_col = need_col(&headers, "comparison", path)?;
    let panel_col = need_col(&headers, "panel", path)?;
    let gene_col = need_col(&headers, "gene_symbol", path)?;
    let effect_col = need_col(&headers, "mean_diff", path)?;
    let t_col = need_col(&headers, "t", path)?;
    let p_col = need_col(&headers, "p_value", path)?;
    let mut out = Vec::new();
    for row in reader.records() {
        let row = row?;
        let effect = match parse_f64(row[effect_col].trim()) {
            Some(v) => v,
            None => continue,
        };
        let t = match parse_f64(row[t_col].trim()) {
            Some(v) if v != 0.0 => v,
            _ => continue,
        };
        let p_value = match parse_f64(row[p_col].trim()) {
            Some(v) => v.clamp(1e-300, 1.0),
            None => continue,
        };
        let se = (effect / t).abs();
        if !(se.is_finite() && se > 0.0) {
            continue;
        }
        out.push((
            Key {
                comparison: row[comparison_col].to_string(),
                panel: row[panel_col].to_string(),
                gene_symbol: row[gene_col].to_string(),
            },
            CohortEffect {
                effect,
                se,
                p_value,
            },
        ));
    }
    Ok(out)
}

fn read_module_effects(path: &Path) -> Result<Vec<(ModuleKey, CohortEffect)>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let program_col = need_col(&headers, "program", path)?;
    let category_col = headers.iter().position(|h| h == "category");
    let effect_col = need_col(&headers, "point_delta", path)?;
    let ci_lo_col = need_col(&headers, "ci_lo", path)?;
    let ci_hi_col = need_col(&headers, "ci_hi", path)?;
    let p_col = headers.iter().position(|h| h == "p_sign_stable");
    let mut out = Vec::new();
    for row in reader.records() {
        let row = row?;
        let effect = match parse_f64(row[effect_col].trim()) {
            Some(v) => v,
            None => continue,
        };
        let ci_lo = match parse_f64(row[ci_lo_col].trim()) {
            Some(v) => v,
            None => continue,
        };
        let ci_hi = match parse_f64(row[ci_hi_col].trim()) {
            Some(v) => v,
            None => continue,
        };
        let se = ((ci_hi - ci_lo).abs() / (2.0 * 1.959_963_984_540_054)).max(1e-12);
        let p_value = p_col
            .and_then(|col| parse_f64(row[col].trim()))
            .unwrap_or_else(|| {
                let normal = Normal::new(0.0, 1.0).expect("normal");
                two_sided_normal_p(&normal, effect / se)
            })
            .clamp(1e-300, 1.0);
        out.push((
            ModuleKey {
                program: row[program_col].to_string(),
                category: category_col
                    .map(|col| row[col].to_string())
                    .unwrap_or_default(),
            },
            CohortEffect {
                effect,
                se,
                p_value,
            },
        ));
    }
    Ok(out)
}

fn summarize(key: Key, effects: &[CohortEffect]) -> Result<MetaRow> {
    let normal = Normal::new(0.0, 1.0).context("building normal distribution")?;
    let weights: Vec<f64> = effects.iter().map(|e| 1.0 / e.se.powi(2)).collect();
    let sum_w: f64 = weights.iter().sum();
    let fixed_effect = weighted_mean(effects, &weights);
    let fixed_se = (1.0 / sum_w).sqrt();
    let fixed_z = fixed_effect / fixed_se;
    let fixed_p = two_sided_normal_p(&normal, fixed_z);
    let q_heterogeneity: f64 = effects
        .iter()
        .zip(weights.iter())
        .map(|(e, w)| w * (e.effect - fixed_effect).powi(2))
        .sum();
    let df = effects.len() as f64 - 1.0;
    let c = sum_w - weights.iter().map(|w| w.powi(2)).sum::<f64>() / sum_w;
    let tau2 = if c > 0.0 {
        ((q_heterogeneity - df) / c).max(0.0)
    } else {
        0.0
    };
    let random_weights: Vec<f64> = effects
        .iter()
        .map(|e| 1.0 / (e.se.powi(2) + tau2))
        .collect();
    let sum_wr: f64 = random_weights.iter().sum();
    let random_effect = weighted_mean(effects, &random_weights);
    let random_se = (1.0 / sum_wr).sqrt();
    let random_z = random_effect / random_se;
    let random_p = two_sided_normal_p(&normal, random_z);
    let i2 = if q_heterogeneity > df {
        ((q_heterogeneity - df) / q_heterogeneity).max(0.0)
    } else {
        0.0
    };
    let stouffer_z = stouffer_z(&normal, effects);
    let stouffer_p = two_sided_normal_p(&normal, stouffer_z);
    let positive = effects.iter().filter(|e| e.effect > 0.0).count();
    let negative = effects.iter().filter(|e| e.effect < 0.0).count();
    let sign_consistency = positive.max(negative) as f64 / effects.len() as f64;
    Ok(MetaRow {
        key,
        n_cohorts: effects.len(),
        fixed_effect,
        fixed_se,
        fixed_z,
        fixed_p,
        random_effect,
        random_se,
        random_z,
        random_p,
        tau2,
        q_heterogeneity,
        i2,
        stouffer_z,
        stouffer_p,
        sign_consistency,
        bh_q: None,
    })
}

fn summarize_module(key: ModuleKey, effects: &[CohortEffect]) -> Result<ModuleRow> {
    let normal = Normal::new(0.0, 1.0).context("building normal distribution")?;
    let weights: Vec<f64> = effects.iter().map(|e| 1.0 / e.se.powi(2)).collect();
    let sum_w: f64 = weights.iter().sum();
    let fixed_effect = weighted_mean(effects, &weights);
    let fixed_se = (1.0 / sum_w).sqrt();
    let fixed_z = fixed_effect / fixed_se;
    let fixed_p = two_sided_normal_p(&normal, fixed_z);
    let q_heterogeneity: f64 = effects
        .iter()
        .zip(weights.iter())
        .map(|(e, w)| w * (e.effect - fixed_effect).powi(2))
        .sum();
    let df = effects.len() as f64 - 1.0;
    let c = sum_w - weights.iter().map(|w| w.powi(2)).sum::<f64>() / sum_w;
    let tau2 = if c > 0.0 {
        ((q_heterogeneity - df) / c).max(0.0)
    } else {
        0.0
    };
    let random_weights: Vec<f64> = effects
        .iter()
        .map(|e| 1.0 / (e.se.powi(2) + tau2))
        .collect();
    let sum_wr: f64 = random_weights.iter().sum();
    let random_effect = weighted_mean(effects, &random_weights);
    let random_se = (1.0 / sum_wr).sqrt();
    let random_z = random_effect / random_se;
    let random_p = two_sided_normal_p(&normal, random_z);
    let i2 = if q_heterogeneity > df {
        ((q_heterogeneity - df) / q_heterogeneity).max(0.0)
    } else {
        0.0
    };
    let n_positive = effects.iter().filter(|e| e.effect > 0.0).count();
    let n_negative = effects.iter().filter(|e| e.effect < 0.0).count();
    let sign_consistency = n_positive.max(n_negative) as f64 / effects.len() as f64;
    let sign_binomial_p = binomial_upper_tail(n_positive.max(n_negative), effects.len());
    Ok(ModuleRow {
        key,
        n_cohorts: effects.len(),
        fixed_effect,
        fixed_se,
        fixed_z,
        fixed_p,
        random_effect,
        random_se,
        random_z,
        random_p,
        tau2,
        q_heterogeneity,
        i2,
        sign_consistency,
        n_positive,
        n_negative,
        sign_binomial_p,
        bh_q: None,
    })
}

fn weighted_mean(effects: &[CohortEffect], weights: &[f64]) -> f64 {
    effects
        .iter()
        .zip(weights.iter())
        .map(|(e, w)| e.effect * w)
        .sum::<f64>()
        / weights.iter().sum::<f64>()
}

fn stouffer_z(normal: &Normal, effects: &[CohortEffect]) -> f64 {
    let z_sum: f64 = effects
        .iter()
        .map(|e| {
            let z_abs = normal.inverse_cdf(1.0 - e.p_value / 2.0);
            if e.effect < 0.0 {
                -z_abs
            } else {
                z_abs
            }
        })
        .sum();
    z_sum / (effects.len() as f64).sqrt()
}

fn two_sided_normal_p(normal: &Normal, z: f64) -> f64 {
    2.0 * (1.0 - normal.cdf(z.abs()))
}

fn parse_f64(input: &str) -> Option<f64> {
    let value = input.parse::<f64>().ok()?;
    value.is_finite().then_some(value)
}


fn write_rows(path: &Path, rows: &[MetaRow]) -> Result<()> {
    let mut out = String::from(
        "comparison\tpanel\tgene_symbol\tn_cohorts\tfixed_effect\tfixed_se\tfixed_z\tfixed_p\tfixed_bh_q\trandom_effect\trandom_se\trandom_z\trandom_p\ttau2\tq_heterogeneity\ti2\tstouffer_z\tstouffer_p\tsign_consistency\n",
    );
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.key.comparison,
            r.key.panel,
            r.key.gene_symbol,
            r.n_cohorts,
            r.fixed_effect,
            r.fixed_se,
            r.fixed_z,
            r.fixed_p,
            r.bh_q.map(|q| q.to_string()).unwrap_or_default(),
            r.random_effect,
            r.random_se,
            r.random_z,
            r.random_p,
            r.tau2,
            r.q_heterogeneity,
            r.i2,
            r.stouffer_z,
            r.stouffer_p,
            r.sign_consistency,
        ));
    }
    atomic_write(path, out.as_bytes())
}

fn write_module_rows(path: &Path, rows: &[ModuleRow]) -> Result<()> {
    let mut out = String::from(
        "program\tcategory\tn_cohorts\tfixed_effect\tfixed_se\tfixed_z\tfixed_p\tfixed_bh_q\trandom_effect\trandom_se\trandom_z\trandom_p\ttau2\tq_heterogeneity\ti2\tsign_consistency\tn_positive\tn_negative\tsign_binomial_p\n",
    );
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.key.program,
            r.key.category,
            r.n_cohorts,
            r.fixed_effect,
            r.fixed_se,
            r.fixed_z,
            r.fixed_p,
            r.bh_q.map(|q| q.to_string()).unwrap_or_default(),
            r.random_effect,
            r.random_se,
            r.random_z,
            r.random_p,
            r.tau2,
            r.q_heterogeneity,
            r.i2,
            r.sign_consistency,
            r.n_positive,
            r.n_negative,
            r.sign_binomial_p,
        ));
    }
    atomic_write(path, out.as_bytes())
}

fn binomial_upper_tail(k: usize, n: usize) -> f64 {
    if n == 0 || k > n {
        return f64::NAN;
    }
    let denom = 2_f64.powi(n as i32);
    (k..=n).map(|i| combination(n, i) / denom).sum()
}

fn combination(n: usize, k: usize) -> f64 {
    let k = k.min(n - k);
    if k == 0 {
        return 1.0;
    }
    (1..=k).fold(1.0, |acc, i| acc * (n + 1 - i) as f64 / i as f64)
}
