use anyhow::{bail, Context, Result};
use clap::{Args as ClapArgs, ValueEnum};
use csv::ReaderBuilder;
use statrs::distribution::{ContinuousCDF, Normal};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use crate::io::atomic_write;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Long TSV with cohort, subject_id, program, and activation columns.
    #[arg(long)]
    activations: PathBuf,

    /// Pair manifest TSV with cohort, pair_label, progs_a, progs_b, hypothesized_sign.
    #[arg(long)]
    pairs: PathBuf,

    /// Correlation method.
    #[arg(long, value_enum, default_value_t = Method::Spearman)]
    method: Method,

    /// Report type. Initial implementation supports sign-consistency.
    #[arg(long, default_value = "sign-consistency")]
    report: String,

    /// Per-pair coupling output TSV.
    #[arg(long)]
    output: PathBuf,

    /// Sign-consistency summary output TSV.
    #[arg(long)]
    output_summary: PathBuf,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum Method {
    Spearman,
}

#[derive(Debug, Clone)]
struct Activation {
    subject_id: String,
    program: String,
    activation: f64,
}

#[derive(Debug, Clone)]
struct PairSpec {
    cohort: String,
    pair_label: String,
    progs_a: Vec<String>,
    progs_b: Vec<String>,
    hypothesized_sign: Sign,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Sign {
    Positive,
    Negative,
}

#[derive(Debug, Clone)]
struct CouplingRow {
    cohort: String,
    pair_label: String,
    progs_a: String,
    progs_b: String,
    hypothesized_sign: Sign,
    n_subjects: usize,
    rho: f64,
    ci_lo: f64,
    ci_hi: f64,
    p_value: f64,
    sign_match: bool,
}

pub fn run(args: Args) -> Result<()> {
    if args.report != "sign-consistency" {
        bail!("coupling supports --report sign-consistency");
    }
    if args.method != Method::Spearman {
        bail!("coupling currently supports --method spearman");
    }

    let activations = read_activations(&args.activations)?;
    let pairs = read_pairs(&args.pairs)?;
    let mut rows = Vec::new();
    for pair in &pairs {
        let cohort_activations = activations
            .get(&pair.cohort)
            .with_context(|| format!("no activations found for cohort {:?}", pair.cohort))?;
        rows.push(compute_pair(pair, cohort_activations)?);
    }

    write_coupling_rows(&args.output, &rows)?;
    write_summary(&args.output_summary, &rows)?;
    eprintln!(
        "coupling: pairs={} output={} summary={}",
        rows.len(),
        args.output.display(),
        args.output_summary.display()
    );
    Ok(())
}

fn read_activations(path: &Path) -> Result<BTreeMap<String, Vec<Activation>>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let cohort_col = need_col(&headers, "cohort", path)?;
    let subject_col = need_col(&headers, "subject_id", path)?;
    let program_col = need_col(&headers, "program", path)?;
    let activation_col = need_col(&headers, "activation", path)?;
    let mut by_cohort: BTreeMap<String, Vec<Activation>> = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let activation: f64 = row[activation_col]
            .parse()
            .with_context(|| format!("parsing activation in {:?}", path))?;
        if !activation.is_finite() {
            continue;
        }
        by_cohort
            .entry(row[cohort_col].to_string())
            .or_default()
            .push(Activation {
                subject_id: row[subject_col].to_string(),
                program: row[program_col].to_string(),
                activation,
            });
    }
    Ok(by_cohort)
}

fn read_pairs(path: &Path) -> Result<Vec<PairSpec>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let cohort_col = need_col(&headers, "cohort", path)?;
    let label_col = need_col(&headers, "pair_label", path)?;
    let progs_a_col = need_col(&headers, "progs_a", path)?;
    let progs_b_col = need_col(&headers, "progs_b", path)?;
    let sign_col = need_col(&headers, "hypothesized_sign", path)?;
    let mut out = Vec::new();
    for row in reader.records() {
        let row = row?;
        out.push(PairSpec {
            cohort: row[cohort_col].to_string(),
            pair_label: row[label_col].to_string(),
            progs_a: parse_programs(&row[progs_a_col])?,
            progs_b: parse_programs(&row[progs_b_col])?,
            hypothesized_sign: parse_sign(&row[sign_col])?,
        });
    }
    if out.is_empty() {
        bail!("no coupling pairs found in {:?}", path);
    }
    Ok(out)
}

fn compute_pair(pair: &PairSpec, activations: &[Activation]) -> Result<CouplingRow> {
    let set_a: BTreeSet<&str> = pair.progs_a.iter().map(String::as_str).collect();
    let set_b: BTreeSet<&str> = pair.progs_b.iter().map(String::as_str).collect();
    let mut per_subject: BTreeMap<&str, (Vec<f64>, Vec<f64>)> = BTreeMap::new();
    for activation in activations {
        let entry = per_subject
            .entry(&activation.subject_id)
            .or_insert_with(|| (Vec::new(), Vec::new()));
        if set_a.contains(activation.program.as_str()) {
            entry.0.push(activation.activation);
        }
        if set_b.contains(activation.program.as_str()) {
            entry.1.push(activation.activation);
        }
    }

    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for (_subject, (a_values, b_values)) in per_subject {
        if a_values.is_empty() || b_values.is_empty() {
            continue;
        }
        xs.push(mean(&a_values));
        ys.push(mean(&b_values));
    }
    if xs.len() < 4 {
        bail!(
            "coupling pair {:?} in cohort {:?} has {} complete subjects; need at least 4",
            pair.pair_label,
            pair.cohort,
            xs.len()
        );
    }

    let rho = spearman(&xs, &ys)?;
    let (ci_lo, ci_hi, p_value) = fisher_interval_and_p(rho, xs.len())?;
    let sign_match = match pair.hypothesized_sign {
        Sign::Positive => rho > 0.0,
        Sign::Negative => rho < 0.0,
    };
    Ok(CouplingRow {
        cohort: pair.cohort.clone(),
        pair_label: pair.pair_label.clone(),
        progs_a: pair.progs_a.join(";"),
        progs_b: pair.progs_b.join(";"),
        hypothesized_sign: pair.hypothesized_sign,
        n_subjects: xs.len(),
        rho,
        ci_lo,
        ci_hi,
        p_value,
        sign_match,
    })
}

fn spearman(xs: &[f64], ys: &[f64]) -> Result<f64> {
    if xs.len() != ys.len() {
        bail!("spearman inputs have different lengths");
    }
    pearson(&ranks(xs), &ranks(ys))
}

fn ranks(values: &[f64]) -> Vec<f64> {
    let mut indexed: Vec<(usize, f64)> = values.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let mut ranks = vec![0.0; values.len()];
    let mut i = 0;
    while i < indexed.len() {
        let mut j = i + 1;
        while j < indexed.len() && indexed[j].1 == indexed[i].1 {
            j += 1;
        }
        let rank = (i + 1 + j) as f64 / 2.0;
        for k in i..j {
            ranks[indexed[k].0] = rank;
        }
        i = j;
    }
    ranks
}

fn pearson(xs: &[f64], ys: &[f64]) -> Result<f64> {
    let mean_x = mean(xs);
    let mean_y = mean(ys);
    let mut sxx = 0.0;
    let mut syy = 0.0;
    let mut sxy = 0.0;
    for (x, y) in xs.iter().zip(ys.iter()) {
        let dx = x - mean_x;
        let dy = y - mean_y;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    if sxx == 0.0 || syy == 0.0 {
        bail!("correlation undefined for constant activation vector");
    }
    Ok((sxy / (sxx.sqrt() * syy.sqrt())).clamp(-1.0, 1.0))
}

fn fisher_interval_and_p(rho: f64, n: usize) -> Result<(f64, f64, f64)> {
    let normal = Normal::new(0.0, 1.0).context("building normal distribution")?;
    let bounded = rho.clamp(-0.999_999, 0.999_999);
    let z = 0.5 * ((1.0 + bounded) / (1.0 - bounded)).ln();
    let se = 1.0 / ((n as f64) - 3.0).sqrt();
    let zcrit = normal.inverse_cdf(0.975);
    let lo = (z - zcrit * se).tanh();
    let hi = (z + zcrit * se).tanh();
    let p = 2.0 * (1.0 - normal.cdf((z / se).abs()));
    Ok((lo, hi, p))
}

fn write_coupling_rows(path: &Path, rows: &[CouplingRow]) -> Result<()> {
    let mut out = String::from(
        "cohort\tpair_label\tprogs_a\tprogs_b\thypothesized_sign\tn_subjects\tspearman_rho\tci_lo\tci_hi\tp_value\tsign_match\n",
    );
    for row in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            row.cohort,
            row.pair_label,
            row.progs_a,
            row.progs_b,
            row.hypothesized_sign.as_str(),
            row.n_subjects,
            row.rho,
            row.ci_lo,
            row.ci_hi,
            row.p_value,
            row.sign_match as u8
        ));
    }
    atomic_write(path, out.as_bytes())
}

fn write_summary(path: &Path, rows: &[CouplingRow]) -> Result<()> {
    let mut out = String::from(
        "scope\tpair_label\tn_tests\tn_sign_matches\tn_positive_observed\tn_negative_observed\tbinomial_p_one_sided\n",
    );
    write_summary_row(&mut out, "all", "all", rows);
    let mut by_label: BTreeMap<&str, Vec<&CouplingRow>> = BTreeMap::new();
    for row in rows {
        by_label.entry(&row.pair_label).or_default().push(row);
    }
    for (label, label_rows) in by_label {
        let owned: Vec<CouplingRow> = label_rows.into_iter().cloned().collect();
        write_summary_row(&mut out, "pair_label", label, &owned);
    }
    atomic_write(path, out.as_bytes())
}

fn write_summary_row(out: &mut String, scope: &str, label: &str, rows: &[CouplingRow]) {
    let n = rows.len();
    let matches = rows.iter().filter(|row| row.sign_match).count();
    let positive = rows.iter().filter(|row| row.rho > 0.0).count();
    let negative = rows.iter().filter(|row| row.rho < 0.0).count();
    let p = binomial_upper_tail(matches, n);
    out.push_str(&format!(
        "{scope}\t{label}\t{n}\t{matches}\t{positive}\t{negative}\t{p}\n"
    ));
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

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn parse_programs(input: &str) -> Result<Vec<String>> {
    let programs: Vec<String> = input
        .split(';')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(ToOwned::to_owned)
        .collect();
    if programs.is_empty() {
        bail!("empty program set in coupling pair");
    }
    Ok(programs)
}

fn parse_sign(input: &str) -> Result<Sign> {
    match input.trim() {
        "positive" | "+" => Ok(Sign::Positive),
        "negative" | "-" => Ok(Sign::Negative),
        other => bail!(
            "invalid hypothesized_sign {:?}; expected positive or negative",
            other
        ),
    }
}

impl Sign {
    fn as_str(self) -> &'static str {
        match self {
            Sign::Positive => "positive",
            Sign::Negative => "negative",
        }
    }
}

fn need_col(headers: &csv::StringRecord, name: &str, path: &Path) -> Result<usize> {
    headers
        .iter()
        .position(|h| h == name)
        .with_context(|| format!("missing column {:?} in {:?}", name, path))
}
