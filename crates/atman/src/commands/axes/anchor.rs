//! `atman axes anchor`: Spearman rho (with subject-level bootstrap CI and the
//! t-approximation p) between subject scores and clinical covariate
//! expressions, per cohort and scope; optional rank-residual partial
//! Spearman given another expression.

use anyhow::{anyhow, bail, Result};
use atman_core::contrast::{percentile, spearman_with_p};
use atman_core::expr::Expr;
use atman_core::stats::{pearson, ranks, spearman};
use atman_core::{derive_sub_seed, SplitMix64};
use clap::Args as ClapArgs;
use serde_json::json;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use super::table::ScoreTable;
use super::{fmt_opt, load_frame, parse_cols, Context};
use crate::io::{atomic_write, hash_labeled_inputs, sidecar_path_for, write_run_sidecar};

#[derive(ClapArgs, Debug)]
pub struct AnchorArgs {
    #[arg(long)]
    scores: PathBuf,
    #[arg(long)]
    cohort_dirs: Option<String>,
    #[arg(long, action = clap::ArgAction::Append)]
    covariates_tsv: Vec<PathBuf>,
    /// Comma-separated cohorts (default: every cohort in the table).
    #[arg(long)]
    cohorts: Option<String>,
    #[arg(long)]
    score_cols: String,
    /// Comma-separated covariate expressions.
    #[arg(long)]
    anchors: String,
    /// Scope (repeatable): `all`, `controls`, `cases`, or `condition=<level>`.
    /// Default `all`.
    #[arg(long, action = clap::ArgAction::Append)]
    scope: Vec<String>,
    /// Expression to partial out (rank-residual partial Spearman rows).
    #[arg(long)]
    partial: Option<String>,
    #[arg(long, default_value_t = 2000)]
    n_bootstrap: usize,
    #[arg(long, default_value_t = 1)]
    seed: u64,
    #[arg(long, default_value_t = 0.95)]
    ci: f64,
    /// Fewer complete pairs than this: report n only.
    #[arg(long, default_value_t = 10)]
    min_n: usize,
    #[arg(long)]
    output: PathBuf,
}

enum Scope {
    All,
    Controls,
    Cases,
    Condition(String),
}

fn parse_scope(text: &str) -> Result<Scope> {
    match text.trim() {
        "all" => Ok(Scope::All),
        "controls" => Ok(Scope::Controls),
        "cases" => Ok(Scope::Cases),
        other => other
            .strip_prefix("condition=")
            .map(|c| Scope::Condition(c.trim().to_string()))
            .ok_or_else(|| {
                anyhow!(
                    "--scope {:?}: expected all, controls, cases, or condition=<level>",
                    other
                )
            }),
    }
}

fn in_scope(scope: &Scope, table: &ScoreTable, i: usize) -> bool {
    match scope {
        Scope::All => true,
        Scope::Controls => table.is_control[i],
        Scope::Cases => !table.is_control[i],
        Scope::Condition(c) => &table.condition[i] == c,
    }
}

fn residual_on(y: &[f64], x: &[f64]) -> Vec<f64> {
    let n = y.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let sxx: f64 = x.iter().map(|v| (v - mx).powi(2)).sum();
    let slope = if sxx > 0.0 {
        x.iter()
            .zip(y)
            .map(|(a, b)| (a - mx) * (b - my))
            .sum::<f64>()
            / sxx
    } else {
        0.0
    };
    let intercept = my - slope * mx;
    y.iter()
        .zip(x)
        .map(|(b, a)| b - (intercept + slope * a))
        .collect()
}

pub fn run(args: AnchorArgs) -> Result<()> {
    let started_at = SystemTime::now();
    if !(0.0..1.0).contains(&args.ci) {
        bail!("--ci must be in (0, 1)");
    }
    let table = ScoreTable::read(&args.scores)?;
    let loaded = load_frame(args.cohort_dirs.as_deref(), &args.covariates_tsv)?;
    let ctx = Context {
        table: &table,
        frame: &loaded.frame,
    };
    let score_cols = parse_cols(&args.score_cols);
    if score_cols.is_empty() {
        bail!("--score-cols is empty");
    }
    for s in &score_cols {
        table.numeric_column(s)?;
    }
    let anchors: Vec<(String, Expr)> = parse_cols(&args.anchors)
        .into_iter()
        .map(|t| {
            Expr::parse(&t)
                .map(|e| (t.clone(), e))
                .map_err(|e| anyhow!("--anchors {:?}: {}", t, e))
        })
        .collect::<Result<_>>()?;
    if anchors.is_empty() {
        bail!("--anchors is empty");
    }
    let partial: Option<(String, Expr)> = match &args.partial {
        Some(t) => Some((
            t.clone(),
            Expr::parse(t).map_err(|e| anyhow!("--partial {:?}: {}", t, e))?,
        )),
        None => None,
    };
    let scope_texts: Vec<String> = if args.scope.is_empty() {
        vec!["all".to_string()]
    } else {
        args.scope.clone()
    };
    let scopes: Vec<(String, Scope)> = scope_texts
        .iter()
        .map(|s| parse_scope(s).map(|sc| (s.trim().to_string(), sc)))
        .collect::<Result<_>>()?;
    let cohorts: Vec<String> = match &args.cohorts {
        Some(c) => parse_cols(c),
        None => table
            .cohort
            .iter()
            .cloned()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect(),
    };
    let alpha = (1.0 - args.ci) / 2.0;

    let mut buf = String::from("cohort\tscope\tgiven\tscore\tanchor\tn\trho\tci_lo\tci_hi\tp\n");
    let mut unit = 0usize;
    for cohort in &cohorts {
        for (scope_name, scope) in &scopes {
            let rows: Vec<usize> = (0..table.len())
                .filter(|&i| &table.cohort[i] == cohort && in_scope(scope, &table, i))
                .collect();
            for score in &score_cols {
                let col = table.numeric_column(score)?;
                for (anchor_text, anchor) in &anchors {
                    let pairs: Vec<(f64, f64)> = rows
                        .iter()
                        .filter_map(|&i| Some((col[i]?, anchor.eval(&|name| ctx.num(i, name))?)))
                        .collect();
                    let n = pairs.len();
                    let mut rng = SplitMix64::new(derive_sub_seed(args.seed, unit));
                    unit += 1;
                    let (rho, p, ci_lo, ci_hi) = if n >= args.min_n && n >= 3 {
                        let xs: Vec<f64> = pairs.iter().map(|p| p.0).collect();
                        let ys: Vec<f64> = pairs.iter().map(|p| p.1).collect();
                        let (rho, p) = spearman_with_p(&xs, &ys)
                            .map(|(r, p)| (Some(r), Some(p)))
                            .unwrap_or((None, None));
                        let mut boots: Vec<f64> = Vec::with_capacity(args.n_bootstrap);
                        for _ in 0..args.n_bootstrap {
                            let mut bx = Vec::with_capacity(n);
                            let mut by = Vec::with_capacity(n);
                            for _ in 0..n {
                                let k = rng.bounded(n);
                                bx.push(xs[k]);
                                by.push(ys[k]);
                            }
                            if let Some(r) = spearman(&bx, &by) {
                                boots.push(r);
                            }
                        }
                        boots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                        let (lo, hi) = if boots.is_empty() {
                            (None, None)
                        } else {
                            (percentile(&boots, alpha), percentile(&boots, 1.0 - alpha))
                        };
                        (rho, p, lo, hi)
                    } else {
                        (None, None, None, None)
                    };
                    buf.push_str(&format!(
                        "{}\t{}\t\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        cohort,
                        scope_name,
                        score,
                        anchor_text,
                        n,
                        fmt_opt(rho),
                        fmt_opt(ci_lo),
                        fmt_opt(ci_hi),
                        fmt_opt(p)
                    ));

                    if let Some((given_text, given)) = &partial {
                        let triples: Vec<(f64, f64, f64)> = rows
                            .iter()
                            .filter_map(|&i| {
                                Some((
                                    col[i]?,
                                    anchor.eval(&|name| ctx.num(i, name))?,
                                    given.eval(&|name| ctx.num(i, name))?,
                                ))
                            })
                            .collect();
                        let n = triples.len();
                        let rho = if n >= args.min_n && n >= 3 {
                            let rx = ranks(&triples.iter().map(|t| t.0).collect::<Vec<_>>());
                            let ry = ranks(&triples.iter().map(|t| t.1).collect::<Vec<_>>());
                            let rz = ranks(&triples.iter().map(|t| t.2).collect::<Vec<_>>());
                            pearson(&residual_on(&rx, &rz), &residual_on(&ry, &rz))
                        } else {
                            None
                        };
                        buf.push_str(&format!(
                            "{}\t{}_partial\t{}\t{}\t{}\t{}\t{}\t\t\t\n",
                            cohort,
                            scope_name,
                            given_text,
                            score,
                            anchor_text,
                            n,
                            fmt_opt(rho)
                        ));
                    }
                }
            }
        }
    }
    atomic_write(&args.output, buf.as_bytes())?;
    eprintln!(
        "axes anchor: cohorts={} scopes={} scores={} anchors={} output={}",
        cohorts.len(),
        scopes.len(),
        score_cols.len(),
        anchors.len(),
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let mut entries: Vec<(String, PathBuf)> = vec![("scores".to_string(), args.scores.clone())];
    entries.extend(loaded.hashed_inputs.iter().cloned());
    let refs: Vec<(&str, &Path)> = entries
        .iter()
        .map(|(l, p)| (l.as_str(), p.as_path()))
        .collect();
    let inputs = hash_labeled_inputs(&refs)?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "axes anchor",
        json!({
            "scores": args.scores.display().to_string(),
            "cohort-dirs": args.cohort_dirs,
            "covariates-tsv": args.covariates_tsv.iter().map(|p| p.display().to_string()).collect::<Vec<_>>(),
            "cohorts": args.cohorts,
            "score-cols": args.score_cols,
            "anchors": args.anchors,
            "scope": scope_texts,
            "partial": args.partial,
            "n-bootstrap": args.n_bootstrap,
            "seed": args.seed,
            "ci": args.ci,
            "min-n": args.min_n,
            "output": args.output.display().to_string(),
        }),
        &inputs,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("axes anchor: sidecar={}", sidecar.display());
    Ok(())
}
