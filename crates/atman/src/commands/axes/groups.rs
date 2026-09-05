//! `atman axes groups`: per-group summaries of subject scores within one
//! cohort (mean with t-based CI, medians of covariate expressions),
//! Kruskal–Wallis across groups, and pooled-SD Cohen d plus Welch p of a
//! reference group against every other group.

use anyhow::{bail, Result};
use atman_core::contrast::{cohen_d_pooled, kruskal_wallis, median, t_ci_mean};
use atman_core::de::{welch_t, PairedTResult};
use atman_core::expr::Expr;
use atman_core::stats::mean;
use clap::Args as ClapArgs;
use serde_json::json;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use super::table::ScoreTable;
use super::{fmt_opt, load_frame, parse_cols, Context};
use crate::io::{atomic_write, hash_labeled_inputs, sidecar_path_for, write_run_sidecar};

#[derive(ClapArgs, Debug)]
pub struct GroupsArgs {
    #[arg(long)]
    scores: PathBuf,
    #[arg(long)]
    cohort_dirs: Option<String>,
    #[arg(long, action = clap::ArgAction::Append)]
    covariates_tsv: Vec<PathBuf>,
    /// Cohort whose rows are summarized.
    #[arg(long)]
    cohort: String,
    /// Column whose values define the groups.
    #[arg(long)]
    group_by: String,
    /// Comma-separated group order (default: all observed values, sorted).
    #[arg(long)]
    groups: Option<String>,
    /// Comma-separated score columns.
    #[arg(long)]
    score_cols: String,
    /// Comma-separated covariate expressions summarized by their median.
    #[arg(long)]
    median_cols: Option<String>,
    /// Reference group for pairwise Cohen d / Welch p rows.
    #[arg(long)]
    reference: Option<String>,
    #[arg(long, default_value_t = 0.95)]
    ci: f64,
    #[arg(long)]
    output: PathBuf,
    /// Kruskal–Wallis and reference-vs-group rows.
    #[arg(long)]
    output_tests: Option<PathBuf>,
}

pub fn run(args: GroupsArgs) -> Result<()> {
    let started_at = SystemTime::now();
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
    let median_exprs: Vec<(String, Expr)> = args
        .median_cols
        .as_deref()
        .map(parse_cols)
        .unwrap_or_default()
        .into_iter()
        .map(|t| {
            Expr::parse(&t)
                .map(|e| (t.clone(), e))
                .map_err(|e| anyhow::anyhow!("--median-cols {:?}: {}", t, e))
        })
        .collect::<Result<_>>()?;

    let cohort_rows: Vec<usize> = (0..table.len())
        .filter(|&i| table.cohort[i] == args.cohort)
        .collect();
    if cohort_rows.is_empty() {
        bail!("no rows for cohort {:?}", args.cohort);
    }
    let observed: BTreeSet<String> = cohort_rows
        .iter()
        .filter_map(|&i| ctx.text(i, &args.group_by))
        .collect();
    let groups: Vec<String> = match &args.groups {
        Some(g) => {
            let listed = parse_cols(g);
            for name in &listed {
                if !observed.contains(name) {
                    bail!(
                        "group {:?} not observed in column {:?} of cohort {:?}",
                        name,
                        args.group_by,
                        args.cohort
                    );
                }
            }
            listed
        }
        None => observed.iter().cloned().collect(),
    };
    if groups.is_empty() {
        bail!("no groups found in column {:?}", args.group_by);
    }
    if let Some(r) = &args.reference {
        if !groups.contains(r) {
            bail!("--reference {:?} is not one of the groups", r);
        }
    }
    let members: Vec<Vec<usize>> = groups
        .iter()
        .map(|g| {
            cohort_rows
                .iter()
                .copied()
                .filter(|&i| ctx.text(i, &args.group_by).as_deref() == Some(g.as_str()))
                .collect()
        })
        .collect();

    let mut header = vec!["group".to_string(), "n".to_string()];
    for (t, _) in &median_exprs {
        header.push(format!("{t}_median"));
    }
    for s in &score_cols {
        header.push(format!("{s}_mean"));
        header.push(format!("{s}_ci_lo"));
        header.push(format!("{s}_ci_hi"));
    }
    let mut buf = header.join("\t");
    buf.push('\n');
    for (g, idx) in groups.iter().zip(&members) {
        let mut cells = vec![g.clone(), idx.len().to_string()];
        for (_, e) in &median_exprs {
            let vals: Vec<f64> = idx
                .iter()
                .filter_map(|&i| e.eval(&|name| ctx.num(i, name)))
                .collect();
            cells.push(fmt_opt(median(&vals)));
        }
        for s in &score_cols {
            let col = table.numeric_column(s)?;
            let vals: Vec<f64> = idx.iter().filter_map(|&i| col[i]).collect();
            let m = if vals.is_empty() {
                None
            } else {
                Some(mean(&vals))
            };
            let ci = if vals.len() >= 3 {
                t_ci_mean(&vals, args.ci)
            } else {
                None
            };
            cells.push(fmt_opt(m));
            cells.push(fmt_opt(ci.map(|c| c.0)));
            cells.push(fmt_opt(ci.map(|c| c.1)));
        }
        buf.push_str(&cells.join("\t"));
        buf.push('\n');
    }
    atomic_write(&args.output, buf.as_bytes())?;
    let mut outputs = vec![args.output.clone()];

    if let Some(path) = &args.output_tests {
        let mut tb = String::from("kind\tscore\treference\tgroup\tn_ref\tn_group\tstatistic\tp\n");
        for s in &score_cols {
            let col = table.numeric_column(s)?;
            let per_group: Vec<Vec<f64>> = members
                .iter()
                .map(|idx| idx.iter().filter_map(|&i| col[i]).collect())
                .collect();
            let kw = kruskal_wallis(&per_group);
            tb.push_str(&format!(
                "kruskal\t{}\t\t\t\t\t{}\t{}\n",
                s,
                fmt_opt(kw.map(|k| k.h)),
                fmt_opt(kw.map(|k| k.p_value))
            ));
            if let Some(r) = &args.reference {
                let ri = groups.iter().position(|g| g == r).expect("validated");
                for (gi, g) in groups.iter().enumerate() {
                    if gi == ri {
                        continue;
                    }
                    let a = &per_group[ri];
                    let b = &per_group[gi];
                    let p = match welch_t(a, b, 2) {
                        PairedTResult::Computed { p_value, .. } => Some(p_value),
                        PairedTResult::Skipped { .. } => None,
                    };
                    tb.push_str(&format!(
                        "reference_vs_group\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        s,
                        r,
                        g,
                        a.len(),
                        b.len(),
                        fmt_opt(cohen_d_pooled(a, b)),
                        fmt_opt(p)
                    ));
                }
            }
        }
        atomic_write(path, tb.as_bytes())?;
        outputs.push(path.clone());
    }
    eprintln!(
        "axes groups: cohort={} groups={} scores={} output={}",
        args.cohort,
        groups.len(),
        score_cols.len(),
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
        "axes groups",
        json!({
            "scores": args.scores.display().to_string(),
            "cohort-dirs": args.cohort_dirs,
            "covariates-tsv": args.covariates_tsv.iter().map(|p| p.display().to_string()).collect::<Vec<_>>(),
            "cohort": args.cohort,
            "group-by": args.group_by,
            "groups": args.groups,
            "score-cols": args.score_cols,
            "median-cols": args.median_cols,
            "reference": args.reference,
            "ci": args.ci,
            "output": args.output.display().to_string(),
            "output-tests": args.output_tests.as_ref().map(|p| p.display().to_string()),
        }),
        &inputs,
        &outputs,
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("axes groups: sidecar={}", sidecar.display());
    Ok(())
}
