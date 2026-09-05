//! `atman axes loco`: leave-one-cohort-out stability of group centroids on
//! globally z-scored scores. For each score column and condition group,
//! drop each cohort in turn, recompute the global z on the remaining
//! subjects, and report the group's mean z before and after.

use anyhow::{anyhow, bail, Result};
use atman_core::contrast::zscore;
use clap::Args as ClapArgs;
use serde_json::json;
use std::path::PathBuf;
use std::time::SystemTime;

use super::table::ScoreTable;
use super::{fmt_opt, parse_cols};
use crate::io::{atomic_write, hash_labeled_inputs, sidecar_path_for, write_run_sidecar};

#[derive(ClapArgs, Debug)]
pub struct LocoArgs {
    #[arg(long)]
    scores: PathBuf,
    /// Comma-separated raw score columns (z is recomputed here).
    #[arg(long)]
    score_cols: String,
    /// Semicolon-separated `Group=cond|cond,...` condition groups.
    #[arg(long)]
    groups: String,
    #[arg(long)]
    output: PathBuf,
}

fn mean_over(values: &[Option<f64>], rows: &[usize]) -> Option<f64> {
    let v: Vec<f64> = rows.iter().filter_map(|&i| values[i]).collect();
    if v.is_empty() {
        None
    } else {
        Some(v.iter().sum::<f64>() / v.len() as f64)
    }
}

pub fn run(args: LocoArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let table = ScoreTable::read(&args.scores)?;
    let score_cols = parse_cols(&args.score_cols);
    if score_cols.is_empty() {
        bail!("--score-cols is empty");
    }
    for s in &score_cols {
        table.numeric_column(s)?;
    }
    let groups: Vec<(String, Vec<String>)> = args
        .groups
        .split(';')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|chunk| {
            let (name, conds) = chunk
                .split_once('=')
                .ok_or_else(|| anyhow!("--groups entry {:?} must be Group=cond|cond", chunk))?;
            Ok((
                name.trim().to_string(),
                conds
                    .split('|')
                    .map(|c| c.trim().to_string())
                    .filter(|c| !c.is_empty())
                    .collect(),
            ))
        })
        .collect::<Result<_>>()?;
    if groups.is_empty() {
        bail!("--groups is empty");
    }
    let mut cohorts: Vec<String> = Vec::new();
    for c in &table.cohort {
        if !cohorts.contains(c) {
            cohorts.push(c.clone());
        }
    }
    let n = table.len();
    let mut buf =
        String::from("score\tgroup\tdropped_cohort\tn_remaining\tfull_z_mean\tloo_z_mean\tdelta\n");
    for score in &score_cols {
        let raw = table.numeric_column(score)?;
        let full_z = zscore(raw);
        for (group, conds) in &groups {
            let members: Vec<usize> = (0..n)
                .filter(|&i| conds.contains(&table.condition[i]))
                .collect();
            let full_mean = mean_over(&full_z, &members);
            for dropped in &cohorts {
                let remaining: Vec<usize> =
                    (0..n).filter(|&i| &table.cohort[i] != dropped).collect();
                let sub_raw: Vec<Option<f64>> = remaining.iter().map(|&i| raw[i]).collect();
                let sub_z = zscore(&sub_raw);
                let group_pos: Vec<usize> = remaining
                    .iter()
                    .enumerate()
                    .filter(|(_, &i)| conds.contains(&table.condition[i]))
                    .map(|(k, _)| k)
                    .collect();
                if group_pos.is_empty() {
                    continue;
                }
                let loo_mean = mean_over(&sub_z, &group_pos);
                buf.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    score,
                    group,
                    dropped,
                    group_pos.len(),
                    fmt_opt(full_mean),
                    fmt_opt(loo_mean),
                    fmt_opt(loo_mean.zip(full_mean).map(|(l, f)| l - f))
                ));
            }
        }
    }
    atomic_write(&args.output, buf.as_bytes())?;
    eprintln!(
        "axes loco: scores={} groups={} cohorts={} output={}",
        score_cols.len(),
        groups.len(),
        cohorts.len(),
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let inputs = hash_labeled_inputs(&[("scores", args.scores.as_path())])?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "axes loco",
        json!({
            "scores": args.scores.display().to_string(),
            "score-cols": args.score_cols,
            "groups": args.groups,
            "output": args.output.display().to_string(),
        }),
        &inputs,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("axes loco: sidecar={}", sidecar.display());
    Ok(())
}
