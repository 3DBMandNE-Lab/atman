//! `atman axes icc`: trait stability of a subject-level score across
//! repeated samples of the same subject — the one-way random-effects
//! intraclass correlation ICC(1) (Shrout & Fleiss 1979), per cohort.
//!
//! With `n` subjects, `k_i` samples for subject `i`, `N = Σ k_i`:
//! `MSB = Σ k_i (ȳ_i − ȳ)² / (n − 1)`, `MSW = Σ_i Σ_j (y_ij − ȳ_i)² / (N − n)`,
//! `k0 = (N − Σ k_i² / N) / (n − 1)` (equals `k` when balanced),
//! `icc1 = (MSB − MSW) / (MSB + (k0 − 1) MSW)`,
//! `between_sd = sqrt(max(0, (MSB − MSW) / k0))`, `within_sd = sqrt(MSW)`.
//! Subjects with a single scored sample are excluded and counted.

use anyhow::{bail, Result};
use clap::Args as ClapArgs;
use serde_json::json;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use super::table::ScoreTable;
use super::{fmt_opt, load_frame, parse_cols, Context};
use crate::io::{atomic_write, hash_labeled_inputs, sidecar_path_for, write_run_sidecar};

#[derive(ClapArgs, Debug)]
pub struct IccArgs {
    #[arg(long)]
    scores: PathBuf,
    #[arg(long)]
    cohort_dirs: Option<String>,
    #[arg(long, action = clap::ArgAction::Append)]
    covariates_tsv: Vec<PathBuf>,
    /// Column holding the subject identifier (score-table text column or a
    /// joined samples.tsv column).
    #[arg(long, default_value = "subject_id")]
    subject_col: String,
    /// Comma-separated cohorts (default: every cohort in the table).
    #[arg(long)]
    cohorts: Option<String>,
    #[arg(long)]
    score_cols: String,
    #[arg(long)]
    output: PathBuf,
}

/// One-way random-effects ICC(1) over `groups` (each a subject's scores).
/// Returns `(icc1, between_sd, within_sd, k0)`; `None` when fewer than two
/// subjects or no within-subject degrees of freedom.
pub fn icc1(groups: &[Vec<f64>]) -> Option<(f64, f64, f64, f64)> {
    let n = groups.len();
    let total: usize = groups.iter().map(Vec::len).sum();
    if n < 2 || total <= n {
        return None;
    }
    let grand = groups.iter().flatten().sum::<f64>() / total as f64;
    let mut ss_between = 0.0;
    let mut ss_within = 0.0;
    let mut sum_k2 = 0.0;
    for g in groups {
        let k = g.len() as f64;
        let mean = g.iter().sum::<f64>() / k;
        ss_between += k * (mean - grand).powi(2);
        ss_within += g.iter().map(|v| (v - mean).powi(2)).sum::<f64>();
        sum_k2 += k * k;
    }
    let msb = ss_between / (n as f64 - 1.0);
    let msw = ss_within / (total - n) as f64;
    let k0 = (total as f64 - sum_k2 / total as f64) / (n as f64 - 1.0);
    let denom = msb + (k0 - 1.0) * msw;
    if !denom.is_finite() || denom <= 0.0 {
        return None;
    }
    let icc = (msb - msw) / denom;
    let between_sd = ((msb - msw) / k0).max(0.0).sqrt();
    Some((icc, between_sd, msw.sqrt(), k0))
}

pub fn run(args: IccArgs) -> Result<()> {
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
    let cohorts: Vec<String> = match &args.cohorts {
        Some(c) => parse_cols(c),
        None => {
            let mut seen = Vec::new();
            for c in &table.cohort {
                if !seen.contains(c) {
                    seen.push(c.clone());
                }
            }
            seen
        }
    };
    let mut buf = String::from(
        "cohort\tscore\tn_subjects\tn_samples\tk_mean\tk0\tn_singletons_excluded\ticc1\tbetween_sd\twithin_sd\n",
    );
    for cohort in &cohorts {
        let rows: Vec<usize> = (0..table.len())
            .filter(|&i| &table.cohort[i] == cohort)
            .collect();
        if rows.is_empty() {
            bail!("no rows for cohort {:?}", cohort);
        }
        for score in &score_cols {
            let col = table.numeric_column(score)?;
            let mut by_subject: BTreeMap<String, Vec<f64>> = BTreeMap::new();
            let mut n_no_subject = 0usize;
            for &i in &rows {
                let Some(v) = col[i] else { continue };
                match ctx.text(i, &args.subject_col) {
                    Some(subject) => by_subject.entry(subject).or_default().push(v),
                    None => n_no_subject += 1,
                }
            }
            if n_no_subject > 0 {
                eprintln!(
                    "axes icc: cohort {} score {}: {} scored samples lack {:?}",
                    cohort, score, n_no_subject, args.subject_col
                );
            }
            let singletons = by_subject.values().filter(|g| g.len() < 2).count();
            let groups: Vec<Vec<f64>> = by_subject.into_values().filter(|g| g.len() >= 2).collect();
            let n_samples: usize = groups.iter().map(Vec::len).sum();
            let k_mean = if groups.is_empty() {
                None
            } else {
                Some(n_samples as f64 / groups.len() as f64)
            };
            let stats = icc1(&groups);
            buf.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                cohort,
                score,
                groups.len(),
                n_samples,
                fmt_opt(k_mean),
                fmt_opt(stats.map(|s| s.3)),
                singletons,
                fmt_opt(stats.map(|s| s.0)),
                fmt_opt(stats.map(|s| s.1)),
                fmt_opt(stats.map(|s| s.2)),
            ));
        }
    }
    atomic_write(&args.output, buf.as_bytes())?;
    eprintln!(
        "axes icc: cohorts={} scores={} output={}",
        cohorts.len(),
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
        "axes icc",
        json!({
            "scores": args.scores.display().to_string(),
            "cohort-dirs": args.cohort_dirs,
            "covariates-tsv": args.covariates_tsv.iter().map(|p| p.display().to_string()).collect::<Vec<_>>(),
            "subject-col": args.subject_col,
            "cohorts": args.cohorts,
            "score-cols": args.score_cols,
            "output": args.output.display().to_string(),
        }),
        &inputs,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("axes icc: sidecar={}", sidecar.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::icc1;

    #[test]
    fn icc1_balanced_matches_hand_computation() {
        // 3 subjects x 2 samples. Subject means 1, 3, 5 (grand 3); within
        // deviations ±0.5 for every subject.
        // SSB = 2*(4+0+4) = 16 ⇒ MSB = 8 ; SSW = 6*0.25 = 1.5 ⇒ MSW = 0.5 ; k0 = 2
        // icc1 = (8 − 0.5) / (8 + 0.5) = 0.8823529411764706
        let groups = vec![vec![0.5, 1.5], vec![2.5, 3.5], vec![4.5, 5.5]];
        let (icc, between, within, k0) = icc1(&groups).unwrap();
        assert!((icc - 7.5 / 8.5).abs() < 1e-12);
        assert!((k0 - 2.0).abs() < 1e-12);
        assert!((within - 0.5f64.sqrt()).abs() < 1e-12);
        assert!((between - (7.5f64 / 2.0).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn icc1_unbalanced_k0() {
        // k = [2, 3]: N = 5, Σk² = 13 ⇒ k0 = (5 − 13/5) / 1 = 2.4
        let groups = vec![vec![1.0, 2.0], vec![4.0, 5.0, 6.0]];
        let (_, _, _, k0) = icc1(&groups).unwrap();
        assert!((k0 - 2.4).abs() < 1e-12);
        assert!(icc1(&[vec![1.0, 2.0]]).is_none());
    }
}
