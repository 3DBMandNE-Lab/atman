//! `atman axes displacement`: per-contrast displacement vectors (case
//! centroid minus control centroid, column-wise over subjects with a value)
//! in one or more score spaces, pairwise cosine with a subject-level
//! bootstrap CI, and vector norms.

use anyhow::{anyhow, bail, Result};
use atman_core::contrast::{percentile, zscore};
use atman_core::stats::cosine;
use atman_core::{derive_sub_seed, SplitMix64};
use clap::Args as ClapArgs;
use serde_json::json;
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use super::manifest::{read_contrast_manifest, select_rows, Selection};
use super::table::ScoreTable;
use super::{fmt_opt, load_frame, parse_cols, resample_within_groups, Context};
use crate::io::{
    atomic_write, format_f64, hash_labeled_inputs, sidecar_path_for, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct DisplacementArgs {
    #[arg(long)]
    scores: PathBuf,
    #[arg(long)]
    cohort_dirs: Option<String>,
    #[arg(long, action = clap::ArgAction::Append)]
    covariates_tsv: Vec<PathBuf>,
    /// Contrast manifest (same format as `axes contrast`).
    #[arg(long)]
    manifest: PathBuf,
    /// Semicolon-separated `name=col,col,...` score spaces.
    #[arg(long)]
    score_sets: String,
    /// `global` (z-score each column over every subject; default) or `none`.
    #[arg(long, default_value = "global")]
    standardize: String,
    #[arg(long, default_value_t = 2000)]
    n_bootstrap: usize,
    #[arg(long, default_value_t = 1)]
    seed: u64,
    #[arg(long, default_value_t = 0.95)]
    ci: f64,
    #[arg(long)]
    output: PathBuf,
    #[arg(long)]
    output_vectors: Option<PathBuf>,
}

/// One score space: `(name, column names, per-column values)`.
type Space = (String, Vec<String>, Vec<Vec<Option<f64>>>);

fn centroid(cols: &[Vec<Option<f64>>], rows: &[usize]) -> Vec<f64> {
    cols.iter()
        .map(|c| {
            let vals: Vec<f64> = rows.iter().filter_map(|&i| c[i]).collect();
            if vals.is_empty() {
                f64::NAN
            } else {
                vals.iter().sum::<f64>() / vals.len() as f64
            }
        })
        .collect()
}

fn displacement(cols: &[Vec<Option<f64>>], sel_case: &[usize], sel_control: &[usize]) -> Vec<f64> {
    centroid(cols, sel_case)
        .iter()
        .zip(centroid(cols, sel_control))
        .map(|(a, b)| a - b)
        .collect()
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

pub fn run(args: DisplacementArgs) -> Result<()> {
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
    let contrasts = read_contrast_manifest(&args.manifest)?;
    let standardize = match args.standardize.as_str() {
        "global" => true,
        "none" => false,
        other => bail!("--standardize {:?}; expected global or none", other),
    };
    let mut spaces: Vec<Space> = Vec::new();
    for chunk in args
        .score_sets
        .split(';')
        .map(str::trim)
        .filter(|s| !s.is_empty())
    {
        let (name, cols) = chunk
            .split_once('=')
            .ok_or_else(|| anyhow!("--score-sets entry {:?} must be name=col,...", chunk))?;
        let cols = parse_cols(cols);
        if cols.is_empty() {
            bail!("--score-sets {:?}: no columns", name);
        }
        let data: Vec<Vec<Option<f64>>> = cols
            .iter()
            .map(|c| {
                table
                    .numeric_column(c)
                    .map(|v| if standardize { zscore(v) } else { v.to_vec() })
            })
            .collect::<Result<_>>()?;
        spaces.push((name.trim().to_string(), cols, data));
    }
    if spaces.is_empty() {
        bail!("--score-sets is empty");
    }
    let selections: Vec<Selection> = contrasts
        .iter()
        .map(|c| select_rows(c, &ctx))
        .collect::<Result<_>>()?;
    let alpha = (1.0 - args.ci) / 2.0;

    let mut vectors = String::from("label\tspace\tscore\tdisplacement\tn_case\tn_control\n");
    for (c, sel) in contrasts.iter().zip(&selections) {
        for (space, cols, data) in &spaces {
            let d = displacement(data, &sel.case, &sel.control);
            for (col, v) in cols.iter().zip(d) {
                vectors.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\n",
                    c.label,
                    space,
                    col,
                    fmt_opt(Some(v)),
                    sel.case.len(),
                    sel.control.len()
                ));
            }
        }
    }

    let mut out = String::from("a\tb\tspace\tcosine\tci_lo\tci_hi\tnorm_a\tnorm_b\n");
    let mut pair_idx = 0usize;
    for i in 0..contrasts.len() {
        for j in (i + 1)..contrasts.len() {
            let mut rng = SplitMix64::new(derive_sub_seed(args.seed, pair_idx));
            pair_idx += 1;
            let (sa, sb) = (&selections[i], &selections[j]);
            for (space, _, data) in &spaces {
                let da = displacement(data, &sa.case, &sa.control);
                let db = displacement(data, &sb.case, &sb.control);
                let cos = cosine(&da, &db);
                let mut boots: Vec<f64> = Vec::with_capacity(args.n_bootstrap);
                for _ in 0..args.n_bootstrap {
                    let ra = resample_within_groups(&mut rng, &[&sa.case, &sa.control]);
                    let rb = resample_within_groups(&mut rng, &[&sb.case, &sb.control]);
                    let (ca, cca) = ra.split_at(sa.case.len());
                    let (cb, ccb) = rb.split_at(sb.case.len());
                    if let Some(v) =
                        cosine(&displacement(data, ca, cca), &displacement(data, cb, ccb))
                    {
                        boots.push(v);
                    }
                }
                boots.sort_by(|p, q| p.partial_cmp(q).unwrap_or(std::cmp::Ordering::Equal));
                let (lo, hi) = if boots.is_empty() {
                    (None, None)
                } else {
                    (percentile(&boots, alpha), percentile(&boots, 1.0 - alpha))
                };
                out.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    contrasts[i].label,
                    contrasts[j].label,
                    space,
                    fmt_opt(cos),
                    fmt_opt(lo),
                    fmt_opt(hi),
                    format_f64(norm(&da)),
                    format_f64(norm(&db))
                ));
            }
        }
    }
    atomic_write(&args.output, out.as_bytes())?;
    let mut outputs = vec![args.output.clone()];
    if let Some(p) = &args.output_vectors {
        atomic_write(p, vectors.as_bytes())?;
        outputs.push(p.clone());
    }
    eprintln!(
        "axes displacement: contrasts={} spaces={} pairs={} output={}",
        contrasts.len(),
        spaces.len(),
        pair_idx,
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let mut entries: Vec<(String, PathBuf)> = vec![
        ("scores".to_string(), args.scores.clone()),
        ("manifest".to_string(), args.manifest.clone()),
    ];
    entries.extend(loaded.hashed_inputs.iter().cloned());
    let refs: Vec<(&str, &Path)> = entries
        .iter()
        .map(|(l, p)| (l.as_str(), p.as_path()))
        .collect();
    let inputs = hash_labeled_inputs(&refs)?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "axes displacement",
        json!({
            "scores": args.scores.display().to_string(),
            "cohort-dirs": args.cohort_dirs,
            "covariates-tsv": args.covariates_tsv.iter().map(|p| p.display().to_string()).collect::<Vec<_>>(),
            "manifest": args.manifest.display().to_string(),
            "score-sets": args.score_sets,
            "standardize": args.standardize,
            "n-bootstrap": args.n_bootstrap,
            "seed": args.seed,
            "ci": args.ci,
            "output": args.output.display().to_string(),
            "output-vectors": args.output_vectors.as_ref().map(|p| p.display().to_string()),
        }),
        &inputs,
        &outputs,
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("axes displacement: sidecar={}", sidecar.display());
    Ok(())
}
