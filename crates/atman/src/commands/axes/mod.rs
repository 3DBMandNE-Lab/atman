//! `atman axes`: subject-level score tables — building axis scores from
//! archetype activations and the score-level statistics the CSF
//! CrossDisease manuscript reports on them. "Axes" means any set of
//! per-subject numeric score columns.

pub mod manifest;
pub mod table;

mod anchor;
mod build;
mod contrast;
mod displacement;
mod groups;
mod icc;
mod loco;

use anyhow::{bail, Result};
use atman_core::SplitMix64;
use clap::{Args as ClapArgs, Subcommand};
use std::collections::HashMap;
use std::path::PathBuf;

use crate::design::{parse_labeled_paths, CovariateFrame};
use crate::io::format_f64;
use table::ScoreTable;

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Select representative archetype columns as axes, orthogonalize
    /// within cohort, and z-score globally.
    Build(build::BuildArgs),
    /// Disease modulation of subject scores: Cohen d, Welch, AUC, adjusted
    /// OLS with nested designs, BH within family, subject bootstrap.
    Contrast(contrast::ContrastArgs),
    /// Per-group score summaries, Kruskal–Wallis, reference-vs-group effects.
    Groups(groups::GroupsArgs),
    /// Spearman anchoring of scores to clinical covariates with bootstrap CI.
    Anchor(anchor::AnchorArgs),
    /// Displacement vectors (case minus control centroid) and pairwise
    /// cosines with bootstrap CI.
    Displacement(displacement::DisplacementArgs),
    /// Leave-one-cohort-out stability of group centroids on global z-scores.
    Loco(loco::LocoArgs),
    /// One-way random-effects ICC(1) of a score across repeated samples per subject.
    Icc(icc::IccArgs),
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Build(a) => build::run(a),
        Command::Contrast(a) => contrast::run(a),
        Command::Groups(a) => groups::run(a),
        Command::Anchor(a) => anchor::run(a),
        Command::Displacement(a) => displacement::run(a),
        Command::Loco(a) => loco::run(a),
        Command::Icc(a) => icc::run(a),
    }
}

/// Row-wise metadata view over a score table plus joined covariates.
pub struct Context<'a> {
    pub table: &'a ScoreTable,
    pub frame: &'a CovariateFrame,
}

impl Context<'_> {
    pub fn text(&self, i: usize, col: &str) -> Option<String> {
        self.table.text_value(i, col).or_else(|| {
            self.frame
                .get(&self.table.sample_id[i], col)
                .map(str::to_string)
        })
    }

    pub fn num(&self, i: usize, col: &str) -> Option<f64> {
        self.text(i, col)
            .and_then(|s| s.trim().parse::<f64>().ok())
            .filter(|v| v.is_finite())
    }
}

pub struct LoadedFrame {
    pub frame: CovariateFrame,
    /// sample_id → cohort label, from `--cohort-dirs`.
    pub cohort_of: HashMap<String, String>,
    /// (label, path) pairs to hash into the sidecar.
    pub hashed_inputs: Vec<(String, PathBuf)>,
}

/// Join every `--cohort-dirs` `samples.tsv` and every `--covariates-tsv`.
pub fn load_frame(cohort_dirs: Option<&str>, covariates_tsv: &[PathBuf]) -> Result<LoadedFrame> {
    let mut frame = CovariateFrame::new();
    let mut cohort_of = HashMap::new();
    let mut hashed_inputs = Vec::new();
    if let Some(spec) = cohort_dirs {
        for (label, dir) in parse_labeled_paths(spec)? {
            let path = dir.join("samples.tsv");
            let mut probe = CovariateFrame::new();
            probe.load_tsv(&path)?;
            for sid in probe.sample_ids() {
                if let Some(prev) = cohort_of.insert(sid.clone(), label.clone()) {
                    if prev != label {
                        bail!(
                            "sample {:?} appears in cohorts {:?} and {:?}",
                            sid,
                            prev,
                            label
                        );
                    }
                }
            }
            frame.load_tsv(&path)?;
            hashed_inputs.push((format!("cohort_dirs/{label}/samples.tsv"), path));
        }
    }
    for path in covariates_tsv {
        frame.load_tsv(path)?;
        let label = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("covariates")
            .to_string();
        hashed_inputs.push((format!("covariates_tsv/{label}"), path.clone()));
    }
    Ok(LoadedFrame {
        frame,
        cohort_of,
        hashed_inputs,
    })
}

/// Comma-split outside parentheses (so `log10(a/b), age` stays two items).
pub fn parse_cols(text: &str) -> Vec<String> {
    atman_core::expr::split_top_level(text, ',')
}

pub fn fmt_opt(v: Option<f64>) -> String {
    match v {
        Some(x) if x.is_finite() => format_f64(x),
        _ => String::new(),
    }
}

/// Resample with replacement inside each group, preserving group sizes;
/// returns the concatenated resampled indices (group order preserved).
pub fn resample_within_groups(rng: &mut SplitMix64, groups: &[&[usize]]) -> Vec<usize> {
    let mut out = Vec::with_capacity(groups.iter().map(|g| g.len()).sum());
    for g in groups {
        for _ in 0..g.len() {
            out.push(g[rng.bounded(g.len())]);
        }
    }
    out
}
