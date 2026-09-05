//! `atman axes tree`: average-linkage tree of contrasts (or groups) in a
//! score space with subject-bootstrap clade support.
//!
//! Leaves are the manifest contrasts. `--leaf displacement` (default) uses
//! case centroid − control centroid with distance `1 − cosine`;
//! `--leaf centroid` uses the case-group centroid in z units with Euclidean
//! distance. Bootstrap replicates resample subjects within each contrast's
//! groups (one SplitMix64 stream seeded by `--seed`), rebuild the tree, and
//! support is the fraction of replicates containing each reference clade.

use anyhow::{bail, Result};
use atman_core::contrast::zscore;
use atman_core::stats::cosine;
use atman_core::SplitMix64;
use clap::Args as ClapArgs;
use serde_json::json;
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use super::manifest::{read_contrast_manifest, select_rows, Selection};
use super::table::ScoreTable;
use super::{load_frame, parse_cols, resample_within_groups, Context};
use crate::io::{hash_labeled_inputs, sidecar_path_for, write_run_sidecar};
use crate::tree::{clade_support, write_linkage, write_support, Tree};

#[derive(ClapArgs, Debug)]
pub struct TreeArgs {
    #[arg(long)]
    scores: PathBuf,
    #[arg(long)]
    cohort_dirs: Option<String>,
    #[arg(long, action = clap::ArgAction::Append)]
    covariates_tsv: Vec<PathBuf>,
    /// Contrast manifest (same format as `axes contrast`); each row is a leaf.
    #[arg(long)]
    manifest: PathBuf,
    /// Comma-separated score columns spanning the space.
    #[arg(long)]
    score_cols: String,
    /// `displacement` (case − control, cosine distance) or `centroid`
    /// (case-group centroid, Euclidean distance).
    #[arg(long, default_value = "displacement")]
    leaf: String,
    /// `cosine` or `euclidean`; defaults to the natural choice for `--leaf`.
    #[arg(long)]
    distance: Option<String>,
    /// `global` (z-score each column over every subject; default) or `none`.
    #[arg(long, default_value = "global")]
    standardize: String,
    #[arg(long, default_value_t = 1000)]
    n_bootstrap: usize,
    #[arg(long, default_value_t = 1)]
    seed: u64,
    /// Quote Newick labels instead of sanitizing them.
    #[arg(long, default_value_t = false)]
    quote_labels: bool,
    /// Linkage table (`left, right, distance, n`).
    #[arg(long)]
    output_linkage: PathBuf,
    /// Clade support table.
    #[arg(long)]
    output_support: Option<PathBuf>,
    /// Newick tree with support labels.
    #[arg(long)]
    output_newick: Option<PathBuf>,
}

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

fn leaf_vector(
    cols: &[Vec<Option<f64>>],
    case: &[usize],
    control: &[usize],
    displacement: bool,
) -> Vec<f64> {
    let c = centroid(cols, case);
    if displacement {
        c.iter()
            .zip(centroid(cols, control))
            .map(|(a, b)| a - b)
            .collect()
    } else {
        c
    }
}

fn dissimilarity(vectors: &[Vec<f64>], cosine_distance: bool) -> Vec<Vec<f64>> {
    let n = vectors.len();
    let mut d = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let v = if cosine_distance {
                1.0 - cosine(&vectors[i], &vectors[j]).unwrap_or(0.0)
            } else {
                vectors[i]
                    .iter()
                    .zip(&vectors[j])
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt()
            };
            d[i][j] = v;
            d[j][i] = v;
        }
    }
    d
}

pub fn run(args: TreeArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let displacement = match args.leaf.as_str() {
        "displacement" => true,
        "centroid" => false,
        other => bail!("--leaf {:?}; expected displacement or centroid", other),
    };
    let cosine_distance = match args.distance.as_deref() {
        Some("cosine") => true,
        Some("euclidean") => false,
        None => displacement,
        Some(other) => bail!("--distance {:?}; expected cosine or euclidean", other),
    };
    let standardize = match args.standardize.as_str() {
        "global" => true,
        "none" => false,
        other => bail!("--standardize {:?}; expected global or none", other),
    };
    let table = ScoreTable::read(&args.scores)?;
    let loaded = load_frame(args.cohort_dirs.as_deref(), &args.covariates_tsv)?;
    let ctx = Context {
        table: &table,
        frame: &loaded.frame,
    };
    let contrasts = read_contrast_manifest(&args.manifest)?;
    if contrasts.len() < 2 {
        bail!("a tree needs at least two leaves (manifest rows)");
    }
    let cols: Vec<Vec<Option<f64>>> = parse_cols(&args.score_cols)
        .iter()
        .map(|c| {
            table
                .numeric_column(c)
                .map(|v| if standardize { zscore(v) } else { v.to_vec() })
        })
        .collect::<Result<_>>()?;
    if cols.is_empty() {
        bail!("--score-cols is empty");
    }
    let selections: Vec<Selection> = contrasts
        .iter()
        .map(|c| select_rows(c, &ctx))
        .collect::<Result<_>>()?;
    let labels: Vec<String> = contrasts.iter().map(|c| c.label.clone()).collect();

    let vectors: Vec<Vec<f64>> = selections
        .iter()
        .map(|s| leaf_vector(&cols, &s.case, &s.control, displacement))
        .collect();
    let reference =
        Tree::from_dissimilarity(labels.clone(), &dissimilarity(&vectors, cosine_distance));

    let mut rng = SplitMix64::new(args.seed);
    let mut replicates: Vec<Tree> = Vec::with_capacity(args.n_bootstrap);
    for _ in 0..args.n_bootstrap {
        let vs: Vec<Vec<f64>> = selections
            .iter()
            .map(|s| {
                let draw = resample_within_groups(&mut rng, &[&s.case, &s.control]);
                let (case, control) = draw.split_at(s.case.len());
                leaf_vector(&cols, case, control, displacement)
            })
            .collect();
        replicates.push(Tree::from_dissimilarity(
            labels.clone(),
            &dissimilarity(&vs, cosine_distance),
        ));
    }
    let support = clade_support(&reference, &replicates);

    write_linkage(&args.output_linkage, &reference)?;
    let mut outputs = vec![args.output_linkage.clone()];
    if let Some(p) = &args.output_support {
        write_support(p, &reference, &support, args.n_bootstrap)?;
        outputs.push(p.clone());
    }
    if let Some(p) = &args.output_newick {
        let labelled: Vec<Option<f64>> = support
            .iter()
            .map(|s| if s.is_finite() { Some(*s) } else { None })
            .collect();
        let mut text = reference.newick(&labelled, args.quote_labels);
        text.push('\n');
        crate::io::atomic_write(p, text.as_bytes())?;
        outputs.push(p.clone());
    }
    eprintln!(
        "axes tree: leaves={} leaf={} distance={} n_bootstrap={} output={}",
        labels.len(),
        args.leaf,
        if cosine_distance {
            "cosine"
        } else {
            "euclidean"
        },
        args.n_bootstrap,
        args.output_linkage.display()
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
    let sidecar = sidecar_path_for(&args.output_linkage);
    write_run_sidecar(
        &sidecar,
        "axes tree",
        json!({
            "scores": args.scores.display().to_string(),
            "cohort-dirs": args.cohort_dirs,
            "covariates-tsv": args.covariates_tsv.iter().map(|p| p.display().to_string()).collect::<Vec<_>>(),
            "manifest": args.manifest.display().to_string(),
            "score-cols": args.score_cols,
            "leaf": args.leaf,
            "distance": if cosine_distance { "cosine" } else { "euclidean" },
            "standardize": args.standardize,
            "n-bootstrap": args.n_bootstrap,
            "seed": args.seed,
            "quote-labels": args.quote_labels,
            "output-linkage": args.output_linkage.display().to_string(),
            "output-support": args.output_support.as_ref().map(|p| p.display().to_string()),
            "output-newick": args.output_newick.as_ref().map(|p| p.display().to_string()),
        }),
        &inputs,
        &outputs,
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("axes tree: sidecar={}", sidecar.display());
    Ok(())
}
