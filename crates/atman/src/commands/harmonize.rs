//! `atman harmonize`: cross-cohort harmonisation with a fit-then-apply
//! contract, for benchmarking harmonisation methods by held-out
//! transfer.
//!
//! `fit` reads the training cohorts and writes a model. `apply` reads
//! that model and one cohort, and nothing else — so a held-out
//! evaluation cannot leak, structurally rather than by discipline. The
//! model records the training cohort labels AND a hash of their inputs,
//! because a label can be reused over different data and a record that
//! omits the field cannot disagree with reality.

use anyhow::{bail, Context, Result};
use atman_core::harmonize::{
    apply_with, fit, CohortData, DirectionCentering, HarmonizeMethod, MethodSpec,
};
use clap::{Args as ClapArgs, Subcommand};
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{
    atomic_write, format_float, hash_labeled_inputs, read_measurements_long, read_samples,
    sidecar_path_for, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Fit a harmonisation model and a disease direction on training
    /// cohorts, and write a transportable model artifact.
    Fit(FitArgs),
    /// Apply a fitted model to one cohort it has never seen, producing
    /// a per-subject transfer score.
    Apply(ApplyArgs),
}

#[derive(ClapArgs, Debug)]
pub struct FitArgs {
    /// Training cohorts as `LABEL=path,LABEL=path,...`. Each path is a
    /// canonical directory with `measurements.tsv` and `samples.tsv`.
    #[arg(long, value_delimiter = ',')]
    cohorts: Vec<String>,

    /// Harmonisation method: `zscore` (the null: per-cohort location and
    /// scale only), `rank` (within-sample, learns nothing), `quantile`
    /// (against a reference profile learned on the training cohorts), or
    /// `reference-protein` (per-sample division by the geometric mean of
    /// a low-variance reference set chosen on training cohorts only).
    #[arg(long)]
    method: String,

    /// Number of reference proteins for `--method reference-protein`.
    #[arg(long, default_value_t = 20)]
    reference_k: usize,

    /// Column of `samples.tsv` naming the case/control grouping.
    #[arg(long, default_value = "condition")]
    condition_col: String,

    /// Condition value treated as the case arm. Every other non-empty
    /// value is a control.
    #[arg(long)]
    case_value: String,

    /// Shuffle case/control labels within each training cohort before
    /// learning the direction: the negative-control arm.
    ///
    /// Every method scores something on a held-out cohort, so a
    /// benchmark whose entries all beat zero measures nothing. A method
    /// whose permuted arm still separates cases is detecting structure
    /// unrelated to disease.
    ///
    /// ONE permuted fit is ONE DRAW from the null, not the null. On a
    /// 60-feature synthetic fixture a single draw produced a held-out
    /// separation of −4.45 against a real effect of +11.90, purely by
    /// chance in the direction that was learned. Run several seeds and
    /// compare the real effect against that distribution; a single
    /// permuted run can mislead in either direction.
    #[arg(long, default_value_t = false)]
    permute_labels: bool,

    /// Seed for the label permutation.
    #[arg(long, default_value_t = 20260906)]
    seed: u64,

    /// Output model path (JSON).
    #[arg(long)]
    output_model: PathBuf,
}

#[derive(ClapArgs, Debug)]
pub struct ApplyArgs {
    /// Model written by `harmonize fit`.
    #[arg(long)]
    model: PathBuf,

    /// Held-out cohort as `LABEL=path`.
    #[arg(long)]
    cohort: String,

    #[arg(long, default_value = "condition")]
    condition_col: String,

    #[arg(long)]
    case_value: String,

    /// Subtract the direction's mean over each subject's own used
    /// features before scoring.
    ///
    /// Strongly recommended for every per-sample normalisation, and
    /// close to mandatory for `reference-protein`. Those methods leave a
    /// per-subject term in the harmonised value — for reference-protein
    /// the score literally carries `−anchor × mean(direction)` — so if
    /// the anchor differs between arms in the held-out cohort, ANY
    /// direction separates them, including one learned from shuffled
    /// labels. Measured on real data: a permuted null running +0.542 to
    /// +0.698 against a real effect of +0.611, i.e. the shuffled arm
    /// separated cases better than the real one.
    ///
    /// Off by default only so that runs made before this existed remain
    /// reproducible. Turn it on unless you have a reason not to.
    #[arg(long, default_value_t = false)]
    center_direction: bool,

    /// Output per-subject transfer scores (TSV).
    #[arg(long)]
    output: PathBuf,
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Fit(a) => run_fit(a),
        Command::Apply(a) => run_apply(a),
    }
}

fn parse_cohort_spec(spec: &str) -> Result<(String, PathBuf)> {
    let (label, path) = spec
        .split_once('=')
        .with_context(|| format!("cohort spec {spec:?} must be LABEL=path"))?;
    if label.trim().is_empty() {
        bail!("cohort spec {spec:?} has an empty label");
    }
    Ok((label.trim().to_string(), PathBuf::from(path.trim())))
}

/// Load one cohort's subject × gene matrix and case/control labels.
fn load_cohort(
    label: &str,
    dir: &Path,
    condition_col: &str,
    case_value: &str,
) -> Result<CohortData> {
    let records = read_measurements_long(&dir.join("measurements.tsv"))?;
    let samples = read_samples(&dir.join("samples.tsv"))?;

    // condition_col is `condition` unless the caller names another
    // column; only `condition` is in the canonical Sample struct, so a
    // different column is refused rather than silently ignored.
    if condition_col != "condition" {
        bail!(
            "--condition-col {condition_col:?}: only `condition` is available from samples.tsv \
             in this command; other columns would be silently ignored"
        );
    }
    let mut label_for: BTreeMap<&str, bool> = BTreeMap::new();
    for s in &samples {
        if let Some(c) = &s.condition {
            if c.trim().is_empty() {
                continue;
            }
            label_for.insert(s.sample_id.as_str(), c.trim() == case_value);
        }
    }
    if label_for.is_empty() {
        bail!("cohort {label:?}: no sample has a non-empty {condition_col}");
    }
    if !label_for.values().any(|v| *v) {
        bail!(
            "cohort {label:?}: no sample has {condition_col} == {case_value:?}, so there is no \
             case arm"
        );
    }
    if !label_for.values().any(|v| !*v) {
        bail!("cohort {label:?}: every sample is a case, so there is no control arm");
    }

    let mut genes: BTreeSet<String> = BTreeSet::new();
    let mut cells: BTreeMap<(String, String), f64> = BTreeMap::new();
    for r in &records {
        if r.dropped_by_qc {
            continue;
        }
        let Some(gene) = &r.gene_symbol else { continue };
        let v = r.abundance.as_f64();
        if !v.is_finite() || !label_for.contains_key(r.sample_id.as_str()) {
            continue;
        }
        genes.insert(gene.clone());
        cells.insert((r.sample_id.clone(), gene.clone()), v);
    }
    let subject_ids: Vec<String> = label_for.keys().map(|s| s.to_string()).collect();
    let features: Vec<String> = genes.into_iter().collect();
    let values: Vec<Vec<f64>> = subject_ids
        .iter()
        .map(|s| {
            features
                .iter()
                .map(|g| {
                    cells
                        .get(&(s.clone(), g.clone()))
                        .copied()
                        .unwrap_or(f64::NAN)
                })
                .collect()
        })
        .collect();
    let is_case: Vec<bool> = subject_ids.iter().map(|s| label_for[s.as_str()]).collect();
    Ok(CohortData {
        label: label.to_string(),
        subject_ids,
        features,
        values,
        is_case,
    })
}

fn run_fit(args: FitArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let spec = MethodSpec::parse(&args.method, args.reference_k).with_context(|| {
        format!(
            "--method {:?}; expected zscore, rank, quantile, or reference-protein",
            args.method
        )
    })?;
    if args.cohorts.len() < 2 {
        bail!("harmonize fit needs at least 2 training cohorts");
    }
    let specs: Vec<(String, PathBuf)> = args
        .cohorts
        .iter()
        .map(|s| parse_cohort_spec(s))
        .collect::<Result<_>>()?;

    let mut cohorts = Vec::with_capacity(specs.len());
    for (label, dir) in &specs {
        cohorts.push(load_cohort(
            label,
            dir,
            &args.condition_col,
            &args.case_value,
        )?);
    }
    let model =
        fit(&cohorts, spec, args.permute_labels, args.seed).map_err(|e| anyhow::anyhow!(e))?;

    // Hash the training inputs into the model. Labels alone would let a
    // replay reuse a label over different data and still "prove" the
    // held-out cohort was absent.
    let refs: Vec<(String, &Path)> = specs
        .iter()
        .map(|(l, p)| (l.clone(), p.join("measurements.tsv")))
        .map(|(l, p)| (l, Box::leak(p.into_boxed_path()) as &Path))
        .collect();
    let refs_borrowed: Vec<(&str, &Path)> = refs.iter().map(|(l, p)| (l.as_str(), *p)).collect();
    let inputs_sha256 = hash_labeled_inputs(&refs_borrowed)?;

    let reference_features = match &model.method {
        HarmonizeMethod::ReferenceProtein { reference_features } => {
            Some(reference_features.clone())
        }
        _ => None,
    };
    let doc = json!({
        "schema_version": 1,
        "method": model.method.as_str(),
        "method_is_fitted": model.method.is_fitted(),
        "permuted_labels": model.permuted,
        "seed": args.seed,
        "case_value": args.case_value,
        "fit_cohorts": model.fit_cohorts,
        "fit_inputs_sha256": inputs_sha256,
        "n_features": model.features.len(),
        "features": model.features,
        "direction": model.direction,
        "reference_features": reference_features,
        "quantile_reference": match &model.method {
            HarmonizeMethod::Quantile { reference } => Some(reference.clone()),
            _ => None,
        },
    });
    if let Some(parent) = args.output_model.parent() {
        std::fs::create_dir_all(parent).ok();
    }
    atomic_write(
        &args.output_model,
        serde_json::to_string_pretty(&doc)?.as_bytes(),
    )?;
    eprintln!(
        "harmonize fit: method={} fitted={} permuted={} cohorts={} shared_features={} model={}",
        model.method.as_str(),
        model.method.is_fitted(),
        model.permuted,
        model.fit_cohorts.len(),
        model.features.len(),
        args.output_model.display(),
    );
    if model.permuted {
        eprintln!(
            "harmonize fit: NOTE this is the permuted-label negative control arm. Its scores are \
             not a disease result. This is ONE DRAW from the null, not the null: run several \
             --seed values and compare the real effect against that distribution, because a \
             single permuted draw can produce a large separation in either direction by chance."
        );
    }

    let finished_at = SystemTime::now();
    let sidecar = sidecar_path_for(&args.output_model);
    write_run_sidecar(
        &sidecar,
        "harmonize fit",
        json!({
            "cohorts": args.cohorts,
            "method": args.method,
            "reference-k": args.reference_k,
            "case-value": args.case_value,
            "permute-labels": args.permute_labels,
            "seed": args.seed,
            "output-model": args.output_model.display().to_string(),
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output_model),
        started_at,
        finished_at,
        None,
    )?;
    Ok(())
}

fn run_apply(args: ApplyArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let text = std::fs::read_to_string(&args.model)
        .with_context(|| format!("reading model {:?}", args.model))?;
    let doc: serde_json::Value =
        serde_json::from_str(&text).with_context(|| format!("parsing model {:?}", args.model))?;

    let method_name = doc["method"].as_str().unwrap_or_default().to_string();
    let features: Vec<String> = doc["features"]
        .as_array()
        .map(|a| {
            a.iter()
                .map(|v| v.as_str().unwrap_or_default().to_string())
                .collect()
        })
        .unwrap_or_default();
    let direction: Vec<f64> = doc["direction"]
        .as_array()
        .map(|a| a.iter().map(|v| v.as_f64().unwrap_or(0.0)).collect())
        .unwrap_or_default();
    let fit_cohorts: Vec<String> = doc["fit_cohorts"]
        .as_array()
        .map(|a| {
            a.iter()
                .map(|v| v.as_str().unwrap_or_default().to_string())
                .collect()
        })
        .unwrap_or_default();
    if features.is_empty() || direction.len() != features.len() {
        bail!(
            "model {:?} is malformed: features/direction mismatch",
            args.model
        );
    }
    let method = match method_name.as_str() {
        "zscore" => HarmonizeMethod::ZScore,
        "rank" => HarmonizeMethod::Rank,
        "quantile" => HarmonizeMethod::Quantile {
            reference: doc["quantile_reference"]
                .as_array()
                .map(|a| a.iter().map(|v| v.as_f64().unwrap_or(0.0)).collect())
                .unwrap_or_default(),
        },
        "reference-protein" => HarmonizeMethod::ReferenceProtein {
            reference_features: doc["reference_features"]
                .as_array()
                .map(|a| {
                    a.iter()
                        .map(|v| v.as_str().unwrap_or_default().to_string())
                        .collect()
                })
                .unwrap_or_default(),
        },
        other => bail!("model {:?} names unknown method {other:?}", args.model),
    };
    let model = atman_core::harmonize::HarmonizeModel {
        method,
        features,
        direction,
        fit_cohorts,
        permuted: doc["permuted_labels"].as_bool().unwrap_or(false),
    };

    let (label, dir) = parse_cohort_spec(&args.cohort)?;
    let cohort = load_cohort(&label, &dir, &args.condition_col, &args.case_value)?;
    let centering = if args.center_direction {
        DirectionCentering::PerSubject
    } else {
        DirectionCentering::None
    };
    let scores = apply_with(&model, &cohort, centering).map_err(|e| anyhow::anyhow!(e))?;
    if !args.center_direction && model.method.as_str() != "zscore" {
        eprintln!(
            "harmonize apply: warning: --center-direction is off and {} is a per-sample \
             normalisation, which leaves a per-subject term the direction's mean projects onto. \
             If that term differs between arms in this cohort, a direction learned from shuffled \
             labels will separate cases too. Run the permuted arm and compare, or pass \
             --center-direction.",
            model.method.as_str()
        );
    }

    let mut buf = String::from("subject_id\ttransfer_score\tis_case\tn_features_used\n");
    for s in &scores {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            s.subject_id,
            format_float(s.score),
            u8::from(s.is_case),
            s.n_features_used,
        ));
    }
    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent).ok();
    }
    atomic_write(&args.output, buf.as_bytes())?;

    let n_cases = scores.iter().filter(|s| s.is_case).count();
    eprintln!(
        "harmonize apply: method={} permuted={} cohort={} subjects={} cases={} output={}",
        model.method.as_str(),
        model.permuted,
        label,
        scores.len(),
        n_cases,
        args.output.display(),
    );

    let finished_at = SystemTime::now();
    let model_path: &Path = &args.model;
    let cohort_path = dir.join("measurements.tsv");
    let inputs_sha256 =
        hash_labeled_inputs(&[("model", model_path), ("cohort", cohort_path.as_path())])?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "harmonize apply",
        json!({
            "model": args.model.display().to_string(),
            "cohort": args.cohort,
            "case-value": args.case_value,
            "output": args.output.display().to_string(),
            // Carried forward so a reader of the scores does not have to
            // open the model to learn they are a negative control.
            "permuted_labels": model.permuted,
            "center-direction": args.center_direction,
            "fit_cohorts": model.fit_cohorts,
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    Ok(())
}
