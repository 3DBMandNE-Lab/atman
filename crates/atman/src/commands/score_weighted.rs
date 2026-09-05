//! `atman score weighted`: score every subject of a cohort with a signed
//! per-protein weight vector learned elsewhere (e.g. a DE table's Cohen d).
//! Each protein is z-scored within the scored cohort (ddof = 1, over subjects
//! with a value); `score = Σ w·z / Σ|w|` over the proteins the subject has.
//! With `--groups A-B` a summary row tests the score between the two
//! conditions (pooled-SD Cohen d with CI, Welch p, Mann–Whitney AUC).

use anyhow::{anyhow, bail, Context, Result};
use atman_core::contrast::{auc_mann_whitney, cohen_d_ci, cohen_d_pooled, zscore};
use atman_core::de::{welch_t, PairedTResult};
use clap::Args as ClapArgs;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::design::CollapseGenes;
use crate::io::{
    atomic_write, format_f64, hash_canonical_inputs, hash_labeled_inputs, read_measurements_long,
    read_samples, sidecar_path_for, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct WeightedArgs {
    /// Canonical directory of the cohort to score (or a `residuals --output-canonical-dir`).
    #[arg(long)]
    input_dir: PathBuf,
    /// Weight table (a `de_results.tsv`, a Python DE table, or any TSV).
    #[arg(long)]
    weights: PathBuf,
    /// Feature column in `--weights`.
    #[arg(long, default_value = "gene_symbol")]
    feature_col: String,
    /// Signed weight column in `--weights` (e.g. `effect_size`, `cohen_d`, `mean_diff`).
    #[arg(long)]
    weight_col: String,
    /// Signature label written to the output (default: the weights file stem).
    #[arg(long)]
    signature: Option<String>,
    /// Drop proteins missing in more than this fraction of the cohort's
    /// measured subjects (1.0 = keep all).
    #[arg(long, default_value_t = 1.0)]
    max_missing_fraction: f64,
    /// Minimum proteins with a value for a subject's score to be reported.
    #[arg(long, default_value_t = 10)]
    min_shared: usize,
    /// Optional `A-B` comparison of the score between two condition labels.
    #[arg(long)]
    groups: Option<String>,
    #[arg(long)]
    output: PathBuf,
    /// Summary TSV for `--groups`.
    #[arg(long)]
    output_summary: Option<PathBuf>,
    /// How protein groups sharing a gene symbol are reduced to one value per
    /// sample: `none` (lexically first assay), `mean`, `max-observed`.
    #[arg(long, value_enum, default_value_t = CollapseGenes::None)]
    collapse_genes: CollapseGenes,
}

fn read_weights(path: &Path, feature_col: &str, weight_col: &str) -> Result<BTreeMap<String, f64>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let need = |n: &str| {
        headers
            .iter()
            .position(|h| h.trim() == n)
            .ok_or_else(|| anyhow!("missing column {:?} in {:?}", n, path))
    };
    let (cf, cw) = (need(feature_col)?, need(weight_col)?);
    let mut out = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let f = row.get(cf).unwrap_or("").trim().to_string();
        let w = row
            .get(cw)
            .unwrap_or("")
            .trim()
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite());
        if let (false, Some(w)) = (f.is_empty(), w) {
            if out.insert(f.clone(), w).is_some() {
                bail!("feature {:?} appears twice in {:?}", f, path);
            }
        }
    }
    if out.is_empty() {
        bail!("no finite weights in {:?}", path);
    }
    Ok(out)
}

pub fn run(args: WeightedArgs) -> Result<()> {
    let started_at = SystemTime::now();
    if !(0.0..=1.0).contains(&args.max_missing_fraction) {
        bail!("--max-missing-fraction must be in [0, 1]");
    }
    let signature = args.signature.clone().unwrap_or_else(|| {
        args.weights
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("signature")
            .to_string()
    });
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    let sample_index: HashMap<&str, usize> = samples
        .iter()
        .enumerate()
        .map(|(i, s)| (s.sample_id.as_str(), i))
        .collect();
    let measurements = read_measurements_long(&args.input_dir.join("measurements.tsv"))?;
    let weights = read_weights(&args.weights, &args.feature_col, &args.weight_col)?;

    // gene → sample index → assay → (sum, count); duplicate (sample, assay)
    // cells are averaged, assays sharing a symbol collapse per --collapse-genes.
    type GeneCells = BTreeMap<String, BTreeMap<usize, BTreeMap<String, (f64, usize)>>>;
    let mut acc: GeneCells = BTreeMap::new();
    let mut measured: BTreeSet<usize> = BTreeSet::new();
    for m in &measurements {
        let Some(v) = m.effective_abundance() else {
            continue;
        };
        let Some(&si) = sample_index.get(m.sample_id.as_str()) else {
            continue;
        };
        let gene = m
            .gene_symbol
            .clone()
            .unwrap_or_else(|| m.assay_id.0.clone());
        let e = acc
            .entry(gene)
            .or_default()
            .entry(si)
            .or_default()
            .entry(m.assay_id.0.clone())
            .or_insert((0.0, 0));
        e.0 += v;
        e.1 += 1;
        measured.insert(si);
    }
    let n_measured = measured.len().max(1);
    let collapse = args.collapse_genes;
    let mut n_genes_multi_assay = 0usize;
    let mut shared: Vec<(String, f64, Vec<Option<f64>>)> = Vec::new();
    for (gene, by_sample) in &acc {
        let Some(&w) = weights.get(gene) else {
            continue;
        };
        let mut assay_counts: BTreeMap<String, usize> = BTreeMap::new();
        for assays in by_sample.values() {
            for assay in assays.keys() {
                *assay_counts.entry(assay.clone()).or_default() += 1;
            }
        }
        if assay_counts.len() > 1 {
            n_genes_multi_assay += 1;
        }
        let representative = collapse
            .representative(&assay_counts)
            .cloned()
            .unwrap_or_default();
        let raw: Vec<Option<f64>> = (0..samples.len())
            .map(|i| {
                by_sample.get(&i).and_then(|assays| {
                    let values: Vec<(String, f64)> = assays
                        .iter()
                        .map(|(a, (s, n))| (a.clone(), s / *n as f64))
                        .collect();
                    collapse.collapse(&values, &representative)
                })
            })
            .collect();
        let n_present = raw.iter().filter(|v| v.is_some()).count();
        let missing = 1.0 - n_present as f64 / n_measured as f64;
        if missing > args.max_missing_fraction + 1e-12 {
            continue;
        }
        shared.push((gene.clone(), w, zscore(&raw)));
    }
    if n_genes_multi_assay > 0 {
        eprintln!(
            "score weighted: {} gene symbols are carried by more than one assay; --collapse-genes {}",
            n_genes_multi_assay,
            collapse.as_str()
        );
    }
    if shared.is_empty() {
        bail!(
            "no proteins shared between {:?} and the cohort after filtering",
            args.weights
        );
    }

    let mut scores: Vec<Option<f64>> = Vec::with_capacity(samples.len());
    let mut n_used: Vec<usize> = Vec::with_capacity(samples.len());
    for i in 0..samples.len() {
        let mut num = 0.0;
        let mut den = 0.0;
        let mut used = 0usize;
        for (_, w, z) in &shared {
            if let Some(zv) = z[i] {
                num += w * zv;
                den += w.abs();
                used += 1;
            }
        }
        n_used.push(used);
        scores.push(if used >= args.min_shared && den > 0.0 {
            Some(num / den)
        } else {
            None
        });
    }

    let mut buf = String::from(
        "sample_id\tsubject_id\tcondition\tis_control\tsignature\tn_shared\tn_used\tscore\n",
    );
    for (i, s) in samples.iter().enumerate() {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            s.sample_id,
            s.subject_id.clone().unwrap_or_default(),
            s.condition.clone().unwrap_or_default(),
            s.is_control as u8,
            signature,
            shared.len(),
            n_used[i],
            scores[i].map(format_f64).unwrap_or_default()
        ));
    }
    atomic_write(&args.output, buf.as_bytes())?;
    let mut outputs = vec![args.output.clone()];

    if let Some(groups) = &args.groups {
        let path = args
            .output_summary
            .as_ref()
            .ok_or_else(|| anyhow!("--groups requires --output-summary"))?;
        let (a, b) = groups
            .split_once('-')
            .ok_or_else(|| anyhow!("--groups must be A-B"))?;
        let pick = |cond: &str| -> Vec<f64> {
            samples
                .iter()
                .enumerate()
                .filter(|(_, s)| s.condition.as_deref() == Some(cond))
                .filter_map(|(i, _)| scores[i])
                .collect()
        };
        let (va, vb) = (pick(a.trim()), pick(b.trim()));
        if va.is_empty() || vb.is_empty() {
            bail!("--groups {}: a group has no scored subjects", groups);
        }
        let d = cohen_d_pooled(&va, &vb);
        let ci = d.and_then(|d| cohen_d_ci(d, va.len(), vb.len()));
        let p = match welch_t(&va, &vb, 2) {
            PairedTResult::Computed { p_value, .. } => Some(p_value),
            PairedTResult::Skipped { .. } => None,
        };
        let f = |v: Option<f64>| v.map(format_f64).unwrap_or_default();
        let text = format!(
            "signature\tcomparison\tn_shared_proteins\tn_case\tn_control\tcohen_d\td_ci_lo\td_ci_hi\twelch_p\tauc\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            signature,
            groups,
            shared.len(),
            va.len(),
            vb.len(),
            f(d),
            f(ci.map(|c| c.0)),
            f(ci.map(|c| c.1)),
            f(p),
            f(auc_mann_whitney(&va, &vb))
        );
        atomic_write(path, text.as_bytes())?;
        outputs.push(path.clone());
    }
    eprintln!(
        "score weighted: signature={} samples={} shared_proteins={} output={}",
        signature,
        samples.len(),
        shared.len(),
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let mut inputs = hash_canonical_inputs(
        &args.input_dir,
        &["measurements.tsv", "samples.tsv", "proteins.tsv"],
    )?;
    inputs.extend(hash_labeled_inputs(&[("weights", args.weights.as_path())])?);
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "score weighted",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "weights": args.weights.display().to_string(),
            "feature-col": args.feature_col,
            "weight-col": args.weight_col,
            "signature": signature,
            "max-missing-fraction": args.max_missing_fraction,
            "min-shared": args.min_shared,
            "groups": args.groups,
            "output": args.output.display().to_string(),
            "output-summary": args.output_summary.as_ref().map(|p| p.display().to_string()),
            "collapse-genes": collapse.as_str(),
        }),
        &inputs,
        &outputs,
        started_at,
        finished_at,
        Some({
            let mut extras = serde_json::Map::new();
            extras.insert(
                "gene_symbol_collapse".into(),
                json!({
                    "rule": collapse.as_str(),
                    "n_genes_with_multiple_assays": n_genes_multi_assay,
                }),
            );
            extras
        }),
    )?;
    eprintln!("score weighted: sidecar={}", sidecar.display());
    Ok(())
}
