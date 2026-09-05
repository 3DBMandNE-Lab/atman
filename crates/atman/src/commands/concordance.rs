//! `atman concordance`: agreement between per-feature effect tables
//! (DE results, consensus loading vectors) over their shared features —
//! Spearman rho with a feature-resampling bootstrap CI, sign concordance
//! among features significant in both, and Jaccard of the top-N by
//! absolute effect. Tables carry an optional `stage`; a pair present in
//! two stages also gets a delta-rho row with a joint-resample CI.
//! Features are not independent units: the CIs are descriptive.

use anyhow::{anyhow, bail, Context, Result};
use atman_core::contrast::{percentile, spearman_with_p};
use atman_core::stats::{jaccard_top_n, spearman};
use atman_core::{derive_sub_seed, SplitMix64};
use clap::Args as ClapArgs;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{
    atomic_write, format_f64, hash_labeled_inputs, sidecar_path_for, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Manifest TSV: `label, path, effect_col[, feature_col, q_col, stage]`.
    #[arg(long)]
    manifest: PathBuf,
    /// Comma-separated `a:b` label pairs (default: every pair within a stage).
    #[arg(long)]
    pairs: Option<String>,
    /// Top-N by |effect| for the Jaccard column.
    #[arg(long, default_value_t = 100)]
    top_n: usize,
    /// Significance threshold on `q_col` for the hit-based columns.
    #[arg(long, default_value_t = 0.05)]
    q_threshold: f64,
    #[arg(long, default_value_t = 2000)]
    n_bootstrap: usize,
    #[arg(long, default_value_t = 1)]
    seed: u64,
    #[arg(long, default_value_t = 0.95)]
    ci: f64,
    #[arg(long)]
    output: PathBuf,
    /// Delta-rho rows between stages for pairs present in more than one stage.
    #[arg(long)]
    output_delta: Option<PathBuf>,
}

struct TableSpec {
    label: String,
    path: PathBuf,
    effect_col: String,
    feature_col: String,
    q_col: Option<String>,
    stage: String,
}

type EffectTable = BTreeMap<String, (f64, Option<f64>)>;

fn read_manifest(path: &Path) -> Result<Vec<TableSpec>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let col = |n: &str| headers.iter().position(|h| h.trim() == n);
    let need = |n: &str| col(n).ok_or_else(|| anyhow!("missing column {:?} in {:?}", n, path));
    let (c_label, c_path, c_effect) = (need("label")?, need("path")?, need("effect_col")?);
    let (c_feature, c_q, c_stage) = (col("feature_col"), col("q_col"), col("stage"));
    let mut out = Vec::new();
    for row in reader.records() {
        let row = row?;
        let get = |i: usize| row.get(i).unwrap_or("").trim().to_string();
        let label = get(c_label);
        if label.is_empty() {
            continue;
        }
        let stage = c_stage.map(get).unwrap_or_default();
        if out
            .iter()
            .any(|t: &TableSpec| t.label == label && t.stage == stage)
        {
            bail!(
                "duplicate (label, stage) = ({:?}, {:?}) in {:?}",
                label,
                stage,
                path
            );
        }
        out.push(TableSpec {
            label,
            path: PathBuf::from(get(c_path)),
            effect_col: get(c_effect),
            feature_col: c_feature
                .map(get)
                .filter(|s| !s.is_empty())
                .unwrap_or_else(|| "gene_symbol".to_string()),
            q_col: c_q.map(get).filter(|s| !s.is_empty()),
            stage,
        });
    }
    if out.is_empty() {
        bail!("no tables in {:?}", path);
    }
    Ok(out)
}

fn read_effect_table(spec: &TableSpec) -> Result<EffectTable> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .from_path(&spec.path)
        .with_context(|| format!("opening {:?}", spec.path))?;
    let headers = reader.headers()?.clone();
    let need = |n: &str| {
        headers
            .iter()
            .position(|h| h.trim() == n)
            .ok_or_else(|| anyhow!("missing column {:?} in {:?}", n, spec.path))
    };
    let c_feature = need(&spec.feature_col)?;
    let c_effect = need(&spec.effect_col)?;
    let c_q = match &spec.q_col {
        Some(q) => Some(need(q)?),
        None => None,
    };
    let mut out = EffectTable::new();
    for row in reader.records() {
        let row = row?;
        let feature = row.get(c_feature).unwrap_or("").trim().to_string();
        let effect = row
            .get(c_effect)
            .unwrap_or("")
            .trim()
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite());
        let (Some(effect), false) = (effect, feature.is_empty()) else {
            continue;
        };
        let q = c_q
            .and_then(|i| row.get(i))
            .and_then(|v| v.trim().parse::<f64>().ok())
            .filter(|v| v.is_finite());
        if out.insert(feature.clone(), (effect, q)).is_some() {
            bail!("feature {:?} appears twice in {:?}", feature, spec.path);
        }
    }
    if out.is_empty() {
        bail!("no finite effects in {:?}", spec.path);
    }
    Ok(out)
}

struct PairStats {
    n: usize,
    rho: Option<f64>,
    p: Option<f64>,
    n_hits_both: Option<usize>,
    sign_concordance: Option<f64>,
    jaccard: f64,
    x: Vec<f64>,
    y: Vec<f64>,
}

fn pair_stats(
    a: &EffectTable,
    b: &EffectTable,
    features: &[String],
    q_threshold: f64,
    top_n: usize,
) -> PairStats {
    let x: Vec<f64> = features.iter().map(|f| a[f].0).collect();
    let y: Vec<f64> = features.iter().map(|f| b[f].0).collect();
    let (rho, p) = spearman_with_p(&x, &y)
        .map(|(r, p)| (Some(r), Some(p)))
        .unwrap_or((None, None));
    let both_q = !features.is_empty()
        && features.iter().all(|f| a[f].1.is_some())
        && features.iter().all(|f| b[f].1.is_some());
    let (n_hits_both, sign_concordance) = if both_q {
        let hits: Vec<&String> = features
            .iter()
            .filter(|f| a[*f].1.unwrap() < q_threshold && b[*f].1.unwrap() < q_threshold)
            .collect();
        let agree = hits
            .iter()
            .filter(|f| a[**f].0.signum() == b[**f].0.signum())
            .count();
        (
            Some(hits.len()),
            if hits.is_empty() {
                None
            } else {
                Some(agree as f64 / hits.len() as f64)
            },
        )
    } else {
        (None, None)
    };
    let jaccard = jaccard_top_n(&x, &y, top_n.min(features.len()));
    PairStats {
        n: features.len(),
        rho,
        p,
        n_hits_both,
        sign_concordance,
        jaccard,
        x,
        y,
    }
}

fn fmt(v: Option<f64>) -> String {
    v.filter(|x| x.is_finite())
        .map(format_f64)
        .unwrap_or_default()
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    if !(0.0..1.0).contains(&args.ci) {
        bail!("--ci must be in (0, 1)");
    }
    let specs = read_manifest(&args.manifest)?;
    let mut tables: BTreeMap<(String, String), EffectTable> = BTreeMap::new();
    for s in &specs {
        tables.insert((s.stage.clone(), s.label.clone()), read_effect_table(s)?);
    }
    let mut stages: Vec<String> = Vec::new();
    let mut labels: Vec<String> = Vec::new();
    for s in &specs {
        if !stages.contains(&s.stage) {
            stages.push(s.stage.clone());
        }
        if !labels.contains(&s.label) {
            labels.push(s.label.clone());
        }
    }
    let requested: Option<Vec<(String, String)>> = match &args.pairs {
        Some(p) => Some(
            p.split(',')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(|s| {
                    s.split_once(':')
                        .map(|(a, b)| (a.trim().to_string(), b.trim().to_string()))
                        .ok_or_else(|| anyhow!("--pairs entry {:?} must be a:b", s))
                })
                .collect::<Result<_>>()?,
        ),
        None => None,
    };
    let pairs: Vec<(String, String)> = match requested {
        Some(p) => {
            for (a, b) in &p {
                if !labels.contains(a) || !labels.contains(b) {
                    bail!("--pairs names unknown label in {:?}:{:?}", a, b);
                }
            }
            p
        }
        None => {
            let mut v = Vec::new();
            for (i, a) in labels.iter().enumerate() {
                for b in labels.iter().skip(i + 1) {
                    v.push((a.clone(), b.clone()));
                }
            }
            v
        }
    };
    let alpha = (1.0 - args.ci) / 2.0;

    let mut out = String::from(
        "a\tb\tstage\tn_features\trho\tp\tci_lo\tci_hi\tn_hits_both\tsign_concordance\tjaccard_top_n\n",
    );
    let mut delta_out =
        String::from("a\tb\tstage_from\tstage_to\tn_features\tdelta_rho\tdelta_lo\tdelta_hi\n");
    let mut n_rows = 0usize;
    for (pair_idx, (a, b)) in pairs.iter().enumerate() {
        let present: Vec<&String> = stages
            .iter()
            .filter(|s| {
                tables.contains_key(&((*s).clone(), a.clone()))
                    && tables.contains_key(&((*s).clone(), b.clone()))
            })
            .collect();
        if present.is_empty() {
            continue;
        }
        // Common features across every stage where the pair exists.
        let mut common: Option<BTreeSet<String>> = None;
        for stage in &present {
            let ta = &tables[&((*stage).clone(), a.clone())];
            let tb = &tables[&((*stage).clone(), b.clone())];
            let shared: BTreeSet<String> =
                ta.keys().filter(|f| tb.contains_key(*f)).cloned().collect();
            common = Some(match common {
                Some(c) => c.intersection(&shared).cloned().collect(),
                None => shared,
            });
        }
        let features: Vec<String> = common.unwrap_or_default().into_iter().collect();
        if features.len() < 3 {
            bail!("pair {:?}:{:?} shares fewer than 3 features", a, b);
        }
        let mut rng = SplitMix64::new(derive_sub_seed(args.seed, pair_idx));
        let draws: Vec<Vec<usize>> = (0..args.n_bootstrap)
            .map(|_| {
                (0..features.len())
                    .map(|_| rng.bounded(features.len()))
                    .collect()
            })
            .collect();
        let mut boot_rho: BTreeMap<String, Vec<Option<f64>>> = BTreeMap::new();
        let mut rho_by_stage: BTreeMap<String, Option<f64>> = BTreeMap::new();
        for stage in &present {
            let ta = &tables[&((*stage).clone(), a.clone())];
            let tb = &tables[&((*stage).clone(), b.clone())];
            let st = pair_stats(ta, tb, &features, args.q_threshold, args.top_n);
            let boots: Vec<Option<f64>> = draws
                .iter()
                .map(|d| {
                    let bx: Vec<f64> = d.iter().map(|&i| st.x[i]).collect();
                    let by: Vec<f64> = d.iter().map(|&i| st.y[i]).collect();
                    spearman(&bx, &by)
                })
                .collect();
            let mut sorted: Vec<f64> = boots.iter().filter_map(|v| *v).collect();
            sorted.sort_by(|p, q| p.partial_cmp(q).unwrap_or(std::cmp::Ordering::Equal));
            let (lo, hi) = if sorted.is_empty() {
                (None, None)
            } else {
                (percentile(&sorted, alpha), percentile(&sorted, 1.0 - alpha))
            };
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                a,
                b,
                stage,
                st.n,
                fmt(st.rho),
                fmt(st.p),
                fmt(lo),
                fmt(hi),
                st.n_hits_both.map(|n| n.to_string()).unwrap_or_default(),
                fmt(st.sign_concordance),
                format_f64(st.jaccard)
            ));
            n_rows += 1;
            boot_rho.insert((*stage).clone(), boots);
            rho_by_stage.insert((*stage).clone(), st.rho);
        }
        if present.len() >= 2 {
            for (i, from) in present.iter().enumerate() {
                for to in present.iter().skip(i + 1) {
                    let delta = rho_by_stage[*from]
                        .zip(rho_by_stage[*to])
                        .map(|(f, t)| t - f);
                    let mut diffs: Vec<f64> = boot_rho[*from]
                        .iter()
                        .zip(boot_rho[*to].iter())
                        .filter_map(|(f, t)| Some((*t)? - (*f)?))
                        .collect();
                    diffs.sort_by(|p, q| p.partial_cmp(q).unwrap_or(std::cmp::Ordering::Equal));
                    let (lo, hi) = if diffs.is_empty() {
                        (None, None)
                    } else {
                        (percentile(&diffs, alpha), percentile(&diffs, 1.0 - alpha))
                    };
                    delta_out.push_str(&format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        a,
                        b,
                        from,
                        to,
                        features.len(),
                        fmt(delta),
                        fmt(lo),
                        fmt(hi)
                    ));
                }
            }
        }
    }
    if n_rows == 0 {
        bail!("no label pair is present within a common stage");
    }
    atomic_write(&args.output, out.as_bytes())?;
    let mut outputs = vec![args.output.clone()];
    if let Some(p) = &args.output_delta {
        atomic_write(p, delta_out.as_bytes())?;
        outputs.push(p.clone());
    }
    eprintln!(
        "concordance: tables={} pairs={} rows={} output={}",
        specs.len(),
        pairs.len(),
        n_rows,
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let mut entries: Vec<(String, PathBuf)> = vec![("manifest".to_string(), args.manifest.clone())];
    for s in &specs {
        entries.push((format!("table/{}/{}", s.stage, s.label), s.path.clone()));
    }
    let refs: Vec<(&str, &Path)> = entries
        .iter()
        .map(|(l, p)| (l.as_str(), p.as_path()))
        .collect();
    let inputs = hash_labeled_inputs(&refs)?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "concordance",
        json!({
            "manifest": args.manifest.display().to_string(),
            "pairs": args.pairs,
            "top-n": args.top_n,
            "q-threshold": args.q_threshold,
            "n-bootstrap": args.n_bootstrap,
            "seed": args.seed,
            "ci": args.ci,
            "output": args.output.display().to_string(),
            "output-delta": args.output_delta.as_ref().map(|p| p.display().to_string()),
        }),
        &inputs,
        &outputs,
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("concordance: sidecar={}", sidecar.display());
    Ok(())
}
