//! `atman axes build`: pick representative archetype columns as axes,
//! orthogonalize chosen axes against a reference axis within a grouping
//! column (OLS residual with intercept, fitted per group), and z-score every
//! axis globally (ddof = 1).

use anyhow::{anyhow, bail, Context, Result};
use atman_core::contrast::zscore;
use clap::Args as ClapArgs;
use serde_json::json;
use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use super::table::ScoreTable;
use super::{fmt_opt, load_frame};
use crate::design::parse_labeled_paths;
use crate::io::{atomic_write, hash_labeled_inputs, sidecar_path_for, write_run_sidecar};

#[derive(ClapArgs, Debug)]
pub struct BuildArgs {
    /// One wide activations TSV (`sample_id` plus score columns, optionally
    /// `cohort`, `condition`, `is_control`), or a comma-separated
    /// `label=path` list of `align project` outputs whose columns are
    /// prefixed `label_` and outer-joined on `sample_id`.
    #[arg(long)]
    activations: String,

    /// Comma-separated `cohort=dir` canonical directories supplying
    /// `cohort`, `condition`, `is_control` when the activations lack them.
    #[arg(long)]
    cohort_dirs: Option<String>,

    /// Comma-separated `axis=column` list naming each axis's representative
    /// score column, e.g. `axis1=primary_A0002,axis2=axis_residualized_A0001`.
    #[arg(long)]
    representatives: String,

    /// Comma-separated axis names to orthogonalize against `--against`.
    #[arg(long)]
    orthogonalize: Option<String>,

    /// Reference axis name for `--orthogonalize`.
    #[arg(long)]
    against: Option<String>,

    /// Grouping column for the within-group residualization.
    #[arg(long, default_value = "cohort")]
    within: String,

    /// Output TSV path.
    #[arg(long)]
    output: PathBuf,
}

struct Activations {
    table: ScoreTable,
    hashed_inputs: Vec<(String, PathBuf)>,
}

type WideColumns = (
    Vec<String>,
    Vec<String>,
    HashMap<String, Vec<Option<f64>>>,
    HashMap<String, Vec<String>>,
);

/// Read `sample_id` + numeric columns (+ text key columns) from one wide TSV.
fn read_wide(path: &Path) -> Result<WideColumns> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers: Vec<String> = reader
        .headers()?
        .iter()
        .map(|h| h.trim().to_string())
        .collect();
    let c_sample = headers
        .iter()
        .position(|h| h == "sample_id")
        .ok_or_else(|| anyhow!("missing column \"sample_id\" in {:?}", path))?;
    let mut ids = Vec::new();
    let mut raw: Vec<Vec<String>> = vec![Vec::new(); headers.len()];
    for row in reader.records() {
        let row = row.with_context(|| format!("reading {:?}", path))?;
        let sid = row.get(c_sample).unwrap_or("").trim().to_string();
        if sid.is_empty() {
            continue;
        }
        ids.push(sid);
        for (i, bucket) in raw.iter_mut().enumerate() {
            bucket.push(row.get(i).unwrap_or("").trim().to_string());
        }
    }
    let mut order = Vec::new();
    let mut numeric = HashMap::new();
    let mut text = HashMap::new();
    for (i, name) in headers.iter().enumerate() {
        if i == c_sample || name.is_empty() {
            continue;
        }
        let values = &raw[i];
        let key_column = matches!(name.as_str(), "cohort" | "condition" | "is_control");
        let is_numeric = !key_column
            && values
                .iter()
                .filter(|v| !v.is_empty())
                .all(|v| v.parse::<f64>().is_ok());
        if is_numeric {
            order.push(name.clone());
            numeric.insert(
                name.clone(),
                values
                    .iter()
                    .map(|v| {
                        if v.is_empty() {
                            None
                        } else {
                            v.parse::<f64>().ok().filter(|x| x.is_finite())
                        }
                    })
                    .collect(),
            );
        } else {
            text.insert(name.clone(), values.clone());
        }
    }
    Ok((ids, order, numeric, text))
}

fn load_activations(spec: &str, cohort_dirs: Option<&str>) -> Result<Activations> {
    let labeled = spec.contains('=');
    let mut hashed_inputs = Vec::new();
    let mut ids: Vec<String> = Vec::new();
    let mut id_index: HashMap<String, usize> = HashMap::new();
    let mut order: Vec<String> = Vec::new();
    let mut columns: BTreeMap<String, HashMap<String, Option<f64>>> = BTreeMap::new();
    let mut text: HashMap<String, HashMap<String, String>> = HashMap::new();

    let sources: Vec<(String, PathBuf)> = if labeled {
        parse_labeled_paths(spec)?
    } else {
        vec![(String::new(), PathBuf::from(spec.trim()))]
    };
    for (label, path) in &sources {
        let (file_ids, file_order, numeric, file_text) = read_wide(path)?;
        for sid in &file_ids {
            if !id_index.contains_key(sid) {
                id_index.insert(sid.clone(), ids.len());
                ids.push(sid.clone());
            }
        }
        for name in file_order {
            let full = if label.is_empty() {
                name.clone()
            } else {
                format!("{label}_{name}")
            };
            if columns.contains_key(&full) {
                bail!("activation column {:?} appears twice", full);
            }
            order.push(full.clone());
            let mut by_sample = HashMap::new();
            for (k, sid) in file_ids.iter().enumerate() {
                by_sample.insert(sid.clone(), numeric[&name][k]);
            }
            columns.insert(full, by_sample);
        }
        for (name, values) in file_text {
            let entry = text.entry(name).or_default();
            for (k, sid) in file_ids.iter().enumerate() {
                if !values[k].is_empty() {
                    entry.insert(sid.clone(), values[k].clone());
                }
            }
        }
        let key = if label.is_empty() {
            "activations".to_string()
        } else {
            format!("activations/{label}")
        };
        hashed_inputs.push((key, path.clone()));
    }

    let loaded = load_frame(cohort_dirs, &[])?;
    hashed_inputs.extend(loaded.hashed_inputs);

    let mut table = ScoreTable::default();
    for sid in &ids {
        let cohort = text
            .get("cohort")
            .and_then(|m| m.get(sid).cloned())
            .or_else(|| loaded.cohort_of.get(sid).cloned())
            .ok_or_else(|| {
                anyhow!(
                    "sample {:?} has no cohort; add a `cohort` column or --cohort-dirs",
                    sid
                )
            })?;
        let condition = text
            .get("condition")
            .and_then(|m| m.get(sid).cloned())
            .or_else(|| loaded.frame.get(sid, "condition").map(str::to_string))
            .ok_or_else(|| {
                anyhow!(
                    "sample {:?} has no condition; add a `condition` column or --cohort-dirs",
                    sid
                )
            })?;
        let control_text = text
            .get("is_control")
            .and_then(|m| m.get(sid).cloned())
            .or_else(|| loaded.frame.get(sid, "is_control").map(str::to_string))
            .ok_or_else(|| {
                anyhow!(
                    "sample {:?} has no is_control; add an `is_control` column or --cohort-dirs",
                    sid
                )
            })?;
        let is_control = match control_text.to_ascii_lowercase().as_str() {
            "1" | "true" => true,
            "0" | "false" => false,
            other => bail!(
                "is_control value {:?} for sample {:?} is not 0/1",
                other,
                sid
            ),
        };
        table.sample_id.push(sid.clone());
        table.cohort.push(cohort);
        table.condition.push(condition);
        table.is_control.push(is_control);
    }
    for name in &order {
        let by_sample = &columns[name];
        table.numeric.insert(
            name.clone(),
            ids.iter()
                .map(|sid| by_sample.get(sid).copied().flatten())
                .collect(),
        );
    }
    table.numeric_order = order;
    Ok(Activations {
        table,
        hashed_inputs,
    })
}

fn parse_pairs(spec: &str, what: &str) -> Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    for chunk in spec.split(',') {
        let chunk = chunk.trim();
        if chunk.is_empty() {
            continue;
        }
        let (k, v) = chunk
            .split_once('=')
            .ok_or_else(|| anyhow!("{what}: expected name=column, got {:?}", chunk))?;
        out.push((k.trim().to_string(), v.trim().to_string()));
    }
    if out.is_empty() {
        bail!("{what}: no entries");
    }
    Ok(out)
}

/// Residual of `y` on `[1, x]` fitted separately within each group.
fn within_group_residual(
    y: &[Option<f64>],
    x: &[Option<f64>],
    group: &[String],
) -> Result<Vec<Option<f64>>> {
    let mut out = vec![None; y.len()];
    let mut members: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (i, g) in group.iter().enumerate() {
        members.entry(g.as_str()).or_default().push(i);
    }
    for (g, idx) in members {
        let pairs: Vec<(usize, f64, f64)> = idx
            .iter()
            .filter_map(|&i| match (x[i], y[i]) {
                (Some(xv), Some(yv)) => Some((i, xv, yv)),
                _ => None,
            })
            .collect();
        if pairs.len() < 2 {
            continue;
        }
        let n = pairs.len() as f64;
        let mx = pairs.iter().map(|p| p.1).sum::<f64>() / n;
        let my = pairs.iter().map(|p| p.2).sum::<f64>() / n;
        let sxx: f64 = pairs.iter().map(|p| (p.1 - mx).powi(2)).sum();
        if sxx == 0.0 {
            bail!(
                "reference axis is constant within group {:?}; cannot orthogonalize",
                g
            );
        }
        let sxy: f64 = pairs.iter().map(|p| (p.1 - mx) * (p.2 - my)).sum();
        let slope = sxy / sxx;
        let intercept = my - slope * mx;
        for (i, xv, yv) in pairs {
            out[i] = Some(yv - (intercept + slope * xv));
        }
    }
    Ok(out)
}

pub fn run(args: BuildArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let representatives = parse_pairs(&args.representatives, "--representatives")?;
    let orthogonalize: Vec<String> = args
        .orthogonalize
        .as_deref()
        .map(|s| {
            s.split(',')
                .map(|x| x.trim().to_string())
                .filter(|x| !x.is_empty())
                .collect()
        })
        .unwrap_or_default();
    if !orthogonalize.is_empty() && args.against.is_none() {
        bail!("--orthogonalize requires --against <axis>");
    }
    for axis in &orthogonalize {
        if !representatives.iter().any(|(a, _)| a == axis) {
            bail!("--orthogonalize names unknown axis {:?}", axis);
        }
    }
    if let Some(against) = &args.against {
        if !representatives.iter().any(|(a, _)| a == against) {
            bail!("--against names unknown axis {:?}", against);
        }
    }

    let activations = load_activations(&args.activations, args.cohort_dirs.as_deref())?;
    let table = &activations.table;
    for (axis, column) in &representatives {
        if !table.numeric.contains_key(column) {
            bail!("representative for {axis:?}: column {column:?} not found among numeric activation columns");
        }
    }
    let group: Vec<String> = (0..table.len())
        .map(|i| {
            table.text_value(i, &args.within).ok_or_else(|| {
                anyhow!(
                    "--within column {:?} missing for sample {:?}",
                    args.within,
                    table.sample_id[i]
                )
            })
        })
        .collect::<Result<_>>()?;

    let reference: Option<Vec<Option<f64>>> = args.against.as_ref().map(|a| {
        let col = &representatives
            .iter()
            .find(|(axis, _)| axis == a)
            .expect("validated")
            .1;
        table.numeric[col].clone()
    });
    let mut out_cols: Vec<(String, Vec<Option<f64>>)> = Vec::new();
    for (axis, column) in &representatives {
        let rep = table.numeric[column].clone();
        if orthogonalize.contains(axis) {
            let resid = within_group_residual(
                &rep,
                reference.as_ref().expect("against validated"),
                &group,
            )?;
            out_cols.push((format!("{axis}_unorth_raw"), rep.clone()));
            out_cols.push((format!("{axis}_unorth_z"), zscore(&rep)));
            out_cols.push((format!("{axis}_z"), zscore(&resid)));
            let z_pos = out_cols.len() - 1;
            out_cols.insert(z_pos, (format!("{axis}_raw"), resid));
        } else {
            out_cols.push((format!("{axis}_raw"), rep.clone()));
            out_cols.push((format!("{axis}_z"), zscore(&rep)));
        }
    }

    let mut buf = String::from("sample_id\tcohort\tcondition\tis_control");
    for name in &table.numeric_order {
        buf.push('\t');
        buf.push_str(name);
    }
    for (name, _) in &out_cols {
        buf.push('\t');
        buf.push_str(name);
    }
    buf.push('\n');
    for i in 0..table.len() {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}",
            table.sample_id[i], table.cohort[i], table.condition[i], table.is_control[i] as u8
        ));
        for name in &table.numeric_order {
            buf.push('\t');
            buf.push_str(&fmt_opt(table.numeric[name][i]));
        }
        for (_, values) in &out_cols {
            buf.push('\t');
            buf.push_str(&fmt_opt(values[i]));
        }
        buf.push('\n');
    }
    atomic_write(&args.output, buf.as_bytes())?;
    eprintln!(
        "axes build: samples={} axes={} orthogonalized={} output={}",
        table.len(),
        representatives.len(),
        orthogonalize.len(),
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let entries: Vec<(&str, &Path)> = activations
        .hashed_inputs
        .iter()
        .map(|(l, p)| (l.as_str(), p.as_path()))
        .collect();
    let inputs = hash_labeled_inputs(&entries)?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "axes build",
        json!({
            "activations": args.activations,
            "cohort-dirs": args.cohort_dirs,
            "representatives": args.representatives,
            "orthogonalize": args.orthogonalize,
            "against": args.against,
            "within": args.within,
            "output": args.output.display().to_string(),
        }),
        &inputs,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("axes build: sidecar={}", sidecar.display());
    Ok(())
}
