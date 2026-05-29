use anyhow::{bail, Context, Result};
use atman_core::variance_decomposition::{decompose_archetype_variance, FixedFactor, VarianceRow};
use clap::Args as ClapArgs;
use csv::ReaderBuilder;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{atomic_write, read_samples, sidecar_path_for, write_run_sidecar};

#[derive(ClapArgs, Debug)]
pub struct VarianceArgs {
    /// Per-sample activation TSV produced by `atman decompose ica`:
    /// columns `cohort, subject_id, sample_id, program, activation`.
    #[arg(long)]
    activations: PathBuf,

    /// samples.tsv providing per-sample metadata columns used as
    /// fixed factors or the random-intercept grouping.
    #[arg(long)]
    samples: PathBuf,

    /// Mixed-model formula: `<fixed_terms> [+ (1|group_col)]`.
    /// Fixed terms are comma- or plus-separated column names from
    /// `samples.tsv` (or the built-ins `cohort`, `condition`).
    /// Categorical columns are one-hot encoded (alphabetically
    /// first level dropped as reference); numeric columns are single
    /// coefficient columns. At most one random-intercept term is
    /// supported; the grouping column must exist in samples.tsv.
    /// Example: `--factors "cohort + condition + (1|subject_id)"`.
    #[arg(long)]
    factors: String,

    /// Minimum samples per archetype for the fit to run. Below this
    /// threshold, the archetype emits a skip row with NaN stats.
    #[arg(long, default_value_t = 5)]
    min_samples: usize,

    /// Output TSV path.
    #[arg(long)]
    output: PathBuf,
}

pub(super) fn run_variance(args: VarianceArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let samples = read_samples(&args.samples)?;
    if samples.is_empty() {
        bail!("samples.tsv at {:?} is empty", args.samples);
    }
    let (fixed_terms, random_group) = parse_variance_formula(&args.factors)?;
    if fixed_terms.is_empty() {
        bail!("--factors must include at least one fixed term");
    }

    // Read per-sample extras from the raw samples.tsv file (headers
    // include user-defined covariates beyond the canonical cols).
    let extras = read_samples_extras(&args.samples)?;

    // Activations: one row per (sample_id, program). Re-shape to
    // per-program per-sample.
    let ActivationsTable {
        by_key: program_by_sample,
        programs_order,
        samples_order,
    } = read_activations(&args.activations)?;
    if programs_order.is_empty() {
        bail!("no programs found in {:?}", args.activations);
    }
    if samples_order.is_empty() {
        bail!("no samples found in {:?}", args.activations);
    }

    // Build the design matrix and fixed-factor column ranges. Only
    // samples that appear in both `activations` and `samples.tsv`
    // are kept, and only samples with finite covariates.
    let VarianceDesign {
        design,
        fixed_factors,
        kept_sample_ids,
        group_labels,
    } = build_variance_design(
        &samples,
        &extras,
        &fixed_terms,
        random_group.as_deref(),
        &samples_order,
    )?;
    if design.is_empty() {
        bail!("no samples remain after fixed-covariate complete-case filtering");
    }

    // Build per-program activation rows in the reduced sample order.
    let mut activations_mat: Vec<Vec<f64>> = Vec::with_capacity(programs_order.len());
    for program in &programs_order {
        let mut row = Vec::with_capacity(kept_sample_ids.len());
        for sid in &kept_sample_ids {
            let v = program_by_sample
                .get(&(program.clone(), sid.clone()))
                .copied()
                .unwrap_or(f64::NAN);
            row.push(v);
        }
        activations_mat.push(row);
    }

    let rows: Vec<VarianceRow> = decompose_archetype_variance(
        &programs_order,
        &activations_mat,
        &design,
        &fixed_factors,
        group_labels.as_deref(),
        args.min_samples,
    )
    .map_err(|e| anyhow::anyhow!(e))?;

    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output dir {:?}", parent))?;
    }
    write_variance_tsv(&args.output, &rows, &fixed_factors, random_group.as_deref())?;
    eprintln!(
        "decompose variance: wrote {} archetypes to {}",
        rows.len(),
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = crate::io::hash_labeled_inputs(&[
        ("activations", args.activations.as_path()),
        ("samples", args.samples.as_path()),
    ])?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "decompose variance",
        json!({
            "activations": args.activations.display().to_string(),
            "samples": args.samples.display().to_string(),
            "factors": args.factors,
            "min-samples": args.min_samples,
            "output": args.output.display().to_string(),
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("decompose variance: sidecar={}", sidecar.display());
    Ok(())
}

fn parse_variance_formula(formula: &str) -> Result<(Vec<String>, Option<String>)> {
    let mut fixed = Vec::new();
    let mut random: Option<String> = None;
    for raw in formula.split(['+', ',']) {
        let term = raw.trim();
        if term.is_empty() {
            continue;
        }
        if let Some(rest) = term.strip_prefix('(').and_then(|s| s.strip_suffix(')')) {
            // (1|group_col)
            let inner = rest.trim();
            if let Some(group) = inner.strip_prefix("1|") {
                let g = group.trim();
                if g.is_empty() {
                    bail!("empty random-intercept grouping column in formula");
                }
                if random.is_some() {
                    bail!("at most one random-intercept term is supported");
                }
                random = Some(g.to_string());
            } else {
                bail!(
                    "unsupported random-effect term {:?}; expected `(1|<column>)`",
                    term
                );
            }
        } else {
            fixed.push(term.to_string());
        }
    }
    Ok((fixed, random))
}

/// Load the user-visible samples.tsv extra columns (beyond the
/// canonical `sample_id, subject_id, condition, is_control,
/// sample_type, ingest_order` that `read_samples` returns). Needed
/// for variance-decomposition formulas that reference cohort / sex
/// / batch / etc.
fn read_samples_extras(path: &Path) -> Result<BTreeMap<String, BTreeMap<String, String>>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers: Vec<String> = reader
        .headers()
        .with_context(|| format!("reading headers from {:?}", path))?
        .iter()
        .map(|s| s.to_string())
        .collect();
    let sample_col = headers
        .iter()
        .position(|h| h == "sample_id")
        .ok_or_else(|| anyhow::anyhow!("samples.tsv {:?} missing sample_id", path))?;
    let mut out: BTreeMap<String, BTreeMap<String, String>> = BTreeMap::new();
    for row in reader.records() {
        let row = row.with_context(|| format!("reading {:?}", path))?;
        let sample = row.get(sample_col).unwrap_or_default().to_string();
        if sample.is_empty() {
            continue;
        }
        let mut extras: BTreeMap<String, String> = BTreeMap::new();
        for (i, h) in headers.iter().enumerate() {
            if i == sample_col {
                continue;
            }
            let v = row.get(i).unwrap_or_default().to_string();
            extras.insert(h.clone(), v);
        }
        out.insert(sample, extras);
    }
    Ok(out)
}

struct ActivationsTable {
    /// `(program, sample) → activation value`.
    by_key: BTreeMap<(String, String), f64>,
    /// Stable alphabetical program ordering.
    programs_order: Vec<String>,
    /// Sample insertion order as first seen in the TSV.
    samples_order: Vec<String>,
}

fn read_activations(path: &Path) -> Result<ActivationsTable> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers: Vec<String> = reader
        .headers()
        .with_context(|| format!("reading headers from {:?}", path))?
        .iter()
        .map(|s| s.to_string())
        .collect();
    let col = |name: &str| -> Result<usize> {
        headers
            .iter()
            .position(|h| h == name)
            .ok_or_else(|| anyhow::anyhow!("activations.tsv {:?} missing column {:?}", path, name))
    };
    let c_sample = col("sample_id")?;
    let c_program = col("program")?;
    let c_activation = col("activation")?;
    let mut by_key: BTreeMap<(String, String), f64> = BTreeMap::new();
    let mut programs: BTreeSet<String> = BTreeSet::new();
    let mut samples_seen: BTreeSet<String> = BTreeSet::new();
    let mut samples_order: Vec<String> = Vec::new();
    for row in reader.records() {
        let row = row.with_context(|| format!("reading {:?}", path))?;
        let sample = row.get(c_sample).unwrap_or_default().to_string();
        let program = row.get(c_program).unwrap_or_default().to_string();
        let value: f64 = row
            .get(c_activation)
            .unwrap_or_default()
            .parse()
            .unwrap_or(f64::NAN);
        if sample.is_empty() || program.is_empty() {
            continue;
        }
        if samples_seen.insert(sample.clone()) {
            samples_order.push(sample.clone());
        }
        programs.insert(program.clone());
        by_key.insert((program, sample), value);
    }
    let programs_order: Vec<String> = programs.into_iter().collect();
    Ok(ActivationsTable {
        by_key,
        programs_order,
        samples_order,
    })
}

/// Build the design matrix + per-factor column ranges + per-sample
/// random-group labels in the order of `samples_order`. Drops samples
/// whose fixed covariates can't be resolved.
struct VarianceDesign {
    /// `[sample][column]` design matrix over retained samples.
    design: Vec<Vec<f64>>,
    /// Per-fixed-factor column ranges.
    fixed_factors: Vec<FixedFactor>,
    /// Retained sample IDs in row order.
    kept_sample_ids: Vec<String>,
    /// Per-retained-sample random-group label; `None` when
    /// `random_group` is unset.
    group_labels: Option<Vec<String>>,
}

fn build_variance_design(
    samples: &[atman_core::Sample],
    extras: &BTreeMap<String, BTreeMap<String, String>>,
    fixed_terms: &[String],
    random_group: Option<&str>,
    samples_order: &[String],
) -> Result<VarianceDesign> {
    let sample_by_id: BTreeMap<&str, &atman_core::Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    // For each fixed term, determine whether it's numeric or
    // categorical by probing values across samples_order.
    fn value_of(
        s: &atman_core::Sample,
        extras: &BTreeMap<String, String>,
        term: &str,
    ) -> Option<String> {
        match term {
            "condition" => s.condition.clone(),
            "is_control" => Some(if s.is_control { "1".into() } else { "0".into() }),
            "sample_type" => s.sample_type.clone(),
            "ingest_order" => Some(s.ingest_order.to_string()),
            "subject_id" => s.subject_id.clone(),
            other => extras.get(other).cloned(),
        }
    }

    // Pass 1: collect all values per term; decide type (numeric if
    // every non-empty value parses as f64, else categorical).
    struct TermInfo {
        name: String,
        numeric: bool,
        // For categorical: alphabetically sorted non-reference levels.
        levels: Vec<String>,
    }
    let mut term_info: Vec<TermInfo> = Vec::new();
    for term in fixed_terms {
        let mut seen: BTreeSet<String> = BTreeSet::new();
        let mut all_numeric = true;
        for sid in samples_order {
            let s = match sample_by_id.get(sid.as_str()) {
                Some(s) => *s,
                None => continue,
            };
            let extras_for_s = extras.get(sid).cloned().unwrap_or_default();
            if let Some(v) = value_of(s, &extras_for_s, term) {
                if v.is_empty() {
                    continue;
                }
                if v.parse::<f64>().is_err() {
                    all_numeric = false;
                }
                seen.insert(v);
            }
        }
        let levels: Vec<String> = seen.into_iter().collect();
        if levels.is_empty() {
            bail!("fixed term {:?} has no observed values", term);
        }
        term_info.push(TermInfo {
            name: term.clone(),
            numeric: all_numeric,
            levels,
        });
    }

    // Pass 2: determine column layout. Column 0 is intercept. Each
    // numeric term is 1 column. Each categorical term drops the
    // alphabetically-first level as reference and uses `k-1` columns.
    let mut fixed_factors: Vec<FixedFactor> = Vec::new();
    let mut total_cols = 1; // intercept
                            // Intercept isn't a named factor row in the output; factor rows
                            // are for the user-requested terms.
    for t in &term_info {
        let width = if t.numeric {
            1
        } else {
            t.levels.len().saturating_sub(1)
        };
        if width == 0 {
            bail!(
                "fixed term {:?} has only one level; cannot contribute variance",
                t.name
            );
        }
        fixed_factors.push(FixedFactor {
            name: t.name.clone(),
            columns: total_cols..(total_cols + width),
        });
        total_cols += width;
    }

    let mut design: Vec<Vec<f64>> = Vec::new();
    let mut kept: Vec<String> = Vec::new();
    let mut groups: Vec<String> = Vec::new();
    for sid in samples_order {
        let s = match sample_by_id.get(sid.as_str()) {
            Some(s) => *s,
            None => continue,
        };
        let extras_for_s = extras.get(sid).cloned().unwrap_or_default();
        let mut row = vec![0.0_f64; total_cols];
        row[0] = 1.0; // intercept
        let mut complete = true;
        for (ti, info) in term_info.iter().enumerate() {
            let factor = &fixed_factors[ti];
            let v = value_of(s, &extras_for_s, &info.name).unwrap_or_default();
            if v.is_empty() {
                complete = false;
                break;
            }
            if info.numeric {
                let parsed: f64 = match v.parse() {
                    Ok(x) => x,
                    Err(_) => {
                        complete = false;
                        break;
                    }
                };
                row[factor.columns.start] = parsed;
            } else {
                // Reference level = levels[0]; other levels get a 1 in
                // their respective column.
                let ref_level = &info.levels[0];
                if v == *ref_level {
                    // All zeros is correct.
                } else {
                    let position = info.levels[1..].iter().position(|l| l == &v);
                    match position {
                        Some(pos) => row[factor.columns.start + pos] = 1.0,
                        None => {
                            complete = false;
                            break;
                        }
                    }
                }
            }
        }
        if !complete {
            continue;
        }
        if let Some(rg) = random_group {
            let v = value_of(s, &extras_for_s, rg).unwrap_or_default();
            if v.is_empty() {
                continue;
            }
            groups.push(v);
        }
        design.push(row);
        kept.push(sid.clone());
    }
    let group_labels = if random_group.is_some() {
        Some(groups)
    } else {
        None
    };
    Ok(VarianceDesign {
        design,
        fixed_factors,
        kept_sample_ids: kept,
        group_labels,
    })
}

/// 12-decimal float formatter for variance-decomposition outputs so
/// parity tests against R lm + Type III can assert sub-1e-8 drift.
/// The shared `format_float` truncates to 6 decimals which is too
/// coarse for parity on statistics that run into the hundreds (F,
/// SS) when scaled.
fn format_stat(v: f64) -> String {
    if !v.is_finite() {
        if v.is_nan() {
            "NaN".into()
        } else if v > 0.0 {
            "Inf".into()
        } else {
            "-Inf".into()
        }
    } else if v == 0.0 {
        "0".into()
    } else {
        format!("{v:.12}")
    }
}

fn write_variance_tsv(
    path: &Path,
    rows: &[VarianceRow],
    factors: &[FixedFactor],
    random_group: Option<&str>,
) -> Result<()> {
    let mut header = String::from("archetype_id\ttotal_var\tvar_residual");
    if let Some(g) = random_group {
        header.push_str(&format!("\tvar_random_{g}\ticc_random_{g}"));
    }
    for f in factors {
        // Type III suite per factor: SS_III, F, df_num, df_den, p,
        // plus the per-coefficient max|t|/min p summary for quick
        // inspection.
        header.push_str(&format!(
            "\tss_type3_{name}\tf_statistic_{name}\tdf_num_{name}\tdf_den_{name}\tp_value_{name}\tmax_abs_t_{name}\tmin_p_{name}",
            name = f.name
        ));
    }
    header.push_str("\tskip_reason\n");

    let mut buf = header;
    for r in rows {
        buf.push_str(&r.archetype_id);
        buf.push('\t');
        buf.push_str(&format_stat(r.total_var));
        buf.push('\t');
        buf.push_str(&format_stat(r.var_residual));
        if random_group.is_some() {
            buf.push('\t');
            buf.push_str(&format_stat(r.var_random.unwrap_or(f64::NAN)));
            buf.push('\t');
            buf.push_str(&format_stat(r.icc_random.unwrap_or(f64::NAN)));
        }
        for (i, _) in factors.iter().enumerate() {
            let f = &r.per_factor[i];
            buf.push('\t');
            buf.push_str(&format_stat(f.ss_type3));
            buf.push('\t');
            buf.push_str(&format_stat(f.f_statistic));
            buf.push('\t');
            buf.push_str(&f.df_num.to_string());
            buf.push('\t');
            buf.push_str(&format_stat(f.df_den));
            buf.push('\t');
            buf.push_str(&format_stat(f.p_value));
            buf.push('\t');
            buf.push_str(&format_stat(f.max_abs_t));
            buf.push('\t');
            buf.push_str(&format_stat(f.min_p));
        }
        buf.push('\t');
        buf.push_str(r.skip_reason.as_deref().unwrap_or(""));
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}
