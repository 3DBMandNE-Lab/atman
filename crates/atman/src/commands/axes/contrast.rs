//! `atman axes contrast`: disease modulation of subject-level scores.
//!
//! Per manifest contrast × score column × design: unadjusted statistics on
//! every subject with a score (means, pooled-SD Cohen d with large-sample
//! CI, Welch t, Mann–Whitney AUC), an OLS fit of the score on the design
//! (complete-case over the design covariates) with the `case` coefficient,
//! its t-based CI and p, BH q within (family, design), and an optional
//! subject-level bootstrap of the adjusted coefficient (resampling case and
//! control subjects separately, preserving group sizes; one SplitMix64
//! stream per contrast seeded by `derive_sub_seed(seed, contrast_index)`).

use anyhow::{bail, Context as _, Result};
use atman_core::contrast::{auc_mann_whitney, cohen_d_ci, cohen_d_pooled, percentile, sample_sd};
use atman_core::de::{bh_fdr, ols, omnibus_f_test, welch_t, OlsOutcome, PairedTResult};
use atman_core::stats::mean;
use atman_core::{derive_sub_seed, SplitMix64};
use clap::Args as ClapArgs;
use serde_json::json;
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::collections::{BTreeMap, HashSet};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use super::manifest::{read_contrast_manifest, select_rows};
use super::table::ScoreTable;
use super::{fmt_opt, load_frame, parse_cols, resample_within_groups, Context};
use crate::design::{build_design, parse_design, DegenerateFactor, DesignMatrix, DesignSpec, Term};
use crate::io::{atomic_write, hash_labeled_inputs, sidecar_path_for, write_run_sidecar};

#[derive(ClapArgs, Debug)]
pub struct ContrastArgs {
    /// Wide score table (`sample_id`, `cohort`, `condition`, `is_control`, scores).
    #[arg(long)]
    scores: PathBuf,
    /// Comma-separated `cohort=dir` canonical directories whose `samples.tsv`
    /// supply covariates.
    #[arg(long)]
    cohort_dirs: Option<String>,
    /// Extra `sample_id`-keyed covariate TSVs (repeatable).
    #[arg(long, action = clap::ArgAction::Append)]
    covariates_tsv: Vec<PathBuf>,
    /// Contrast manifest TSV (`label, cohort, case, control, family[, condition_col, subset]`).
    #[arg(long)]
    manifest: PathBuf,
    /// Comma-separated score columns to test.
    #[arg(long)]
    score_cols: String,
    /// Semicolon-separated designs, each containing `case`,
    /// e.g. `~ case; ~ case + z(age) + sex`.
    #[arg(long, default_value = "~ case")]
    designs: String,
    /// Subject-level bootstrap resamples of the adjusted coefficient (0 = off).
    #[arg(long, default_value_t = 0)]
    n_bootstrap: usize,
    /// Bootstrap seed.
    #[arg(long, default_value_t = 1)]
    seed: u64,
    /// Confidence level for the OLS and bootstrap intervals.
    #[arg(long, default_value_t = 0.95)]
    ci: f64,
    /// Output TSV.
    #[arg(long)]
    output: PathBuf,
    /// Optional per-covariate coefficient TSV.
    #[arg(long)]
    output_covariates: Option<PathBuf>,
    /// Optional long TSV of every bootstrap coefficient.
    #[arg(long)]
    output_bootstrap: Option<PathBuf>,
    /// Categorical design column(s) to test jointly (F-test over their
    /// one-hot columns), reported as `omnibus_f` rows in `--output-covariates`.
    #[arg(long, action = clap::ArgAction::Append)]
    omnibus_factor: Vec<String>,
}

#[derive(Debug, Clone)]
struct Row {
    label: String,
    cohort: String,
    family: String,
    score: String,
    design: String,
    n_case: usize,
    n_control: usize,
    mean_case: Option<f64>,
    mean_control: Option<f64>,
    mean_diff: Option<f64>,
    cohen_d: Option<f64>,
    d_ci_lo: Option<f64>,
    d_ci_hi: Option<f64>,
    welch_t: Option<f64>,
    welch_df: Option<f64>,
    welch_p: Option<f64>,
    welch_q: Option<f64>,
    auc: Option<f64>,
    n_fit: usize,
    beta: Option<f64>,
    se: Option<f64>,
    ci_lo: Option<f64>,
    ci_hi: Option<f64>,
    t: Option<f64>,
    p: Option<f64>,
    q: Option<f64>,
    beta_over_resid_sd: Option<f64>,
    boot_n: Option<usize>,
    boot_mean: Option<f64>,
    boot_sd: Option<f64>,
    boot_ci_lo: Option<f64>,
    boot_ci_hi: Option<f64>,
    boot_frac_gt0: Option<f64>,
    boot_n_skipped: Option<usize>,
}

struct CovRow {
    label: String,
    cohort: String,
    score: String,
    design: String,
    kind: &'static str,
    term: String,
    beta: Option<f64>,
    se: Option<f64>,
    t: Option<f64>,
    f: Option<f64>,
    df_num: Option<f64>,
    df_den: Option<f64>,
    p: f64,
    n_fit: usize,
}

struct BootRow {
    label: String,
    score: String,
    design: String,
    resample: usize,
    beta: f64,
}

/// OLS fit with the coefficient of interest at `col` pulled out.
struct Fit {
    beta: f64,
    se: f64,
    t: f64,
    p: f64,
    sigma2: f64,
    df: f64,
    n: usize,
    all_beta: Vec<f64>,
    all_se: Vec<f64>,
    all_t: Vec<f64>,
    all_p: Vec<f64>,
}

fn fit_ols(rows: &[Vec<f64>], y: &[f64], col: usize) -> Option<Fit> {
    match ols(rows, y, 2) {
        OlsOutcome::Computed(f) => Some(Fit {
            beta: f.beta[col],
            se: f.se[col],
            t: f.t[col],
            p: f.p_value[col],
            sigma2: f.sigma2,
            df: f.df,
            n: f.n,
            all_beta: f.beta,
            all_se: f.se,
            all_t: f.t,
            all_p: f.p_value,
        }),
        OlsOutcome::Skipped { .. } => None,
    }
}

/// Design over the case ∪ control rows of one contrast.
fn design_for(
    spec: &DesignSpec,
    ctx: &Context,
    case: &HashSet<usize>,
    control: &HashSet<usize>,
) -> Result<DesignMatrix> {
    build_design(spec, ctx.table.len(), &|i, col| ctx.text(i, col), &|i| {
        if case.contains(&i) {
            Some(1.0)
        } else if control.contains(&i) {
            Some(0.0)
        } else {
            None
        }
    })
}

/// Design over an explicit list of table rows (used by the bootstrap so
/// `z()` is re-standardized on the resample).
fn design_for_rows(
    spec: &DesignSpec,
    ctx: &Context,
    rows: &[usize],
    case: &HashSet<usize>,
) -> Result<DesignMatrix> {
    build_design(spec, rows.len(), &|k, col| ctx.text(rows[k], col), &|k| {
        Some(if case.contains(&rows[k]) { 1.0 } else { 0.0 })
    })
}

pub fn run(args: ContrastArgs) -> Result<()> {
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
    let score_cols = parse_cols(&args.score_cols);
    if score_cols.is_empty() {
        bail!("--score-cols is empty");
    }
    for s in &score_cols {
        table.numeric_column(s)?;
    }
    let designs: Vec<DesignSpec> = args
        .designs
        .split(';')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|d| parse_design(d, Some("case")))
        .collect::<Result<_>>()?;
    if designs.is_empty() {
        bail!("--designs is empty");
    }

    let mut rows: Vec<Row> = Vec::new();
    let mut cov_rows: Vec<CovRow> = Vec::new();
    let mut boot_rows: Vec<BootRow> = Vec::new();

    for (ci_idx, c) in contrasts.iter().enumerate() {
        let sel = select_rows(c, &ctx)?;
        let case: HashSet<usize> = sel.case.iter().copied().collect();
        let control: HashSet<usize> = sel.control.iter().copied().collect();
        let mut rng = SplitMix64::new(derive_sub_seed(args.seed, ci_idx));

        for spec in &designs {
            let dm = design_for(spec, &ctx, &case, &control)
                .with_context(|| format!("contrast {:?} design {:?}", c.label, spec.formula))?;
            let ind_col = dm.indicator_col.expect("design has case");
            let kept_case: Vec<usize> = dm
                .kept
                .iter()
                .copied()
                .filter(|i| case.contains(i))
                .collect();
            let kept_control: Vec<usize> = dm
                .kept
                .iter()
                .copied()
                .filter(|i| control.contains(i))
                .collect();
            if kept_case.is_empty() || kept_control.is_empty() {
                bail!(
                    "contrast {:?} design {:?}: no complete-case rows in one group",
                    c.label,
                    spec.formula
                );
            }
            let draws: Vec<Vec<usize>> = (0..args.n_bootstrap)
                .map(|_| resample_within_groups(&mut rng, &[&kept_case, &kept_control]))
                .collect();

            for score in &score_cols {
                let values = table.numeric_column(score)?;
                let a: Vec<f64> = sel.case.iter().filter_map(|&i| values[i]).collect();
                let b: Vec<f64> = sel.control.iter().filter_map(|&i| values[i]).collect();
                let cohen_d = cohen_d_pooled(&a, &b);
                let (d_ci_lo, d_ci_hi) = cohen_d
                    .and_then(|d| cohen_d_ci(d, a.len(), b.len()))
                    .map(|(lo, hi)| (Some(lo), Some(hi)))
                    .unwrap_or((None, None));
                let (welch_t_v, welch_df, welch_p) = match welch_t(&a, &b, 2) {
                    PairedTResult::Computed { t, df, p_value, .. } => {
                        (Some(t), Some(df), Some(p_value))
                    }
                    PairedTResult::Skipped { .. } => (None, None, None),
                };
                let mean_case = if a.is_empty() { None } else { Some(mean(&a)) };
                let mean_control = if b.is_empty() { None } else { Some(mean(&b)) };

                // OLS on design rows with a finite score.
                let mut x: Vec<Vec<f64>> = Vec::new();
                let mut y: Vec<f64> = Vec::new();
                for (k, &i) in dm.kept.iter().enumerate() {
                    if let Some(v) = values[i] {
                        x.push(dm.rows[k].clone());
                        y.push(v);
                    }
                }
                let fit = fit_ols(&x, &y, ind_col);
                let (ci_lo, ci_hi) = match &fit {
                    Some(f) => match StudentsT::new(0.0, 1.0, f.df) {
                        Ok(dist) => {
                            let crit = dist.inverse_cdf(0.5 + args.ci / 2.0);
                            (Some(f.beta - crit * f.se), Some(f.beta + crit * f.se))
                        }
                        Err(_) => (None, None),
                    },
                    None => (None, None),
                };
                if let Some(f) = &fit {
                    for (col, label) in dm.labels.iter().enumerate() {
                        if col == 0 || col == ind_col {
                            continue;
                        }
                        cov_rows.push(CovRow {
                            label: c.label.clone(),
                            cohort: c.cohort.clone(),
                            score: score.clone(),
                            design: spec.formula.clone(),
                            kind: "coef",
                            term: label.clone(),
                            beta: Some(f.all_beta[col]),
                            se: Some(f.all_se[col]),
                            t: Some(f.all_t[col]),
                            f: None,
                            df_num: None,
                            df_den: None,
                            p: f.all_p[col],
                            n_fit: f.n,
                        });
                    }
                    for factor in &args.omnibus_factor {
                        let Some(t_idx) = spec
                            .terms
                            .iter()
                            .position(|t| matches!(t, Term::Column(c) if c == factor))
                        else {
                            continue;
                        };
                        let factor_cols: Vec<usize> = dm
                            .term_index
                            .iter()
                            .enumerate()
                            .filter(|(_, t)| **t == Some(t_idx))
                            .map(|(i, _)| i)
                            .collect();
                        if factor_cols.is_empty() {
                            continue;
                        }
                        if let Some(om) =
                            omnibus_f_test(&x, &f.all_beta, &factor_cols, f.sigma2, f.df)
                        {
                            cov_rows.push(CovRow {
                                label: c.label.clone(),
                                cohort: c.cohort.clone(),
                                score: score.clone(),
                                design: spec.formula.clone(),
                                kind: "omnibus_f",
                                term: factor.clone(),
                                beta: None,
                                se: None,
                                t: None,
                                f: Some(om.f_statistic),
                                df_num: Some(om.df_num as f64),
                                df_den: Some(om.df_den),
                                p: om.p_value,
                                n_fit: f.n,
                            });
                        }
                    }
                }

                // Bootstrap. A replicate whose resample leaves a categorical
                // term with a single level (or a singular fit) is skipped and
                // counted in boot_n_skipped rather than aborting the run.
                let mut betas: Vec<f64> = Vec::with_capacity(draws.len());
                let mut n_skipped = 0usize;
                for (b_idx, draw) in draws.iter().enumerate() {
                    let bdm = match design_for_rows(spec, &ctx, draw, &case) {
                        Ok(d) => d,
                        Err(e) if e.downcast_ref::<DegenerateFactor>().is_some() => {
                            n_skipped += 1;
                            continue;
                        }
                        Err(e) => return Err(e),
                    };
                    let bcol = bdm.indicator_col.expect("design has case");
                    let mut bx = Vec::new();
                    let mut by = Vec::new();
                    for (k, &pos) in bdm.kept.iter().enumerate() {
                        if let Some(v) = values[draw[pos]] {
                            bx.push(bdm.rows[k].clone());
                            by.push(v);
                        }
                    }
                    match fit_ols(&bx, &by, bcol) {
                        Some(bf) => {
                            betas.push(bf.beta);
                            boot_rows.push(BootRow {
                                label: c.label.clone(),
                                score: score.clone(),
                                design: spec.formula.clone(),
                                resample: b_idx,
                                beta: bf.beta,
                            });
                        }
                        None => n_skipped += 1,
                    }
                }
                let (boot_n, boot_mean, boot_sd, boot_ci_lo, boot_ci_hi, boot_frac_gt0) =
                    if args.n_bootstrap > 0 && !betas.is_empty() {
                        let mut sorted = betas.clone();
                        sorted
                            .sort_by(|p, q| p.partial_cmp(q).unwrap_or(std::cmp::Ordering::Equal));
                        let alpha = (1.0 - args.ci) / 2.0;
                        (
                            Some(betas.len()),
                            Some(mean(&betas)),
                            sample_sd(&betas),
                            percentile(&sorted, alpha),
                            percentile(&sorted, 1.0 - alpha),
                            Some(
                                betas.iter().filter(|v| **v > 0.0).count() as f64
                                    / betas.len() as f64,
                            ),
                        )
                    } else {
                        (None, None, None, None, None, None)
                    };

                rows.push(Row {
                    label: c.label.clone(),
                    cohort: c.cohort.clone(),
                    family: c.family.clone(),
                    score: score.clone(),
                    design: spec.formula.clone(),
                    n_case: a.len(),
                    n_control: b.len(),
                    mean_case,
                    mean_control,
                    mean_diff: mean_case.zip(mean_control).map(|(x, y)| x - y),
                    cohen_d,
                    d_ci_lo,
                    d_ci_hi,
                    welch_t: welch_t_v,
                    welch_df,
                    welch_p,
                    welch_q: None,
                    auc: auc_mann_whitney(&a, &b),
                    n_fit: fit.as_ref().map(|f| f.n).unwrap_or(0),
                    beta: fit.as_ref().map(|f| f.beta),
                    se: fit.as_ref().map(|f| f.se),
                    ci_lo,
                    ci_hi,
                    t: fit.as_ref().map(|f| f.t),
                    p: fit.as_ref().map(|f| f.p),
                    q: None,
                    beta_over_resid_sd: fit.as_ref().map(|f| f.beta / f.sigma2.sqrt()),
                    boot_n,
                    boot_mean,
                    boot_sd,
                    boot_ci_lo,
                    boot_ci_hi,
                    boot_frac_gt0,
                    boot_n_skipped: if args.n_bootstrap > 0 {
                        Some(n_skipped)
                    } else {
                        None
                    },
                });
            }
        }
    }

    // BH within (family, design).
    let mut groups: BTreeMap<(String, String), Vec<usize>> = BTreeMap::new();
    for (i, r) in rows.iter().enumerate() {
        groups
            .entry((r.family.clone(), r.design.clone()))
            .or_default()
            .push(i);
    }
    for idx in groups.values() {
        let ps: Vec<Option<f64>> = idx.iter().map(|&i| rows[i].p).collect();
        for (k, q) in bh_fdr(&ps).into_iter().enumerate() {
            rows[idx[k]].q = q;
        }
        let wps: Vec<Option<f64>> = idx.iter().map(|&i| rows[i].welch_p).collect();
        for (k, q) in bh_fdr(&wps).into_iter().enumerate() {
            rows[idx[k]].welch_q = q;
        }
    }

    write_rows(&args.output, &rows)?;
    let mut outputs = vec![args.output.clone()];
    if let Some(p) = &args.output_covariates {
        write_cov_rows(p, &cov_rows)?;
        outputs.push(p.clone());
    }
    if let Some(p) = &args.output_bootstrap {
        write_boot_rows(p, &boot_rows)?;
        outputs.push(p.clone());
    }
    eprintln!(
        "axes contrast: contrasts={} scores={} designs={} rows={} n_bootstrap={} output={}",
        contrasts.len(),
        score_cols.len(),
        designs.len(),
        rows.len(),
        args.n_bootstrap,
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
        "axes contrast",
        json!({
            "scores": args.scores.display().to_string(),
            "cohort-dirs": args.cohort_dirs,
            "covariates-tsv": args.covariates_tsv.iter().map(|p| p.display().to_string()).collect::<Vec<_>>(),
            "manifest": args.manifest.display().to_string(),
            "score-cols": args.score_cols,
            "designs": args.designs,
            "n-bootstrap": args.n_bootstrap,
            "seed": args.seed,
            "ci": args.ci,
            "output": args.output.display().to_string(),
            "output-covariates": args.output_covariates.as_ref().map(|p| p.display().to_string()),
            "output-bootstrap": args.output_bootstrap.as_ref().map(|p| p.display().to_string()),
            "omnibus-factor": args.omnibus_factor,
        }),
        &inputs,
        &outputs,
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("axes contrast: sidecar={}", sidecar.display());
    Ok(())
}

const HEADER: &str = "label\tcohort\tscore\tdesign\tn_case\tn_control\tmean_case\tmean_control\tmean_diff\tcohen_d\td_ci_lo\td_ci_hi\twelch_t\twelch_df\twelch_p\twelch_q\tauc\tn_fit\tbeta\tse\tci_lo\tci_hi\tt\tp\tq\tbeta_over_resid_sd\tboot_n\tboot_mean\tboot_sd\tboot_ci_lo\tboot_ci_hi\tboot_frac_gt0\tboot_n_skipped\n";

fn write_rows(path: &Path, rows: &[Row]) -> Result<()> {
    let mut buf = String::from(HEADER);
    for r in rows {
        let cells: Vec<String> = vec![
            r.label.clone(),
            r.cohort.clone(),
            r.score.clone(),
            r.design.clone(),
            r.n_case.to_string(),
            r.n_control.to_string(),
            fmt_opt(r.mean_case),
            fmt_opt(r.mean_control),
            fmt_opt(r.mean_diff),
            fmt_opt(r.cohen_d),
            fmt_opt(r.d_ci_lo),
            fmt_opt(r.d_ci_hi),
            fmt_opt(r.welch_t),
            fmt_opt(r.welch_df),
            fmt_opt(r.welch_p),
            fmt_opt(r.welch_q),
            fmt_opt(r.auc),
            r.n_fit.to_string(),
            fmt_opt(r.beta),
            fmt_opt(r.se),
            fmt_opt(r.ci_lo),
            fmt_opt(r.ci_hi),
            fmt_opt(r.t),
            fmt_opt(r.p),
            fmt_opt(r.q),
            fmt_opt(r.beta_over_resid_sd),
            r.boot_n.map(|n| n.to_string()).unwrap_or_default(),
            fmt_opt(r.boot_mean),
            fmt_opt(r.boot_sd),
            fmt_opt(r.boot_ci_lo),
            fmt_opt(r.boot_ci_hi),
            fmt_opt(r.boot_frac_gt0),
            r.boot_n_skipped.map(|n| n.to_string()).unwrap_or_default(),
        ];
        buf.push_str(&cells.join("\t"));
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn write_cov_rows(path: &Path, rows: &[CovRow]) -> Result<()> {
    let mut buf = String::from(
        "label\tcohort\tscore\tdesign\tkind\tterm\tbeta\tse\tt\tf\tdf_num\tdf_den\tp\tn_fit\n",
    );
    for r in rows {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.label,
            r.cohort,
            r.score,
            r.design,
            r.kind,
            r.term,
            fmt_opt(r.beta),
            fmt_opt(r.se),
            fmt_opt(r.t),
            fmt_opt(r.f),
            fmt_opt(r.df_num),
            fmt_opt(r.df_den),
            fmt_opt(Some(r.p)),
            r.n_fit
        ));
    }
    atomic_write(path, buf.as_bytes())
}

fn write_boot_rows(path: &Path, rows: &[BootRow]) -> Result<()> {
    let mut buf = String::from("label\tscore\tdesign\tresample\tbeta\n");
    for r in rows {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            r.label,
            r.score,
            r.design,
            r.resample,
            fmt_opt(Some(r.beta))
        ));
    }
    atomic_write(path, buf.as_bytes())
}
