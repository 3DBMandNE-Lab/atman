//! `--test msqrob`: peptide-level ridge mixed model per protein.
//! One linear mixed model fit over each protein's peptide-level
//! observations in a two-group comparison, with random intercept per
//! peptide and ridge penalty on non-intercept fixed effects. After
//! per-protein fitting, residual variances are empirically Bayes
//! shrunk across proteins (`squeeze_variance`) and the test statistic
//! is recomputed under the shrunk degrees of freedom.
//!
//! When `--adjust-for` covariates are provided, they are appended to
//! the design matrix after the group-b indicator column, matching the
//! convention in `limma.rs`. The condition contrast is `[0, -1, 0, …, 0]`
//! so covariate columns do not enter the reported effect.

use anyhow::{Context, Result};
use atman_core::msqrob::{fit_msqrob, squeeze_variance, MsqrobFit, MsqrobOutcome};
use atman_core::Sample;
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::collections::{BTreeMap, BTreeSet, HashMap};

use super::{limma::ExternalCovariates, Args};
use crate::io::{read_peptide_measurements, read_peptides, DeReportRow, DeResultRow};

/// Dispatch for `--test msqrob`. One ridge-regularized linear mixed
/// model per protein, fit over its peptide-level observations in the
/// two-group comparison. The base design matrix is `[intercept, group_b]`;
/// when `external_covariates` is provided, one column per covariate is
/// appended in deterministic `BTreeMap` order. The reported effect is
/// `mean_a − mean_b`, matching the sign convention used by the paired-t,
/// welch-t, and ols dispatches (contrast = −β_group_b).
pub(super) fn run_msqrob(
    args: &Args,
    samples: &[Sample],
    proteins: &[atman_core::ProteinIdentity],
    comparisons: &[(String, String)],
    external_covariates: Option<&ExternalCovariates>,
) -> Result<(Vec<DeResultRow>, Vec<DeReportRow>)> {
    let pep_meas_path = args
        .peptide_measurements
        .as_ref()
        .expect("validated above: peptide_measurements required for msqrob");
    let pep_meta_path = args
        .peptide_metadata
        .as_ref()
        .expect("validated above: peptide_metadata required for msqrob");
    let peptides = read_peptides(pep_meta_path)
        .with_context(|| format!("reading peptide metadata from {:?}", pep_meta_path))?;
    let pep_measurements = read_peptide_measurements(pep_meas_path)
        .with_context(|| format!("reading peptide measurements from {:?}", pep_meas_path))?;

    // Validated and resolved in `super::run` so the sidecar can record
    // the penalty that actually applied; `auto` is refused there.
    let ridge_lambda: f64 = args.ridge_lambda.trim().parse::<f64>().with_context(|| {
        format!(
            "--ridge-lambda must be a number, got {:?}",
            args.ridge_lambda
        )
    })?;

    // assay_id → (panel, gene_symbol, uniprot).
    let mut protein_meta: BTreeMap<String, (String, String, Vec<String>)> = BTreeMap::new();
    for p in proteins {
        let panel = p.panel.clone().unwrap_or_default();
        let gene = p.gene_symbol.clone().unwrap_or_default();
        protein_meta
            .entry(p.assay_id.0.clone())
            .or_insert((panel, gene, p.uniprot.clone()));
    }

    // peptide_id → assay_id lookup (parent protein).
    let peptide_parent: BTreeMap<String, String> = peptides
        .iter()
        .map(|p| (p.peptide_id.clone(), p.assay_id.0.clone()))
        .collect();
    // assay_id → ordered, deduped peptide ids, for deterministic iteration.
    let mut peptides_by_protein: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for p in &peptides {
        peptides_by_protein
            .entry(p.assay_id.0.clone())
            .or_default()
            .insert(p.peptide_id.clone());
    }

    // (sample_id, peptide_id) → abundance, QC-masked and non-finite rows
    // already dropped via `effective_abundance`.
    let mut pep_cells: HashMap<(String, String), f64> = HashMap::new();
    for m in &pep_measurements {
        if let Some(v) = m.effective_abundance() {
            pep_cells.insert((m.sample_id.clone(), m.peptide_id.clone()), v);
        }
    }
    let _ = peptide_parent; // kept for future peptide-level joins; silence unused warn.

    let sample_by_id: HashMap<&str, &Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    let mut all_result_rows: Vec<DeResultRow> = Vec::new();
    let mut all_report_rows: Vec<DeReportRow> = Vec::new();

    for (a, b) in comparisons {
        let comparison_label = format!("{}-{}", a, b);

        // Sample filter: non-control, condition ∈ {a, b}; assign group_b
        // indicator (0 for a, 1 for b).
        let mut sample_group: Vec<(String, f64)> = Vec::new();
        for s in samples {
            if s.is_control {
                continue;
            }
            let cond = match s.condition.as_deref() {
                Some(c) => c,
                None => continue,
            };
            let group_b = if cond == a {
                0.0
            } else if cond == b {
                1.0
            } else {
                continue;
            };
            sample_group.push((s.sample_id.clone(), group_b));
        }
        if sample_group.is_empty() {
            continue;
        }

        // Per-panel report accumulator.
        let mut per_panel: BTreeMap<String, PanelAcc> = BTreeMap::new();
        // Two-pass layout: (i) collect each protein's fit or skip reason
        // into parallel Vecs; (ii) variance-squeeze Computed fits across
        // all proteins (empirical-Bayes shrinkage matching the
        // `squeezeVarRob` step in msqrob2); (iii) emit DeResultRows
        // using the shrunk SE / t / p / df.
        let mut fit_contexts: Vec<MsqrobFitContext> = Vec::new();
        let mut computed_fits: Vec<MsqrobFit> = Vec::new();
        let mut computed_ctx_index: Vec<usize> = Vec::new();
        let mut skipped_rows: Vec<(String, MsqrobRow)> = Vec::new();

        // Determine covariate column names in stable order from the first
        // sample's covariate map. All samples are guaranteed to have the
        // same names because join_external_covariates validated this at
        // call-site (in mod.rs). We compute the name list once per
        // comparison and reuse it for every protein.
        let cov_names: Vec<String> = if let Some(ext) = external_covariates {
            sample_group
                .first()
                .and_then(|(sid, _)| ext.get(sid))
                .map(|m| m.keys().cloned().collect())
                .unwrap_or_default()
        } else {
            Vec::new()
        };

        for (assay_id, peptide_ids) in &peptides_by_protein {
            let (panel, gene, uniprot) = protein_meta
                .get(assay_id)
                .cloned()
                .unwrap_or_else(|| (String::new(), String::new(), vec![]));
            // Assign dense peptide indices for this protein based on
            // sorted iteration order.
            let mut pep_index: BTreeMap<String, usize> = BTreeMap::new();
            for (idx, pid) in peptide_ids.iter().enumerate() {
                pep_index.insert(pid.clone(), idx);
            }

            let mut design: Vec<Vec<f64>> = Vec::new();
            let mut y: Vec<f64> = Vec::new();
            let mut peptide_obs: Vec<usize> = Vec::new();
            let mut used_peptides: BTreeSet<String> = BTreeSet::new();
            for (sample_id, group_b) in &sample_group {
                if !sample_by_id.contains_key(sample_id.as_str()) {
                    continue;
                }
                for pid in peptide_ids {
                    let key = (sample_id.clone(), pid.clone());
                    if let Some(&abund) = pep_cells.get(&key) {
                        // Base design row: [intercept, group_b_indicator].
                        // Append one column per external covariate (same order
                        // as cov_names, which is stable BTreeMap key order).
                        let mut row = vec![1.0, *group_b];
                        if let Some(ext) = external_covariates {
                            if let Some(covs) = ext.get(sample_id) {
                                for name in &cov_names {
                                    row.push(*covs.get(name).unwrap_or(&0.0));
                                }
                            }
                        }
                        design.push(row);
                        y.push(abund);
                        peptide_obs.push(pep_index[pid]);
                        used_peptides.insert(pid.clone());
                    }
                }
            }

            // Group-wise raw means (on log-scale peptide abundances).
            let (mean_a, mean_b) = msqrob_raw_group_means(&y, &design);

            let outcome = fit_msqrob(
                &design,
                &y,
                &peptide_obs,
                ridge_lambda,
                args.min_peptides,
                args.min_pairs,
                args.robust,
            );

            match outcome {
                MsqrobOutcome::Computed(fit) => {
                    let ctx = MsqrobFitContext {
                        panel: panel.clone(),
                        assay_id: assay_id.clone(),
                        gene: gene.clone(),
                        uniprot: uniprot.clone(),
                        mean_a,
                        mean_b,
                    };
                    computed_ctx_index.push(fit_contexts.len());
                    fit_contexts.push(ctx);
                    computed_fits.push(fit);
                }
                MsqrobOutcome::Skipped { reason, n } => {
                    let acc = per_panel.entry(panel.clone()).or_default();
                    acc.n_skipped += 1;
                    let reason_str = match reason {
                        atman_core::SkipReason::InsufficientPairs => {
                            if used_peptides.len() < args.min_peptides {
                                "insufficient_peptides"
                            } else {
                                "insufficient_observations"
                            }
                        }
                        atman_core::SkipReason::ZeroVariance => "singular_design",
                        atman_core::SkipReason::NonFiniteInput => "non_finite_input",
                    };
                    skipped_rows.push((
                        panel.clone(),
                        MsqrobRow {
                            de: DeResultRow {
                                panel: panel.clone(),
                                assay_id: assay_id.clone(),
                                gene_symbol: gene.clone(),
                                uniprot: uniprot.join(","),
                                comparison: comparison_label.clone(),
                                n_pairs: n,
                                mean_a: None,
                                mean_b: None,
                                mean_diff: None,
                                t: None,
                                df: None,
                                p_value: None,
                                bh_q: None,
                                effect_size: None,
                                effect_size_method: "msqrob-ridge".to_string(),
                                ci_low: None,
                                ci_high: None,
                                wilcoxon_p: None,
                                wilcoxon_method: String::new(),
                                median_diff: None,
                                trimmed_mean_diff: None,
                                skip_reason: reason_str.to_string(),
                                s2_trend: None,
                                s2_prior: None,
                                s2_posterior: None,
                                df_prior: None,
                                df_total: None,
                                f_statistic: None,
                                f_p_value: None,
                                f_bh_q: None,
                                lfc_threshold: None,
                                n_peptides_observed: Some(used_peptides.len()),
                                peptide_variance_ratio: None,
                                ridge_lambda: Some(ridge_lambda),
                                method: "msqrob".into(),
                                posthoc_method: String::new(),
                                posthoc_p: None,
                                posthoc_adj_p: None,
                                n_a: None,
                                n_b: None,
                            },
                        },
                    ));
                }
            }
        }

        // Empirical-Bayes variance squeeze across all Computed fits in
        // this comparison. Matches msqrob2's `squeezeVarRob` step; no-op
        // when fewer than two fits produced a finite positive variance.
        let squeeze = squeeze_variance(&mut computed_fits);
        let (prior_df, prior_s2) = match squeeze {
            Some((df_prior, s2_prior)) => (Some(df_prior), Some(s2_prior)),
            None => (None, None),
        };

        let mut fits: Vec<(String, MsqrobRow)> = Vec::new();
        for (fit, &ctx_idx) in computed_fits.iter().zip(computed_ctx_index.iter()) {
            let ctx = &fit_contexts[ctx_idx];
            // atman convention: effect = mean_a − mean_b = −β_groupB.
            let effect = -fit.beta[1];
            let se = fit.se[1];
            let t = -fit.t[1];
            let p_value = fit.p_value[1];
            // CI critical value must match the t-based p-value: use the
            // two-sided 97.5% quantile of Student's t at the fit's df,
            // not a fixed 1.96 z-multiplier (which disagrees at small df).
            let t_crit = StudentsT::new(0.0, 1.0, fit.df)
                .ok()
                .map(|d| d.inverse_cdf(0.975))
                .unwrap_or(1.959963985);
            let ci_low = effect - t_crit * se;
            let ci_high = effect + t_crit * se;
            let acc = per_panel.entry(ctx.panel.clone()).or_default();
            acc.n_tests += 1;
            acc.collect_effect(effect);
            fits.push((
                ctx.panel.clone(),
                MsqrobRow {
                    de: DeResultRow {
                        panel: ctx.panel.clone(),
                        assay_id: ctx.assay_id.clone(),
                        gene_symbol: ctx.gene.clone(),
                        uniprot: ctx.uniprot.join(","),
                        comparison: comparison_label.clone(),
                        n_pairs: fit.n,
                        mean_a: Some(ctx.mean_a),
                        mean_b: Some(ctx.mean_b),
                        mean_diff: Some(effect),
                        t: Some(t),
                        df: Some(fit.df),
                        p_value: Some(p_value),
                        bh_q: None,
                        effect_size: Some(effect),
                        effect_size_method: "msqrob-ridge".to_string(),
                        ci_low: Some(ci_low),
                        ci_high: Some(ci_high),
                        wilcoxon_p: None,
                        wilcoxon_method: String::new(),
                        median_diff: None,
                        trimmed_mean_diff: None,
                        skip_reason: String::new(),
                        s2_trend: None,
                        s2_prior: prior_s2,
                        s2_posterior: Some(fit.sigma2),
                        df_prior: prior_df,
                        df_total: Some(fit.df),
                        f_statistic: None,
                        f_p_value: None,
                        f_bh_q: None,
                        lfc_threshold: None,
                        n_peptides_observed: Some(fit.n_peptides),
                        peptide_variance_ratio: Some(fit.peptide_variance_ratio),
                        ridge_lambda: Some(ridge_lambda),
                        method: "msqrob".into(),
                        posthoc_method: String::new(),
                        posthoc_p: None,
                        posthoc_adj_p: None,
                        n_a: None,
                        n_b: None,
                    },
                },
            ));
        }
        fits.extend(skipped_rows);

        // BH-FDR per panel within this comparison.
        let mut by_panel: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (i, (panel, _)) in fits.iter().enumerate() {
            by_panel.entry(panel.clone()).or_default().push(i);
        }
        for indices in by_panel.values() {
            let ps: Vec<Option<f64>> = indices.iter().map(|&i| fits[i].1.de.p_value).collect();
            let qs = atman_core::bh_fdr(&ps);
            for (local, &global) in indices.iter().enumerate() {
                fits[global].1.de.bh_q = qs[local];
                if let Some(q) = qs[local] {
                    let panel = fits[global].0.clone();
                    let acc = per_panel.entry(panel).or_default();
                    if q < args.report_q_strict {
                        acc.n_q_strict += 1;
                    }
                    if q < args.report_q_relaxed {
                        acc.n_q_relaxed += 1;
                    }
                    if acc.min_q.map(|m| q < m).unwrap_or(true) {
                        acc.min_q = Some(q);
                    }
                }
            }
        }

        for (_panel, row) in fits {
            all_result_rows.push(row.de);
        }
        for (panel, acc) in per_panel {
            all_report_rows.push(DeReportRow {
                comparison: comparison_label.clone(),
                panel,
                n_tests: acc.n_tests,
                n_skipped: acc.n_skipped,
                n_q_strict: acc.n_q_strict,
                n_q_relaxed: acc.n_q_relaxed,
                min_q: acc.min_q,
                max_abs_effect: acc.max_abs_effect,
                limma_trend_fallback_used: None,
            });
        }
    }

    Ok((all_result_rows, all_report_rows))
}

#[derive(Default)]
struct PanelAcc {
    n_tests: usize,
    n_skipped: usize,
    n_q_strict: usize,
    n_q_relaxed: usize,
    min_q: Option<f64>,
    max_abs_effect: Option<f64>,
}

impl PanelAcc {
    fn collect_effect(&mut self, effect: f64) {
        let abs_e = effect.abs();
        if self.max_abs_effect.map(|m| abs_e > m).unwrap_or(true) {
            self.max_abs_effect = Some(abs_e);
        }
    }
}

struct MsqrobRow {
    de: DeResultRow,
}

struct MsqrobFitContext {
    panel: String,
    assay_id: String,
    gene: String,
    uniprot: Vec<String>,
    mean_a: f64,
    mean_b: f64,
}

/// Compute raw group means of `y` grouped by the group-b indicator column
/// (col 1) of the design matrix. Returns `(mean_a, mean_b)`.
fn msqrob_raw_group_means(y: &[f64], design: &[Vec<f64>]) -> (f64, f64) {
    let mut sum_a = 0.0;
    let mut n_a = 0usize;
    let mut sum_b = 0.0;
    let mut n_b = 0usize;
    for (row, &yi) in design.iter().zip(y.iter()) {
        if row[1] < 0.5 {
            sum_a += yi;
            n_a += 1;
        } else {
            sum_b += yi;
            n_b += 1;
        }
    }
    let mean_a = if n_a > 0 {
        sum_a / n_a as f64
    } else {
        f64::NAN
    };
    let mean_b = if n_b > 0 {
        sum_b / n_b as f64
    } else {
        f64::NAN
    };
    (mean_a, mean_b)
}
