//! `--test limma`: parametric mean-variance trend with empirical Bayes
//! variance shrinkage (Smyth 2004). Optionally swaps the trend
//! covariate to `log(peptide_count + 1)` (DEqMS, Zhu et al. 2020) when
//! `--peptide-metadata` is provided. Constructs a features × samples
//! `y` matrix per (comparison, panel) and dispatches to
//! [`atman_core::limma::limma_fit`].

use anyhow::{Context, Result};
use atman_core::limma::{limma_fit, LimmaOptions};
use atman_core::{de::bh_fdr, MeasurementRecord, ProteinIdentity, Sample};
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::collections::{BTreeMap, BTreeSet, HashMap};

use super::{build_two_group_design, Args};
use crate::io::{read_peptides, DeReportRow, DeResultRow};

/// External covariates map: sample_id → {covariate_name → value}.
pub(super) type ExternalCovariates = BTreeMap<String, BTreeMap<String, f64>>;

pub(super) fn run_limma(
    args: &Args,
    samples: &[Sample],
    proteins: &[ProteinIdentity],
    measurements: &[MeasurementRecord],
    comparisons: &[(String, String)],
    external_covariates: Option<&ExternalCovariates>,
) -> Result<(Vec<DeResultRow>, Vec<DeReportRow>)> {
    let mut all_result_rows: Vec<DeResultRow> = Vec::new();
    let mut all_report_rows: Vec<DeReportRow> = Vec::new();

    // Optional DEqMS input: peptides.tsv provides peptide→parent
    // assay_id mapping. Count peptides per assay_id once; per-comparison
    // lookups then pull the relevant counts.
    let peptides_per_assay: Option<BTreeMap<String, u32>> =
        if let Some(pep_path) = args.peptide_metadata.as_ref() {
            let peptides = read_peptides(pep_path)
                .with_context(|| format!("reading peptide metadata from {:?}", pep_path))?;
            let mut map: BTreeMap<String, u32> = BTreeMap::new();
            for p in &peptides {
                *map.entry(p.assay_id.0.clone()).or_insert(0) += 1;
            }
            Some(map)
        } else {
            None
        };
    let deqms_active = peptides_per_assay.is_some();

    let effect_size_method = match (args.trend, args.robust, deqms_active) {
        (_, true, true) => "limma-DEqMS-robust-trend",
        (_, false, true) => "limma-DEqMS-trend",
        (true, true, false) => "limma-eBayes-robust-trend",
        (true, false, false) => "limma-eBayes-trend",
        (false, true, false) => "limma-eBayes-robust",
        (false, false, false) => "limma-eBayes",
    };

    // Lookup: (panel, gene_symbol) → (first_assay_id, uniprot) for reporting.
    // Mirrors the lookup built in `fn run`; proteins without a panel or
    // gene_symbol cannot be joined back to measurement rows via cells and
    // are excluded from the feature set.
    let mut gene_meta: BTreeMap<(String, String), (String, Vec<String>)> = BTreeMap::new();
    for p in proteins {
        if let (Some(gene), Some(panel)) = (p.gene_symbol.as_ref(), p.panel.as_ref()) {
            gene_meta
                .entry((panel.clone(), gene.clone()))
                .or_insert_with(|| (p.assay_id.0.clone(), p.uniprot.clone()));
        }
    }

    // Lookup: sample_id → &Sample, so we can gate control status / condition
    // per measurement without scanning `samples` each time.
    let sample_by_id: HashMap<&str, &Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    // Build per-(panel, gene, sample_id) → abundance. QC-masked cells drop
    // out of the map and become NaN in the y matrix for their feature.
    // Also track the set of (panel, gene_symbol) features actually
    // measured so the feature list for each comparison is keyed off
    // measurements (matching welch-t / ols dispatch behaviour).
    let mut cells_by_sample: HashMap<(String, String, String), f64> = HashMap::new();
    let mut measured_features: BTreeSet<(String, String)> = BTreeSet::new();
    for m in measurements {
        let abundance = match m.effective_abundance() {
            Some(a) => a,
            None => continue,
        };
        let s = match sample_by_id.get(m.sample_id.as_str()) {
            Some(s) => *s,
            None => continue,
        };
        if s.is_control {
            continue;
        }
        if s.condition.is_none() {
            continue;
        }
        let panel = match &m.panel {
            Some(p) => p.clone(),
            None => continue,
        };
        let gene = match &m.gene_symbol {
            Some(g) => g.clone(),
            None => continue,
        };
        measured_features.insert((panel.clone(), gene.clone()));
        cells_by_sample.insert((panel, gene, m.sample_id.clone()), abundance);
    }

    for (a, b) in comparisons {
        let comparison_label = format!("{}-{}", a, b);

        // (i) Design matrix: rows = non-control samples with condition in {a, b}.
        // Base columns: [intercept, group_b_indicator].
        // When external_covariates is provided, append one column per
        // covariate (in deterministic BTreeMap order) after the group column.
        let (base_design, sample_ids) = build_two_group_design(samples, a, b);
        if base_design.is_empty() {
            // No samples in either group — skip this comparison silently;
            // no rows, no report.
            continue;
        }
        // Augment with external covariates if provided.
        let design: Vec<Vec<f64>> = if let Some(ext) = external_covariates {
            // Determine covariate names in stable order from the first sample's map.
            // All samples are guaranteed to have the same covariate names because
            // join_external_covariates validated this at call-site.
            let cov_names: Vec<String> = sample_ids
                .first()
                .and_then(|sid| ext.get(sid))
                .map(|m| m.keys().cloned().collect())
                .unwrap_or_default();
            base_design
                .iter()
                .zip(sample_ids.iter())
                .map(|(row, sid)| {
                    let mut augmented = row.clone();
                    if let Some(covs) = ext.get(sid) {
                        for name in &cov_names {
                            augmented.push(*covs.get(name).unwrap_or(&0.0));
                        }
                    }
                    augmented
                })
                .collect()
        } else {
            base_design
        };

        // Feature set = (panel, gene_symbol) pairs observed in
        // measurements that have at least one retained-sample value in
        // this comparison. `gene_meta` lookups fall back to empty
        // (assay_id, uniprot) when a measured gene isn't listed in
        // proteins.tsv — matching the welch-t / ols behaviour. Iteration
        // is deterministic via BTreeSet.
        let mut features: Vec<(String, String, String, Vec<String>)> = Vec::new();
        for (panel, gene) in &measured_features {
            let has_any = sample_ids.iter().any(|sid| {
                cells_by_sample.contains_key(&(panel.clone(), gene.clone(), sid.clone()))
            });
            if !has_any {
                continue;
            }
            let (assay_id, uniprot) = gene_meta
                .get(&(panel.clone(), gene.clone()))
                .cloned()
                .unwrap_or_else(|| (String::new(), vec![]));
            features.push((panel.clone(), gene.clone(), assay_id, uniprot));
        }

        if features.is_empty() {
            // Nothing to fit — emit a singleton report with zero tests so
            // users see the comparison was considered.
            all_report_rows.push(DeReportRow {
                comparison: comparison_label.clone(),
                panel: String::new(),
                n_tests: 0,
                n_skipped: 0,
                n_q_strict: 0,
                n_q_relaxed: 0,
                min_q: None,
                max_abs_effect: None,
                limma_trend_fallback_used: None,
            });
            continue;
        }

        // (ii) y matrix: rows = features, cols = samples in `sample_ids`
        // order. Missing cells → NaN so `limma_fit`'s internal ols skips
        // the feature.
        let mut y: Vec<Vec<f64>> = Vec::with_capacity(features.len());
        for (panel, gene, _, _) in &features {
            let mut row = Vec::with_capacity(sample_ids.len());
            for sid in &sample_ids {
                let v = cells_by_sample
                    .get(&(panel.clone(), gene.clone(), sid.clone()))
                    .copied()
                    .unwrap_or(f64::NAN);
                row.push(v);
            }
            y.push(row);
        }

        // (iii) Single-contrast matrix: test (A − B) direction.
        // Design column 1 is the group-B indicator, so β_1 = mean_b − mean_a;
        // atman's `A-B` comparison semantics are "mean_a minus mean_b"
        // (matching welch-t and paired-t paths). Use contrast −β_1 to
        // produce `mean_a − mean_b`. TREAT p-values use |t|, so the sign
        // flip doesn't affect p-values or f-statistics.
        // When external covariates are present, the design has p > 2
        // columns; the contrast for the non-group columns is 0.
        let n_design_cols = design.first().map(|r| r.len()).unwrap_or(2);
        let mut contrast_col = vec![0.0_f64; n_design_cols];
        if n_design_cols >= 2 {
            contrast_col[1] = -1.0; // negate group-B coefficient
        }
        let contrast_matrix: Vec<Vec<f64>> = contrast_col.into_iter().map(|v| vec![v]).collect();

        // (iv) Limma fit. When --peptide-metadata is supplied on the
        // limma path, count peptides per feature and switch the trend
        // covariate from mean-log2-abundance to `log(peptide_count +
        // 1)` (DEqMS, Zhu et al. 2020). The peptide→parent mapping
        // comes from the ingested `peptides.tsv`; features without
        // peptide records get count 0.
        let peptide_counts: Option<Vec<u32>> = if let Some(counts_map) = peptides_per_assay.as_ref()
        {
            let mut per_feature = Vec::with_capacity(features.len());
            for (_, _, assay_id, _) in &features {
                let c = counts_map.get(assay_id).copied().unwrap_or(0);
                per_feature.push(c);
            }
            Some(per_feature)
        } else {
            None
        };
        let deqms_enabled = peptide_counts.is_some();
        let options = LimmaOptions {
            trend: args.trend || deqms_enabled,
            robust: args.robust,
            winsor_lower: args.limma_winsor_lower,
            winsor_upper: args.limma_winsor_upper,
            lfc_threshold: args.lfc_threshold,
            peptide_counts,
        };
        let output = match limma_fit(&design, &y, &contrast_matrix, options) {
            Some(o) => o,
            None => {
                // Entire comparison degenerate (singular design, no
                // fittable feature, or df_residual <= 0). Emit skipped
                // result rows for each feature plus a report row per
                // panel so the output shape is consistent with other
                // dispatches.
                let mut per_panel: BTreeMap<String, usize> = BTreeMap::new();
                for (panel, gene, assay_id, uniprot) in &features {
                    *per_panel.entry(panel.clone()).or_insert(0) += 1;
                    all_result_rows.push(DeResultRow {
                        panel: panel.clone(),
                        assay_id: assay_id.clone(),
                        gene_symbol: gene.clone(),
                        uniprot: uniprot.join(","),
                        comparison: comparison_label.clone(),
                        n_pairs: design.len(),
                        mean_a: None,
                        mean_b: None,
                        mean_diff: None,
                        t: None,
                        df: None,
                        p_value: None,
                        bh_q: None,
                        effect_size: None,
                        effect_size_method: effect_size_method.to_string(),
                        ci_low: None,
                        ci_high: None,
                        wilcoxon_p: None,
                        wilcoxon_method: String::new(),
                        median_diff: None,
                        trimmed_mean_diff: None,
                        skip_reason: "limma: degenerate-design".to_string(),
                        s2_trend: None,
                        s2_prior: None,
                        s2_posterior: None,
                        df_prior: None,
                        df_total: None,
                        f_statistic: None,
                        f_p_value: None,
                        f_bh_q: None,
                        lfc_threshold: Some(args.lfc_threshold),
                        n_peptides_observed: peptides_per_assay
                            .as_ref()
                            .and_then(|m| m.get(assay_id).copied().map(|c| c as usize)),
                        peptide_variance_ratio: None,
                        ridge_lambda: None,
                        method: "limma".into(),
                        posthoc_method: String::new(),
                        posthoc_p: None,
                        posthoc_adj_p: None,
                    });
                }
                for (panel, n) in per_panel {
                    all_report_rows.push(DeReportRow {
                        comparison: comparison_label.clone(),
                        panel,
                        n_tests: 0,
                        n_skipped: n,
                        n_q_strict: 0,
                        n_q_relaxed: 0,
                        min_q: None,
                        max_abs_effect: None,
                        limma_trend_fallback_used: None,
                    });
                }
                continue;
            }
        };

        // Precompute raw per-comparison group means (over retained
        // samples) per feature for the `mean_a` / `mean_b` columns. These
        // are over the same samples that went into the limma fit.
        let mut group_a_idx: Vec<usize> = Vec::new();
        let mut group_b_idx: Vec<usize> = Vec::new();
        for (i, row) in design.iter().enumerate() {
            if row[1] > 0.5 {
                group_b_idx.push(i);
            } else {
                group_a_idx.push(i);
            }
        }

        // (v) BH within (comparison, panel). Build per-panel p-value
        // vectors plus back-index so we can splice adjusted q-values
        // back onto each feature's row.
        let n_features = features.len();
        let panel_by_feature: Vec<String> = features.iter().map(|(p, _, _, _)| p.clone()).collect();
        let mut p_by_panel: BTreeMap<String, Vec<Option<f64>>> = BTreeMap::new();
        let mut panel_positions: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (i, panel) in panel_by_feature.iter().enumerate() {
            let row = &output.rows[i];
            let p = if row.skipped || !row.p_values[0].is_finite() {
                None
            } else {
                Some(row.p_values[0])
            };
            p_by_panel.entry(panel.clone()).or_default().push(p);
            panel_positions.entry(panel.clone()).or_default().push(i);
        }
        let mut q_values: Vec<Option<f64>> = vec![None; n_features];
        for (panel, ps) in &p_by_panel {
            let qs = bh_fdr(ps);
            for (j, pos) in panel_positions[panel].iter().enumerate() {
                q_values[*pos] = qs[j];
            }
        }

        // (vi) Map LimmaRow → DeResultRow and accumulate per-panel
        // summary state for the report.
        // When df_total = Inf (prior fully pooled), use the z-score 1.96 as
        // the critical value (normal approximation).
        let t_crit: Option<f64> = if output.df_total.is_infinite() {
            Some(1.959963985)  // qnorm(0.975)
        } else {
            StudentsT::new(0.0, 1.0, output.df_total)
                .ok()
                .map(|d| d.inverse_cdf(0.975))
        };

        #[derive(Default)]
        struct PanelAcc {
            n_tests: usize,
            n_skipped: usize,
            n_q_strict: usize,
            n_q_relaxed: usize,
            min_q: Option<f64>,
            max_abs_effect: Option<f64>,
        }
        let mut panel_accs: BTreeMap<String, PanelAcc> = BTreeMap::new();

        for (i, (panel, gene, assay_id, uniprot)) in features.iter().enumerate() {
            let row = &output.rows[i];
            let y_row = &y[i];
            let mean_a = if group_a_idx.is_empty() {
                None
            } else {
                let vals: Vec<f64> = group_a_idx
                    .iter()
                    .map(|k| y_row[*k])
                    .filter(|v| v.is_finite())
                    .collect();
                if vals.is_empty() {
                    None
                } else {
                    Some(vals.iter().sum::<f64>() / vals.len() as f64)
                }
            };
            let mean_b = if group_b_idx.is_empty() {
                None
            } else {
                let vals: Vec<f64> = group_b_idx
                    .iter()
                    .map(|k| y_row[*k])
                    .filter(|v| v.is_finite())
                    .collect();
                if vals.is_empty() {
                    None
                } else {
                    Some(vals.iter().sum::<f64>() / vals.len() as f64)
                }
            };

            let (effect, ci_low, ci_high) = if row.skipped {
                (None, None, None)
            } else {
                let e = row.effects[0];
                let s = row.ses[0];
                let (cl, ch) = match t_crit {
                    Some(c) if s.is_finite() && c.is_finite() => (Some(e - c * s), Some(e + c * s)),
                    _ => (None, None),
                };
                (Some(e), cl, ch)
            };

            let p_value = if row.skipped || !row.p_values[0].is_finite() {
                None
            } else {
                Some(row.p_values[0])
            };
            let t = if row.skipped || !row.t_stats[0].is_finite() {
                None
            } else {
                Some(row.t_stats[0])
            };
            let q = q_values[i];

            let acc = panel_accs.entry(panel.clone()).or_default();
            if p_value.is_some() {
                acc.n_tests += 1;
                if let Some(qv) = q {
                    if qv < args.report_q_strict {
                        acc.n_q_strict += 1;
                    }
                    if qv < args.report_q_relaxed {
                        acc.n_q_relaxed += 1;
                    }
                    if acc.min_q.map(|m| qv < m).unwrap_or(true) {
                        acc.min_q = Some(qv);
                    }
                }
            } else {
                acc.n_skipped += 1;
            }
            if let Some(e) = effect {
                let abs_e = e.abs();
                if acc.max_abs_effect.map(|m| abs_e > m).unwrap_or(true) {
                    acc.max_abs_effect = Some(abs_e);
                }
            }

            all_result_rows.push(DeResultRow {
                panel: panel.clone(),
                assay_id: assay_id.clone(),
                gene_symbol: gene.clone(),
                uniprot: uniprot.join(","),
                comparison: comparison_label.clone(),
                n_pairs: design.len(),
                mean_a,
                mean_b,
                // Per plan Step 5e: `mean_diff: effect` — for the two-group
                // limma design (column 1 = indicator-for-b), the estimated
                // effect equals `mean_b - mean_a`. This is the opposite sign
                // of the paired-t / welch-t / ols dispatches, which report
                // `mean_a - mean_b`. The plan is explicit; downstream
                // consumers must treat the sign convention as test-specific.
                mean_diff: effect,
                t,
                df: Some(output.df_total),
                p_value,
                bh_q: q,
                effect_size: effect,
                effect_size_method: effect_size_method.to_string(),
                ci_low,
                ci_high,
                wilcoxon_p: None,
                wilcoxon_method: String::new(),
                median_diff: None,
                trimmed_mean_diff: None,
                skip_reason: if row.skipped {
                    "limma: ols-skipped".to_string()
                } else {
                    String::new()
                },
                s2_trend: row.s2_trend,
                s2_prior: if output.s2_prior.is_finite() {
                    Some(output.s2_prior)
                } else {
                    None
                },
                s2_posterior: if row.skipped || !row.s2_posterior.is_finite() {
                    None
                } else {
                    Some(row.s2_posterior)
                },
                df_prior: if output.df_prior.is_finite() {
                    Some(output.df_prior)
                } else {
                    None
                },
                df_total: Some(output.df_total),
                f_statistic: row.f_statistic,
                f_p_value: row.f_p_value,
                f_bh_q: None,
                lfc_threshold: Some(args.lfc_threshold),
                n_peptides_observed: peptides_per_assay
                    .as_ref()
                    .and_then(|m| m.get(assay_id).copied().map(|c| c as usize)),
                peptide_variance_ratio: None,
                ridge_lambda: None,
                method: "limma".into(),
                posthoc_method: String::new(),
                posthoc_p: None,
                posthoc_adj_p: None,
            });
        }

        // (vii) One DeReportRow per (comparison, panel).
        for (panel, acc) in panel_accs {
            all_report_rows.push(DeReportRow {
                comparison: comparison_label.clone(),
                panel,
                n_tests: acc.n_tests,
                n_skipped: acc.n_skipped,
                n_q_strict: acc.n_q_strict,
                n_q_relaxed: acc.n_q_relaxed,
                min_q: acc.min_q,
                max_abs_effect: acc.max_abs_effect,
                limma_trend_fallback_used: Some(output.trend_fallback_used),
            });
        }
    }

    Ok((all_result_rows, all_report_rows))
}
