//! Multiple-comparison-corrected post-hoc tests for multi-level
//! OLS designs: Sidak (step-down single-step), Tukey HSD (via
//! from-scratch studentized range), and Dunnett (via equicorrelated
//! multivariate-t; unbalanced uses deterministic MC on the observed
//! correlation matrix). All three build the complete-case full design
//! once, then emit one `de_results.tsv`-like table per comparison set
//! plus the design/covariate/omnibus sidecars.

use anyhow::{Context, Result};
use atman_core::de::SkipReason;
use atman_core::Sample;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;
use std::time::SystemTime;

use super::{
    classify_covariates, covariate_names_from_terms, parse_design_terms,
    push_encoded_covariate, read_covariate_columns, Args, CovKind, DesignReportRow, DesignTerm,
    OlsSetup,
};
use crate::io::{
    atomic_write, hash_canonical_inputs, read_measurements_long, read_proteins, read_samples,
    sidecar_path_for, write_de_results, write_run_sidecar, DeResultRow,
};

fn parse_contrast_list(list: &str) -> Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    for raw in list.split(',') {
        let token = raw.trim();
        if token.is_empty() {
            anyhow::bail!("invalid empty contrast in --contrast-list {:?}", list);
        }
        let mut parts = token.split('-').map(str::trim);
        let a = parts.next().unwrap_or_default();
        let b = parts.next().unwrap_or_default();
        if a.is_empty() || b.is_empty() || parts.next().is_some() || a == b {
            anyhow::bail!(
                "invalid contrast {:?}; expected `<level_a>-<level_b>` with \
                 two distinct non-empty levels",
                token
            );
        }
        out.push((a.to_string(), b.to_string()));
    }
    if out.is_empty() {
        anyhow::bail!("--contrast-list is empty");
    }
    Ok(out)
}

/// Dispatch for `--post-hoc sidak`: one OLS fit per protein, then
/// one contrast per entry in `--contrast-list`, Sidak-adjusted
/// within the list. Emits one `DeResultRow` per (protein, contrast).
///
/// Samples are the full set of non-control samples with finite
/// values for every fixed covariate in `--design`; no `--groups`
/// subsetting. The factor for contrast evaluation comes from
/// `--post-hoc-factor` or `--omnibus-factor` (whichever is set),
/// must be categorical, and must have all levels referenced by
/// the contrast list in its observed levels.
pub(super) fn run_posthoc_sidak(args: Args, started_at: SystemTime) -> Result<()> {
    use atman_core::de::{contrast_inference, ols, OlsOutcome};
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let contrast_list = parse_contrast_list(
        args.contrast_list
            .as_deref()
            .expect("validated non-empty above"),
    )?;
    let factor_name = args
        .post_hoc_factor
        .clone()
        .or_else(|| args.omnibus_factor.clone())
        .expect("validated non-empty above");

    // Read inputs.
    let measurements = read_measurements_long(&args.input_dir.join("qc_measurements.tsv"))?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    let proteins = read_proteins(&args.input_dir.join("proteins.tsv"))?;

    // Build the covariate setup from --design. The standard
    // `build_ols_setup` requires `condition` in the formula (to
    // anchor the binary contrast); post-hoc doesn't use that
    // machinery, so we parse the formula directly and skip the
    // condition check.
    let formula = args.design.as_deref().expect("validated above");
    let terms = parse_design_terms(formula)?;
    if terms.iter().any(|t| matches!(t, DesignTerm::Condition)) {
        anyhow::bail!(
            "--post-hoc sidak does not support `condition` in --design; \
             include the factor {:?} directly (e.g. \"~ {} + age + sex\")",
            factor_name,
            factor_name
        );
    }
    let cov_names = covariate_names_from_terms(&terms);
    let cov_raw = if cov_names.is_empty() {
        HashMap::new()
    } else {
        read_covariate_columns(&args.input_dir.join("samples.tsv"), &cov_names)?
    };
    let setup = OlsSetup {
        terms,
        cov_names,
        cov_raw,
        contrast: None,
        label: formula.to_string(),
    };
    if !setup
        .cov_names
        .iter()
        .any(|n| n == &factor_name)
    {
        anyhow::bail!(
            "--post-hoc-factor {:?} is not in --design covariates {:?}",
            factor_name,
            setup.cov_names
        );
    }

    // Build the full design: every non-control sample with finite
    // values for every covariate in `setup`.
    let PosthocDesign {
        rows: design_rows,
        labels: design_labels,
        cov_kinds,
    } = build_posthoc_full_design(&samples, &setup, &factor_name)?;
    if design_rows.is_empty() {
        anyhow::bail!(
            "no samples retained after complete-case filter on --design covariates"
        );
    }

    // Resolve factor-column-span + ref level for contrast-vector
    // construction.
    let factor_cov_idx = setup
        .cov_names
        .iter()
        .position(|n| n == &factor_name)
        .expect("validated present");
    let factor_levels: Vec<String> = match &cov_kinds[factor_cov_idx] {
        CovKind::Categorical { levels } => levels.clone(),
        _ => anyhow::bail!(
            "--post-hoc-factor {:?} must be categorical; looks numeric in the data",
            factor_name
        ),
    };
    if factor_levels.len() < 2 {
        anyhow::bail!(
            "factor {:?} has only {} observed level(s); need at least 2",
            factor_name,
            factor_levels.len()
        );
    }
    for (a, b) in &contrast_list {
        if !factor_levels.contains(a) {
            anyhow::bail!(
                "contrast level {:?} not among observed factor levels {:?}",
                a,
                factor_levels
            );
        }
        if !factor_levels.contains(b) {
            anyhow::bail!(
                "contrast level {:?} not among observed factor levels {:?}",
                b,
                factor_levels
            );
        }
    }
    // Reference level = alphabetically-first (matches CovKind::col_labels).
    let ref_level = factor_levels[0].clone();
    // Map each non-reference level to its column index in the design.
    let factor_col_by_level: HashMap<String, usize> = factor_levels
        .iter()
        .skip(1)
        .map(|lvl| {
            let label = format!("{factor_name}{lvl}");
            let idx = design_labels
                .iter()
                .position(|l| l == &label)
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "factor column {:?} missing from design labels {:?}",
                        label,
                        design_labels
                    )
                })?;
            Ok::<_, anyhow::Error>((lvl.clone(), idx))
        })
        .collect::<Result<_>>()?;

    // Build a sample→abundance lookup per (panel, gene).
    let sample_by_id: HashMap<&str, &Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();
    let mut gene_meta: BTreeMap<(String, String), (String, Vec<String>)> = BTreeMap::new();
    for p in &proteins {
        if let (Some(gene), Some(panel)) = (p.gene_symbol.as_ref(), p.panel.as_ref()) {
            gene_meta
                .entry((panel.clone(), gene.clone()))
                .or_insert_with(|| (p.assay_id.0.clone(), p.uniprot.clone()));
        }
    }
    let mut cells_by_sample: HashMap<(String, String, String), f64> = HashMap::new();
    let mut measured_features: BTreeSet<(String, String)> = BTreeSet::new();
    for m in &measurements {
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

    let m_contrasts = contrast_list.len();
    let mut all_rows: Vec<DeResultRow> = Vec::new();
    for (panel, gene) in &measured_features {
        let (assay_id, uniprot) = gene_meta
            .get(&(panel.clone(), gene.clone()))
            .cloned()
            .unwrap_or_else(|| (String::new(), vec![]));

        // Complete-case per protein: subset design + y to samples
        // with finite abundance for this (panel, gene).
        let mut design_cc: Vec<Vec<f64>> = Vec::new();
        let mut y_cc: Vec<f64> = Vec::new();
        for (sid, row) in &design_rows {
            if let Some(&abund) =
                cells_by_sample.get(&(panel.clone(), gene.clone(), sid.clone()))
            {
                if abund.is_finite() {
                    design_cc.push(row.clone());
                    y_cc.push(abund);
                }
            }
        }
        let n = y_cc.len();
        let contrast_label_for = |(a, b): &(String, String)| format!("{a}-{b}");
        match ols(&design_cc, &y_cc, args.min_pairs) {
            OlsOutcome::Computed(fit) => {
                let mut raw_ps: Vec<f64> = Vec::with_capacity(m_contrasts);
                let mut ests: Vec<f64> = Vec::with_capacity(m_contrasts);
                let mut ses: Vec<f64> = Vec::with_capacity(m_contrasts);
                let mut ts: Vec<f64> = Vec::with_capacity(m_contrasts);
                for (a, b) in &contrast_list {
                    let c = build_posthoc_contrast_vector(
                        a,
                        b,
                        &ref_level,
                        &factor_col_by_level,
                        design_labels.len(),
                    );
                    match contrast_inference(&design_cc, &fit.beta, &c, fit.sigma2, fit.df) {
                        Some(r) => {
                            raw_ps.push(r.p_value);
                            ests.push(r.estimate);
                            ses.push(r.se);
                            ts.push(r.t);
                        }
                        None => {
                            raw_ps.push(f64::NAN);
                            ests.push(f64::NAN);
                            ses.push(f64::NAN);
                            ts.push(f64::NAN);
                        }
                    }
                }
                let m_f = m_contrasts as f64;
                for (i, contrast_pair) in contrast_list.iter().enumerate() {
                    let p = raw_ps[i];
                    let adj = if p.is_finite() {
                        1.0 - (1.0 - p.clamp(0.0, 1.0)).powf(m_f)
                    } else {
                        f64::NAN
                    };
                    all_rows.push(DeResultRow {
                        panel: panel.clone(),
                        assay_id: assay_id.clone(),
                        gene_symbol: gene.clone(),
                        uniprot: uniprot.join(","),
                        comparison: contrast_label_for(contrast_pair),
                        n_pairs: n,
                        mean_a: None,
                        mean_b: None,
                        mean_diff: if ests[i].is_finite() { Some(ests[i]) } else { None },
                        t: if ts[i].is_finite() { Some(ts[i]) } else { None },
                        df: if fit.df.is_finite() { Some(fit.df) } else { None },
                        p_value: if raw_ps[i].is_finite() { Some(raw_ps[i]) } else { None },
                        bh_q: None,
                        effect_size: if ests[i].is_finite() { Some(ests[i]) } else { None },
                        effect_size_method: "ols-posthoc-sidak".into(),
                        ci_low: if ests[i].is_finite() && ses[i].is_finite() {
                            Some(ests[i] - 1.96 * ses[i])
                        } else {
                            None
                        },
                        ci_high: if ests[i].is_finite() && ses[i].is_finite() {
                            Some(ests[i] + 1.96 * ses[i])
                        } else {
                            None
                        },
                        wilcoxon_p: None,
                        wilcoxon_method: String::new(),
                        median_diff: None,
                        trimmed_mean_diff: None,
                        skip_reason: String::new(),
                        s2_trend: None,
                        s2_prior: None,
                        s2_posterior: None,
                        df_prior: None,
                        df_total: None,
                        f_statistic: None,
                        f_p_value: None,
                        f_bh_q: None,
                        lfc_threshold: None,
                        n_peptides_observed: None,
                        peptide_variance_ratio: None,
                        ridge_lambda: None,
                        method: "ols".into(),
                        posthoc_method: "sidak".into(),
                        posthoc_p: if raw_ps[i].is_finite() { Some(raw_ps[i]) } else { None },
                        posthoc_adj_p: if adj.is_finite() { Some(adj) } else { None },
                    });
                }
            }
            OlsOutcome::Skipped { reason, n } => {
                let reason_str = match reason {
                    SkipReason::InsufficientPairs => "insufficient_samples",
                    SkipReason::ZeroVariance => "zero_variance",
                    SkipReason::NonFiniteInput => "non_finite_input",
                };
                for contrast_pair in &contrast_list {
                    all_rows.push(DeResultRow {
                        panel: panel.clone(),
                        assay_id: assay_id.clone(),
                        gene_symbol: gene.clone(),
                        uniprot: uniprot.join(","),
                        comparison: contrast_label_for(contrast_pair),
                        n_pairs: n,
                        mean_a: None,
                        mean_b: None,
                        mean_diff: None,
                        t: None,
                        df: None,
                        p_value: None,
                        bh_q: None,
                        effect_size: None,
                        effect_size_method: "ols-posthoc-sidak".into(),
                        ci_low: None,
                        ci_high: None,
                        wilcoxon_p: None,
                        wilcoxon_method: String::new(),
                        median_diff: None,
                        trimmed_mean_diff: None,
                        skip_reason: reason_str.into(),
                        s2_trend: None,
                        s2_prior: None,
                        s2_posterior: None,
                        df_prior: None,
                        df_total: None,
                        f_statistic: None,
                        f_p_value: None,
                        f_bh_q: None,
                        lfc_threshold: None,
                        n_peptides_observed: None,
                        peptide_variance_ratio: None,
                        ridge_lambda: None,
                        method: "ols".into(),
                        posthoc_method: "sidak".into(),
                        posthoc_p: None,
                        posthoc_adj_p: None,
                    });
                }
            }
        }
    }

    // BH-q per (comparison, panel) across proteins in a contrast family.
    let mut by_family: BTreeMap<(String, String), Vec<usize>> = BTreeMap::new();
    for (i, r) in all_rows.iter().enumerate() {
        by_family
            .entry((r.comparison.clone(), r.panel.clone()))
            .or_default()
            .push(i);
    }
    for indices in by_family.values() {
        let ps: Vec<Option<f64>> = indices.iter().map(|&i| all_rows[i].p_value).collect();
        let qs = atman_core::bh_fdr(&ps);
        for (j, &i) in indices.iter().enumerate() {
            all_rows[i].bh_q = qs[j];
        }
    }

    // Sort for deterministic output.
    all_rows.sort_by(|a, b| {
        a.comparison
            .cmp(&b.comparison)
            .then_with(|| a.panel.cmp(&b.panel))
            .then_with(|| match (a.posthoc_adj_p, b.posthoc_adj_p) {
                (Some(x), Some(y)) => x.partial_cmp(&y).unwrap_or(std::cmp::Ordering::Equal),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            })
            .then_with(|| a.gene_symbol.cmp(&b.gene_symbol))
    });

    let results_path = args.output_dir.join("de_results.tsv");
    write_de_results(&results_path, &all_rows)?;
    let outputs = vec![results_path.clone()];

    eprintln!(
        "de posthoc sidak: factor={:?} contrasts={} proteins={} rows={}",
        factor_name,
        m_contrasts,
        measured_features.len(),
        all_rows.len()
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &[
            "qc_measurements.tsv",
            "measurements.tsv",
            "samples.tsv",
            "proteins.tsv",
        ],
    )?;
    let sidecar = sidecar_path_for(&results_path);
    write_run_sidecar(
        &sidecar,
        "de",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output-dir": args.output_dir.display().to_string(),
            "test": "ols",
            "design": args.design,
            "groups": args.groups,
            "min-pairs": args.min_pairs,
            "post-hoc": args.post_hoc,
            "post-hoc-factor": factor_name,
            "contrast-list": args.contrast_list,
            "alpha": args.alpha,
            "omnibus-factor": args.omnibus_factor,
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        None,
)?;
    eprintln!("de posthoc sidak: sidecar={}", sidecar.display());
    Ok(())
}

/// Dispatch for `--post-hoc tukey`: one OLS fit per protein, then
/// every ordered pair `(level_a, level_b)` of the post-hoc factor,
/// with p-values adjusted via Tukey's studentized range. Compared
/// to sidak the only differences are (a) the contrast list defaults
/// to every ordered pair of factor levels and (b) the adjustment
/// uses `1 − ptukey(|estimate|·√2 / se, nmeans=k, df=residual)`
/// rather than `1 − (1 − p)^m`.
///
/// Why `|estimate|·√2 / se`: the covariate-adjusted contrast SE
/// from `contrast_inference` is `√(σ² · c'·(X'X)⁻¹·c)`; the
/// Tukey–Kramer SE is that divided by `√2`, so the studentized
/// range statistic `q = diff / SE_tukey = diff · √2 / SE_contrast`.
/// This matches `emmeans(..., adjust = "tukey")` under arbitrary
/// covariate adjustment and unbalanced `n_i`.
pub(super) fn run_posthoc_tukey(args: Args, started_at: SystemTime) -> Result<()> {
    use atman_core::de::{contrast_inference, ols, OlsOutcome};
    use atman_core::studentized_range::ptukey;
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let factor_name = args
        .post_hoc_factor
        .clone()
        .or_else(|| args.omnibus_factor.clone())
        .expect("validated non-empty above");

    let measurements = read_measurements_long(&args.input_dir.join("qc_measurements.tsv"))?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    let proteins = read_proteins(&args.input_dir.join("proteins.tsv"))?;

    let formula = args.design.as_deref().expect("validated above");
    let terms = parse_design_terms(formula)?;
    if terms.iter().any(|t| matches!(t, DesignTerm::Condition)) {
        anyhow::bail!(
            "--post-hoc tukey does not support `condition` in --design; \
             include the factor {:?} directly (e.g. \"~ {} + age + sex\")",
            factor_name,
            factor_name
        );
    }
    let cov_names = covariate_names_from_terms(&terms);
    let cov_raw = if cov_names.is_empty() {
        HashMap::new()
    } else {
        read_covariate_columns(&args.input_dir.join("samples.tsv"), &cov_names)?
    };
    let setup = OlsSetup {
        terms,
        cov_names,
        cov_raw,
        contrast: None,
        label: formula.to_string(),
    };
    if !setup.cov_names.iter().any(|n| n == &factor_name) {
        anyhow::bail!(
            "--post-hoc-factor {:?} is not in --design covariates {:?}",
            factor_name,
            setup.cov_names
        );
    }

    let PosthocDesign {
        rows: design_rows,
        labels: design_labels,
        cov_kinds,
    } = build_posthoc_full_design(&samples, &setup, &factor_name)?;
    if design_rows.is_empty() {
        anyhow::bail!(
            "no samples retained after complete-case filter on --design covariates"
        );
    }

    let factor_cov_idx = setup
        .cov_names
        .iter()
        .position(|n| n == &factor_name)
        .expect("validated present");
    let factor_levels: Vec<String> = match &cov_kinds[factor_cov_idx] {
        CovKind::Categorical { levels } => levels.clone(),
        _ => anyhow::bail!(
            "--post-hoc-factor {:?} must be categorical; looks numeric in the data",
            factor_name
        ),
    };
    if factor_levels.len() < 2 {
        anyhow::bail!(
            "factor {:?} has only {} observed level(s); need at least 2",
            factor_name,
            factor_levels.len()
        );
    }

    // Resolve the contrast list: either user-supplied subset, or
    // every ordered pair of levels (the canonical Tukey HSD family).
    let contrast_list: Vec<(String, String)> = match args.contrast_list.as_deref() {
        Some(s) if !s.is_empty() => {
            let list = parse_contrast_list(s)?;
            for (a, b) in &list {
                if !factor_levels.contains(a) {
                    anyhow::bail!(
                        "contrast level {:?} not among observed factor levels {:?}",
                        a, factor_levels
                    );
                }
                if !factor_levels.contains(b) {
                    anyhow::bail!(
                        "contrast level {:?} not among observed factor levels {:?}",
                        b, factor_levels
                    );
                }
            }
            list
        }
        _ => {
            let mut out = Vec::new();
            for i in 0..factor_levels.len() {
                for j in (i + 1)..factor_levels.len() {
                    out.push((factor_levels[i].clone(), factor_levels[j].clone()));
                }
            }
            out
        }
    };

    let ref_level = factor_levels[0].clone();
    let factor_col_by_level: HashMap<String, usize> = factor_levels
        .iter()
        .skip(1)
        .map(|lvl| {
            let label = format!("{factor_name}{lvl}");
            let idx = design_labels
                .iter()
                .position(|l| l == &label)
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "factor column {:?} missing from design labels {:?}",
                        label, design_labels
                    )
                })?;
            Ok::<_, anyhow::Error>((lvl.clone(), idx))
        })
        .collect::<Result<_>>()?;

    // Abundance lookup by (panel, gene, sample).
    let sample_by_id: HashMap<&str, &Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();
    let mut gene_meta: BTreeMap<(String, String), (String, Vec<String>)> = BTreeMap::new();
    for p in &proteins {
        if let (Some(gene), Some(panel)) = (p.gene_symbol.as_ref(), p.panel.as_ref()) {
            gene_meta
                .entry((panel.clone(), gene.clone()))
                .or_insert_with(|| (p.assay_id.0.clone(), p.uniprot.clone()));
        }
    }
    let mut cells_by_sample: HashMap<(String, String, String), f64> = HashMap::new();
    let mut measured_features: BTreeSet<(String, String)> = BTreeSet::new();
    for m in &measurements {
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

    let k_nmeans = factor_levels.len();
    let mut all_rows: Vec<DeResultRow> = Vec::new();
    for (panel, gene) in &measured_features {
        let (assay_id, uniprot) = gene_meta
            .get(&(panel.clone(), gene.clone()))
            .cloned()
            .unwrap_or_else(|| (String::new(), vec![]));

        let mut design_cc: Vec<Vec<f64>> = Vec::new();
        let mut y_cc: Vec<f64> = Vec::new();
        for (sid, row) in &design_rows {
            if let Some(&abund) =
                cells_by_sample.get(&(panel.clone(), gene.clone(), sid.clone()))
            {
                if abund.is_finite() {
                    design_cc.push(row.clone());
                    y_cc.push(abund);
                }
            }
        }
        let n = y_cc.len();
        let contrast_label_for = |(a, b): &(String, String)| format!("{a}-{b}");
        match ols(&design_cc, &y_cc, args.min_pairs) {
            OlsOutcome::Computed(fit) => {
                for (a, b) in &contrast_list {
                    let c = build_posthoc_contrast_vector(
                        a, b, &ref_level, &factor_col_by_level, design_labels.len(),
                    );
                    let (raw_p, est, se, t_stat, adj_p) = match contrast_inference(
                        &design_cc, &fit.beta, &c, fit.sigma2, fit.df,
                    ) {
                        Some(r) => {
                            let q = if r.se > 0.0 && r.se.is_finite() {
                                r.estimate.abs() * std::f64::consts::SQRT_2 / r.se
                            } else {
                                f64::NAN
                            };
                            let p_tukey = if q.is_finite() && fit.df.is_finite() {
                                (1.0 - ptukey(q, k_nmeans, fit.df)).clamp(0.0, 1.0)
                            } else {
                                f64::NAN
                            };
                            (r.p_value, r.estimate, r.se, r.t, p_tukey)
                        }
                        None => (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN),
                    };
                    all_rows.push(DeResultRow {
                        panel: panel.clone(),
                        assay_id: assay_id.clone(),
                        gene_symbol: gene.clone(),
                        uniprot: uniprot.join(","),
                        comparison: contrast_label_for(&(a.clone(), b.clone())),
                        n_pairs: n,
                        mean_a: None,
                        mean_b: None,
                        mean_diff: if est.is_finite() { Some(est) } else { None },
                        t: if t_stat.is_finite() { Some(t_stat) } else { None },
                        df: if fit.df.is_finite() { Some(fit.df) } else { None },
                        p_value: if raw_p.is_finite() { Some(raw_p) } else { None },
                        bh_q: None,
                        effect_size: if est.is_finite() { Some(est) } else { None },
                        effect_size_method: "ols-posthoc-tukey".into(),
                        ci_low: if est.is_finite() && se.is_finite() {
                            Some(est - 1.96 * se)
                        } else {
                            None
                        },
                        ci_high: if est.is_finite() && se.is_finite() {
                            Some(est + 1.96 * se)
                        } else {
                            None
                        },
                        wilcoxon_p: None,
                        wilcoxon_method: String::new(),
                        median_diff: None,
                        trimmed_mean_diff: None,
                        skip_reason: String::new(),
                        s2_trend: None,
                        s2_prior: None,
                        s2_posterior: None,
                        df_prior: None,
                        df_total: None,
                        f_statistic: None,
                        f_p_value: None,
                        f_bh_q: None,
                        lfc_threshold: None,
                        n_peptides_observed: None,
                        peptide_variance_ratio: None,
                        ridge_lambda: None,
                        method: "ols".into(),
                        posthoc_method: "tukey".into(),
                        posthoc_p: if raw_p.is_finite() { Some(raw_p) } else { None },
                        posthoc_adj_p: if adj_p.is_finite() { Some(adj_p) } else { None },
                    });
                }
            }
            OlsOutcome::Skipped { reason, n } => {
                let reason_str = match reason {
                    SkipReason::InsufficientPairs => "insufficient_samples",
                    SkipReason::ZeroVariance => "zero_variance",
                    SkipReason::NonFiniteInput => "non_finite_input",
                };
                for (a, b) in &contrast_list {
                    all_rows.push(DeResultRow {
                        panel: panel.clone(),
                        assay_id: assay_id.clone(),
                        gene_symbol: gene.clone(),
                        uniprot: uniprot.join(","),
                        comparison: contrast_label_for(&(a.clone(), b.clone())),
                        n_pairs: n,
                        mean_a: None, mean_b: None, mean_diff: None,
                        t: None, df: None, p_value: None, bh_q: None,
                        effect_size: None,
                        effect_size_method: "ols-posthoc-tukey".into(),
                        ci_low: None, ci_high: None,
                        wilcoxon_p: None, wilcoxon_method: String::new(),
                        median_diff: None, trimmed_mean_diff: None,
                        skip_reason: reason_str.into(),
                        s2_trend: None, s2_prior: None, s2_posterior: None,
                        df_prior: None, df_total: None,
                        f_statistic: None, f_p_value: None, f_bh_q: None,
                        lfc_threshold: None,
                        n_peptides_observed: None, peptide_variance_ratio: None,
                        ridge_lambda: None,
                        method: "ols".into(),
                        posthoc_method: "tukey".into(),
                        posthoc_p: None, posthoc_adj_p: None,
                    });
                }
            }
        }
    }

    // BH-q per (comparison, panel) across proteins in the family.
    let mut by_family: BTreeMap<(String, String), Vec<usize>> = BTreeMap::new();
    for (i, r) in all_rows.iter().enumerate() {
        by_family
            .entry((r.comparison.clone(), r.panel.clone()))
            .or_default()
            .push(i);
    }
    for indices in by_family.values() {
        let ps: Vec<Option<f64>> = indices.iter().map(|&i| all_rows[i].p_value).collect();
        let qs = atman_core::bh_fdr(&ps);
        for (j, &i) in indices.iter().enumerate() {
            all_rows[i].bh_q = qs[j];
        }
    }
    all_rows.sort_by(|a, b| {
        a.comparison
            .cmp(&b.comparison)
            .then_with(|| a.panel.cmp(&b.panel))
            .then_with(|| match (a.posthoc_adj_p, b.posthoc_adj_p) {
                (Some(x), Some(y)) => x.partial_cmp(&y).unwrap_or(std::cmp::Ordering::Equal),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            })
            .then_with(|| a.gene_symbol.cmp(&b.gene_symbol))
    });

    let results_path = args.output_dir.join("de_results.tsv");
    write_de_results(&results_path, &all_rows)?;
    let outputs = vec![results_path.clone()];
    eprintln!(
        "de posthoc tukey: factor={:?} pairs={} proteins={} rows={}",
        factor_name,
        contrast_list.len(),
        measured_features.len(),
        all_rows.len()
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &[
            "qc_measurements.tsv",
            "measurements.tsv",
            "samples.tsv",
            "proteins.tsv",
        ],
    )?;
    let sidecar = sidecar_path_for(&results_path);
    write_run_sidecar(
        &sidecar,
        "de",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output-dir": args.output_dir.display().to_string(),
            "test": "ols",
            "design": args.design,
            "groups": args.groups,
            "min-pairs": args.min_pairs,
            "post-hoc": args.post_hoc,
            "post-hoc-factor": factor_name,
            "contrast-list": args.contrast_list,
            "alpha": args.alpha,
            "omnibus-factor": args.omnibus_factor,
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        None,
)?;
    eprintln!("de posthoc tukey: sidecar={}", sidecar.display());
    Ok(())
}

/// Dispatch for `--post-hoc dunnett`: one OLS fit per protein, then
/// every `(non_control, control)` pair where `control` is the
/// alphabetically-first observed level of the post-hoc factor.
///
/// Adjustment: balanced designs (`max n_i / min n_i < 1.25`) use the
/// equicorrelated multivariate-t CDF via
/// `atman_core::multivariate_t::pdunnett` at `ρ = 0.5`; unbalanced
/// designs use the Hsu variant with a full correlation matrix
/// computed from per-level `n_i` and evaluated by deterministic
/// Monte Carlo (`atman_core::multivariate_t::pdunnett_hsu` at 50 000
/// draws, byte-equal under fixed `--seed`).
pub(super) fn run_posthoc_dunnett(args: Args, started_at: SystemTime) -> Result<()> {
    use atman_core::de::{contrast_inference, ols, OlsOutcome};
    use atman_core::multivariate_t::{
        dunnett_hsu_correlation_matrix, pdunnett, pdunnett_hsu,
    };
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    let factor_name = args
        .post_hoc_factor
        .clone()
        .or_else(|| args.omnibus_factor.clone())
        .expect("validated non-empty above");

    let measurements = read_measurements_long(&args.input_dir.join("qc_measurements.tsv"))?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    let proteins = read_proteins(&args.input_dir.join("proteins.tsv"))?;

    let formula = args.design.as_deref().expect("validated above");
    let terms = parse_design_terms(formula)?;
    if terms.iter().any(|t| matches!(t, DesignTerm::Condition)) {
        anyhow::bail!(
            "--post-hoc dunnett does not support `condition` in --design; \
             include the factor {:?} directly",
            factor_name
        );
    }
    let cov_names = covariate_names_from_terms(&terms);
    let cov_raw = if cov_names.is_empty() {
        HashMap::new()
    } else {
        read_covariate_columns(&args.input_dir.join("samples.tsv"), &cov_names)?
    };
    let setup = OlsSetup {
        terms,
        cov_names,
        cov_raw,
        contrast: None,
        label: formula.to_string(),
    };
    if !setup.cov_names.iter().any(|n| n == &factor_name) {
        anyhow::bail!(
            "--post-hoc-factor {:?} is not in --design covariates {:?}",
            factor_name,
            setup.cov_names
        );
    }

    let PosthocDesign {
        rows: design_rows,
        labels: design_labels,
        cov_kinds,
    } = build_posthoc_full_design(&samples, &setup, &factor_name)?;
    if design_rows.is_empty() {
        anyhow::bail!(
            "no samples retained after complete-case filter on --design covariates"
        );
    }

    let factor_cov_idx = setup
        .cov_names
        .iter()
        .position(|n| n == &factor_name)
        .expect("validated present");
    let factor_levels: Vec<String> = match &cov_kinds[factor_cov_idx] {
        CovKind::Categorical { levels } => levels.clone(),
        _ => anyhow::bail!(
            "--post-hoc-factor {:?} must be categorical; looks numeric in the data",
            factor_name
        ),
    };
    if factor_levels.len() < 2 {
        anyhow::bail!(
            "factor {:?} has only {} observed level(s); need at least 2",
            factor_name,
            factor_levels.len()
        );
    }

    // Control is the reference level (alphabetically first observed).
    // Contrasts are (treat, control) for each non-reference level.
    let control = factor_levels[0].clone();
    let contrast_list: Vec<(String, String)> = factor_levels
        .iter()
        .skip(1)
        .map(|lvl| (lvl.clone(), control.clone()))
        .collect();
    let m_contrasts = contrast_list.len();

    let factor_col_by_level: HashMap<String, usize> = factor_levels
        .iter()
        .skip(1)
        .map(|lvl| {
            let label = format!("{factor_name}{lvl}");
            let idx = design_labels
                .iter()
                .position(|l| l == &label)
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "factor column {:?} missing from design labels {:?}",
                        label, design_labels
                    )
                })?;
            Ok::<_, anyhow::Error>((lvl.clone(), idx))
        })
        .collect::<Result<_>>()?;

    let sample_by_id: HashMap<&str, &Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();
    let mut gene_meta: BTreeMap<(String, String), (String, Vec<String>)> = BTreeMap::new();
    for p in &proteins {
        if let (Some(gene), Some(panel)) = (p.gene_symbol.as_ref(), p.panel.as_ref()) {
            gene_meta
                .entry((panel.clone(), gene.clone()))
                .or_insert_with(|| (p.assay_id.0.clone(), p.uniprot.clone()));
        }
    }
    let mut cells_by_sample: HashMap<(String, String, String), f64> = HashMap::new();
    let mut measured_features: BTreeSet<(String, String)> = BTreeSet::new();
    for m in &measurements {
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

    const RHO: f64 = 0.5;
    const HSU_N_MC: usize = 50_000;
    const BALANCE_THRESHOLD: f64 = 1.25;
    // Count n_i per level from the full design. Hsu's correlation
    // matrix depends only on these counts.
    let mut level_counts: HashMap<String, usize> = HashMap::new();
    for (sid, _) in &design_rows {
        // Resolve the sample's factor level via the cov_raw lookup.
        let _ = sid;
    }
    // The simpler way: iterate samples and read the factor value via
    // setup.cov_raw keyed by sample_id.
    for sample in &samples {
        if sample.is_control {
            continue;
        }
        let per_cov = match setup.cov_raw.get(&sample.sample_id) {
            Some(v) => v,
            None => continue,
        };
        let val = match per_cov.get(factor_cov_idx).cloned().unwrap_or(None) {
            Some(v) if !v.is_empty() => v,
            _ => continue,
        };
        *level_counts.entry(val).or_insert(0) += 1;
    }
    let n_control = *level_counts.get(&control).unwrap_or(&0);
    let n_treatments: Vec<usize> = factor_levels
        .iter()
        .skip(1)
        .map(|lvl| *level_counts.get(lvl).unwrap_or(&0))
        .collect();
    let max_n = n_treatments.iter().copied().max().unwrap_or(1).max(n_control) as f64;
    let min_n = n_treatments
        .iter()
        .copied()
        .chain(std::iter::once(n_control))
        .filter(|n| *n > 0)
        .min()
        .unwrap_or(1) as f64;
    let is_unbalanced = max_n / min_n.max(1.0) > BALANCE_THRESHOLD;
    let hsu_matrix = if is_unbalanced {
        Some(dunnett_hsu_correlation_matrix(n_control, &n_treatments))
    } else {
        None
    };
    if is_unbalanced {
        eprintln!(
            "de posthoc dunnett: unbalanced design detected (n_i = {:?}, control {}) — \
             using Dunnett-Hsu via Monte Carlo ({} draws)",
            n_treatments, n_control, HSU_N_MC,
        );
    }
    let hsu_mc_seed = 20260420_u64; // fixed seed for byte-equal CDF calls
    let dunnett_cdf = |q: f64, df: f64| -> f64 {
        match (&hsu_matrix, q.is_finite() && df.is_finite()) {
            (_, false) => f64::NAN,
            (Some(r), _) => pdunnett_hsu(q, df, r, HSU_N_MC, hsu_mc_seed),
            (None, _) => pdunnett(q, m_contrasts, df, RHO),
        }
    };
    let mut all_rows: Vec<DeResultRow> = Vec::new();
    for (panel, gene) in &measured_features {
        let (assay_id, uniprot) = gene_meta
            .get(&(panel.clone(), gene.clone()))
            .cloned()
            .unwrap_or_else(|| (String::new(), vec![]));

        let mut design_cc: Vec<Vec<f64>> = Vec::new();
        let mut y_cc: Vec<f64> = Vec::new();
        for (sid, row) in &design_rows {
            if let Some(&abund) =
                cells_by_sample.get(&(panel.clone(), gene.clone(), sid.clone()))
            {
                if abund.is_finite() {
                    design_cc.push(row.clone());
                    y_cc.push(abund);
                }
            }
        }
        let n = y_cc.len();
        let contrast_label_for = |(a, b): &(String, String)| format!("{a}-{b}");
        match ols(&design_cc, &y_cc, args.min_pairs) {
            OlsOutcome::Computed(fit) => {
                for (a, b) in &contrast_list {
                    let c = build_posthoc_contrast_vector(
                        a, b, &control, &factor_col_by_level, design_labels.len(),
                    );
                    let (raw_p, est, se, t_stat, adj_p) = match contrast_inference(
                        &design_cc, &fit.beta, &c, fit.sigma2, fit.df,
                    ) {
                        Some(r) => {
                            let q = r.t.abs();
                            let p_dun = if q.is_finite() && fit.df.is_finite() {
                                (1.0 - dunnett_cdf(q, fit.df)).clamp(0.0, 1.0)
                            } else {
                                f64::NAN
                            };
                            (r.p_value, r.estimate, r.se, r.t, p_dun)
                        }
                        None => (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN),
                    };
                    all_rows.push(DeResultRow {
                        panel: panel.clone(),
                        assay_id: assay_id.clone(),
                        gene_symbol: gene.clone(),
                        uniprot: uniprot.join(","),
                        comparison: contrast_label_for(&(a.clone(), b.clone())),
                        n_pairs: n,
                        mean_a: None,
                        mean_b: None,
                        mean_diff: if est.is_finite() { Some(est) } else { None },
                        t: if t_stat.is_finite() { Some(t_stat) } else { None },
                        df: if fit.df.is_finite() { Some(fit.df) } else { None },
                        p_value: if raw_p.is_finite() { Some(raw_p) } else { None },
                        bh_q: None,
                        effect_size: if est.is_finite() { Some(est) } else { None },
                        effect_size_method: "ols-posthoc-dunnett".into(),
                        ci_low: if est.is_finite() && se.is_finite() {
                            Some(est - 1.96 * se)
                        } else {
                            None
                        },
                        ci_high: if est.is_finite() && se.is_finite() {
                            Some(est + 1.96 * se)
                        } else {
                            None
                        },
                        wilcoxon_p: None,
                        wilcoxon_method: String::new(),
                        median_diff: None,
                        trimmed_mean_diff: None,
                        skip_reason: String::new(),
                        s2_trend: None,
                        s2_prior: None,
                        s2_posterior: None,
                        df_prior: None,
                        df_total: None,
                        f_statistic: None,
                        f_p_value: None,
                        f_bh_q: None,
                        lfc_threshold: None,
                        n_peptides_observed: None,
                        peptide_variance_ratio: None,
                        ridge_lambda: None,
                        method: "ols".into(),
                        posthoc_method: "dunnett".into(),
                        posthoc_p: if raw_p.is_finite() { Some(raw_p) } else { None },
                        posthoc_adj_p: if adj_p.is_finite() { Some(adj_p) } else { None },
                    });
                }
            }
            OlsOutcome::Skipped { reason, n } => {
                let reason_str = match reason {
                    SkipReason::InsufficientPairs => "insufficient_samples",
                    SkipReason::ZeroVariance => "zero_variance",
                    SkipReason::NonFiniteInput => "non_finite_input",
                };
                for (a, b) in &contrast_list {
                    all_rows.push(DeResultRow {
                        panel: panel.clone(),
                        assay_id: assay_id.clone(),
                        gene_symbol: gene.clone(),
                        uniprot: uniprot.join(","),
                        comparison: contrast_label_for(&(a.clone(), b.clone())),
                        n_pairs: n,
                        mean_a: None, mean_b: None, mean_diff: None,
                        t: None, df: None, p_value: None, bh_q: None,
                        effect_size: None,
                        effect_size_method: "ols-posthoc-dunnett".into(),
                        ci_low: None, ci_high: None,
                        wilcoxon_p: None, wilcoxon_method: String::new(),
                        median_diff: None, trimmed_mean_diff: None,
                        skip_reason: reason_str.into(),
                        s2_trend: None, s2_prior: None, s2_posterior: None,
                        df_prior: None, df_total: None,
                        f_statistic: None, f_p_value: None, f_bh_q: None,
                        lfc_threshold: None,
                        n_peptides_observed: None, peptide_variance_ratio: None,
                        ridge_lambda: None,
                        method: "ols".into(),
                        posthoc_method: "dunnett".into(),
                        posthoc_p: None, posthoc_adj_p: None,
                    });
                }
            }
        }
    }

    let mut by_family: BTreeMap<(String, String), Vec<usize>> = BTreeMap::new();
    for (i, r) in all_rows.iter().enumerate() {
        by_family
            .entry((r.comparison.clone(), r.panel.clone()))
            .or_default()
            .push(i);
    }
    for indices in by_family.values() {
        let ps: Vec<Option<f64>> = indices.iter().map(|&i| all_rows[i].p_value).collect();
        let qs = atman_core::bh_fdr(&ps);
        for (j, &i) in indices.iter().enumerate() {
            all_rows[i].bh_q = qs[j];
        }
    }
    all_rows.sort_by(|a, b| {
        a.comparison
            .cmp(&b.comparison)
            .then_with(|| a.panel.cmp(&b.panel))
            .then_with(|| match (a.posthoc_adj_p, b.posthoc_adj_p) {
                (Some(x), Some(y)) => x.partial_cmp(&y).unwrap_or(std::cmp::Ordering::Equal),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            })
            .then_with(|| a.gene_symbol.cmp(&b.gene_symbol))
    });

    let results_path = args.output_dir.join("de_results.tsv");
    write_de_results(&results_path, &all_rows)?;
    let outputs = vec![results_path.clone()];
    eprintln!(
        "de posthoc dunnett: factor={:?} control={:?} m={} proteins={} rows={}",
        factor_name,
        control,
        m_contrasts,
        measured_features.len(),
        all_rows.len()
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &[
            "qc_measurements.tsv",
            "measurements.tsv",
            "samples.tsv",
            "proteins.tsv",
        ],
    )?;
    let sidecar = sidecar_path_for(&results_path);
    write_run_sidecar(
        &sidecar,
        "de",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output-dir": args.output_dir.display().to_string(),
            "test": "ols",
            "design": args.design,
            "groups": args.groups,
            "min-pairs": args.min_pairs,
            "post-hoc": args.post_hoc,
            "post-hoc-factor": factor_name,
            "dunnett-control": control,
            "dunnett-rho": if is_unbalanced { serde_json::Value::Null } else { serde_json::Value::from(RHO) },
            "dunnett-unbalanced": is_unbalanced,
            "dunnett-n-control": n_control,
            "dunnett-n-treatments": n_treatments,
            "dunnett-hsu-n-mc": if is_unbalanced { serde_json::Value::from(HSU_N_MC) } else { serde_json::Value::Null },
            "alpha": args.alpha,
            "omnibus-factor": args.omnibus_factor,
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        None,
)?;
    eprintln!("de posthoc dunnett: sidecar={}", sidecar.display());
    Ok(())
}

struct PosthocDesign {
    /// One row per retained sample: `(sample_id, design_row)`.
    rows: Vec<(String, Vec<f64>)>,
    /// Column labels aligned with each `design_row`.
    labels: Vec<String>,
    /// Covariate classification aligned with `setup.cov_names`.
    cov_kinds: Vec<CovKind>,
}

fn build_posthoc_full_design(
    samples: &[Sample],
    setup: &OlsSetup,
    _factor_name: &str,
) -> Result<PosthocDesign> {
    // `cov_raw` is keyed by `sample_id`; the per-key value is a
    // `Vec<Option<String>>` aligned with `setup.cov_names` order.
    let mut kept: Vec<(&Sample, usize, Vec<String>)> = Vec::new();
    for (idx, s) in samples.iter().enumerate() {
        if s.is_control {
            continue;
        }
        let per_cov = match setup.cov_raw.get(&s.sample_id) {
            Some(v) => v,
            None => continue,
        };
        let mut vals: Vec<String> = Vec::with_capacity(setup.cov_names.len());
        let mut complete = true;
        for (ci, _name) in setup.cov_names.iter().enumerate() {
            match per_cov.get(ci).cloned().unwrap_or(None) {
                Some(v) if !v.is_empty() => vals.push(v),
                _ => {
                    complete = false;
                    break;
                }
            }
        }
        if complete {
            kept.push((s, idx, vals));
        }
    }
    if kept.is_empty() {
        anyhow::bail!(
            "no samples with complete values for every --design covariate"
        );
    }

    let raw_vals_opt: Vec<Vec<Option<String>>> = kept
        .iter()
        .map(|(_, _, v)| v.iter().map(|s| Some(s.clone())).collect())
        .collect();
    let cov_kinds = classify_covariates(&setup.cov_names, &raw_vals_opt);

    let mut design_labels = vec!["(Intercept)".to_string()];
    for (idx, name) in setup.cov_names.iter().enumerate() {
        design_labels.extend(cov_kinds[idx].col_labels(name));
    }

    let mut rows: Vec<(String, Vec<f64>)> = Vec::with_capacity(kept.len());
    for (s, _idx, vals) in kept {
        let mut row: Vec<f64> = Vec::with_capacity(design_labels.len());
        row.push(1.0);
        for (i, name) in setup.cov_names.iter().enumerate() {
            let val_opt = Some(vals[i].clone());
            push_encoded_covariate(&mut row, name, &cov_kinds[i], &val_opt)?;
        }
        rows.push((s.sample_id.clone(), row));
    }
    Ok(PosthocDesign { rows, labels: design_labels, cov_kinds })
}

/// Build the contrast weight vector `c` (length `p`) for the
/// contrast `level_a − level_b` under ref-level encoding.
fn build_posthoc_contrast_vector(
    level_a: &str,
    level_b: &str,
    ref_level: &str,
    factor_col_by_level: &HashMap<String, usize>,
    p: usize,
) -> Vec<f64> {
    let mut c = vec![0.0; p];
    if level_a == ref_level && level_b != ref_level {
        // estimate = mean(ref) − mean(b) = 0 − β_b
        c[factor_col_by_level[level_b]] = -1.0;
    } else if level_b == ref_level && level_a != ref_level {
        // estimate = mean(a) − mean(ref) = β_a − 0
        c[factor_col_by_level[level_a]] = 1.0;
    } else if level_a != ref_level && level_b != ref_level {
        c[factor_col_by_level[level_a]] = 1.0;
        c[factor_col_by_level[level_b]] = -1.0;
    }
    // level_a == level_b == ref_level is impossible — parse_contrast_list
    // rejects equal levels — but if reached, all-zero c returns
    // estimate=0 gracefully.
    c
}

pub(super) fn write_design_rows(path: &Path, rows: &[DesignReportRow]) -> Result<()> {
    let mut buf =
        String::from("comparison\tsample_id\tcondition\tincluded\tdrop_reason\tcolumns\tvalues\n");
    for r in rows {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.comparison,
            r.sample_id,
            r.condition,
            if r.included { 1 } else { 0 },
            r.drop_reason,
            r.columns,
            r.values,
        ));
    }
    atomic_write(path, buf.as_bytes())
}
