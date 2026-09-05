//! Differential abundance command.
//!
//! `mod.rs` is the dispatch + orchestration layer: it parses/validates
//! `Args`, routes to the per-method submodules, builds the result rows,
//! and writes the run sidecar. The cohesive pieces live alongside:
//!
//! * [`args`] — the clap `Args` struct and `--per-subject-proxy` rewrite.
//! * [`ols_design`] — design-matrix construction, formula parsing,
//!   comparison resolution, covariate encoding, and the OLS bridge.
//! * [`output_rows`] — sidecar row shapes + TSV writers and the
//!   `--paired-by` subject-id override.

mod args;
mod ensemble;
mod limma;
mod moderated;
mod msqrob;
mod ols_design;
mod output_rows;
mod posthoc;
mod preflight;
mod robust_stats;

pub use args::Args;

use args::apply_per_subject_proxy;
use ensemble::run_ensemble;
use limma::run_limma;
use moderated::apply_moderated_shrinkage;
use msqrob::run_msqrob;
use ols_design::{
    augment_ols_design_with_external, build_ols_design, build_ols_design_continuous,
    build_ols_setup, ols_to_paired_t_result, paired_t_covariate_adjusted, raw_group_means,
    resolve_comparisons, validate_random_intercept_design, OlsDesign,
};
use output_rows::{
    apply_subset, override_condition, override_subject_id, write_covariate_rows,
    write_omnibus_rows, write_proxy_summary, CovariateRow, OmnibusRow, ReportAccumulator,
};
use posthoc::{run_posthoc_dunnett, run_posthoc_sidak, run_posthoc_tukey, write_design_rows};
use robust_stats::{median, robust_paired, robust_unpaired, RobustStats};

// Re-exports so sibling submodules can keep referring to these via
// `super::…`. These items now live in `ols_design`/`output_rows` but are
// part of the shared `de`-internal surface. A plain `use` is sufficient:
// private items in `mod.rs` are visible to descendant submodules.
use ols_design::{
    build_two_group_design, classify_covariates, covariate_names_from_terms, parse_design_terms,
    push_encoded_covariate, read_covariate_columns, CovKind, DesignTerm, OlsSetup,
};
use output_rows::DesignReportRow;

use anyhow::{Context, Result};
use atman_core::de::{
    bh_fdr, mixed_random_intercept, ols, paired_t, welch_t, OlsOutcome, PairedTResult, SkipReason,
};
use serde_json::{json, Map as JsonMap};
use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;
use std::time::SystemTime;

use super::parse_comparisons;
use crate::io::{
    hash_canonical_inputs, hash_labeled_inputs, read_measurements_long, read_proteins,
    read_samples, sidecar_path_for, write_de_report, write_de_results, write_run_sidecar,
    DeReportRow, DeResultRow,
};

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    let mut args = args;
    apply_per_subject_proxy(&mut args)?;
    preflight::check_below_lod_gate(
        &args.input_dir.join("measurements.tsv"),
        args.max_below_lod_fraction,
        args.allow_censored,
    )?;
    if args.test != "paired-t"
        && args.test != "moderated"
        && args.test != "welch-t"
        && args.test != "ols"
        && args.test != "mixed"
        && args.test != "limma"
        && args.test != "msqrob"
        && args.test != "ensemble"
    {
        anyhow::bail!(
            "test {:?} not supported; use paired-t, moderated, welch-t, ols, mixed, limma, msqrob, or ensemble",
            args.test
        );
    }
    let is_unpaired = args.test == "welch-t"
        || args.test == "ols"
        || args.test == "mixed"
        || args.test == "limma"
        || args.test == "msqrob"
        || args.test == "ensemble";
    // `--paired-by` accepts any column in samples.tsv; the value of that column
    // overrides the in-memory `subject_id` field below. The canonical default
    // remains `subject_id`. The legacy alias `participant` is handled the same
    // way when that column is present.
    if args.min_pairs < 2 {
        anyhow::bail!("min-pairs must be >= 2");
    }
    if !args.report_q_strict.is_finite() || !(0.0..1.0).contains(&args.report_q_strict) {
        anyhow::bail!(
            "--report-q-strict {} must be in (0, 1)",
            args.report_q_strict
        );
    }
    if !args.report_q_relaxed.is_finite() || !(0.0..1.0).contains(&args.report_q_relaxed) {
        anyhow::bail!(
            "--report-q-relaxed {} must be in (0, 1)",
            args.report_q_relaxed
        );
    }
    if args.report_q_strict >= args.report_q_relaxed {
        anyhow::bail!(
            "--report-q-strict {} must be < --report-q-relaxed {}",
            args.report_q_strict,
            args.report_q_relaxed
        );
    }
    if !args.proxy_p_threshold.is_finite() || !(0.0..1.0).contains(&args.proxy_p_threshold) {
        anyhow::bail!(
            "--proxy-p-threshold {} must be in (0, 1)",
            args.proxy_p_threshold
        );
    }
    if args.test == "moderated"
        && (!args.moderation_prior_df.is_finite() || args.moderation_prior_df <= 0.0)
    {
        anyhow::bail!("moderation-prior-df must be > 0 for moderated test");
    }
    if args.test == "limma" && args.robust {
        if !(0.0..0.5).contains(&args.limma_winsor_lower) || !args.limma_winsor_lower.is_finite() {
            anyhow::bail!(
                "--limma-winsor-lower {} must be in (0, 0.5)",
                args.limma_winsor_lower
            );
        }
        if !(0.0..0.5).contains(&args.limma_winsor_upper) || !args.limma_winsor_upper.is_finite() {
            anyhow::bail!(
                "--limma-winsor-upper {} must be in (0, 0.5)",
                args.limma_winsor_upper
            );
        }
        if args.limma_winsor_lower + args.limma_winsor_upper >= 1.0 {
            anyhow::bail!(
                "--limma-winsor-lower + --limma-winsor-upper must be < 1.0; got {} + {} = {}",
                args.limma_winsor_lower,
                args.limma_winsor_upper,
                args.limma_winsor_lower + args.limma_winsor_upper
            );
        }
    }
    if args.test != "ols"
        && args.test != "mixed"
        && (args.covariates.is_some()
            || args.design.is_some()
            || args.contrast.is_some()
            || args.per_subject_proxy.is_some()
            || args.fixed.is_some()
            || args.random.is_some())
    {
        anyhow::bail!(
            "--covariates, --design, --contrast, --per-subject-proxy, --fixed, and --random require --test ols or mixed"
        );
    }
    if !args.adjust_for.is_empty()
        && args.test != "limma"
        && args.test != "msqrob"
        && args.test != "ols"
        && args.test != "mixed"
        && args.test != "welch-t"
        && args.test != "paired-t"
    {
        anyhow::bail!(
            "--adjust-for requires --test limma, msqrob, ols, mixed, welch-t, or paired-t; \
             got --test {:?}",
            args.test
        );
    }
    if args.covariates.is_some() && args.design.is_some() {
        anyhow::bail!("use either --covariates or --design, not both");
    }
    if args.test == "mixed" {
        if args.fixed.is_none() {
            anyhow::bail!("--fixed is required with --test mixed");
        }
        if args.random.as_deref() != Some("1|subject_id") {
            anyhow::bail!("--test mixed currently supports only --random '1|subject_id'");
        }
        if args.covariates.is_some() || args.design.is_some() {
            anyhow::bail!("use --fixed with --test mixed; --covariates and --design are for ols");
        }
    } else if args.fixed.is_some() || args.random.is_some() {
        anyhow::bail!("--fixed and --random are only valid with --test mixed");
    }
    if args.test == "msqrob" {
        if args.peptide_measurements.is_none() {
            anyhow::bail!("--peptide-measurements is required with --test msqrob");
        }
        if args.peptide_metadata.is_none() {
            anyhow::bail!("--peptide-metadata is required with --test msqrob");
        }
        if args.min_peptides < 2 {
            anyhow::bail!("--min-peptides must be >= 2");
        }
        // `ridge_lambda` is validated inside run_msqrob.
    } else if args.peptide_measurements.is_some() && args.test != "ensemble" {
        anyhow::bail!("--peptide-measurements requires --test msqrob");
    } else if args.peptide_metadata.is_some() && args.test != "limma" && args.test != "ensemble" {
        anyhow::bail!(
            "--peptide-metadata is accepted for --test msqrob, --test limma (DEqMS), or --test ensemble; \
             got --test {:?}",
            args.test
        );
    }
    if args.omnibus_factor.is_some() {
        if args.test != "ols" {
            anyhow::bail!(
                "--omnibus-factor currently requires --test ols; got --test {:?}",
                args.test
            );
        }
        if args.design.is_none() {
            anyhow::bail!(
                "--omnibus-factor requires --design (so the factor has an encoded column group)"
            );
        }
    }
    match args.post_hoc.as_str() {
        "" => {}
        "sidak" => {
            if args.test != "ols" {
                anyhow::bail!(
                    "--post-hoc sidak currently requires --test ols; got --test {:?}",
                    args.test
                );
            }
            if args.design.is_none() {
                anyhow::bail!("--post-hoc sidak requires --design");
            }
            if args.contrast_list.as_deref().unwrap_or("").is_empty() {
                anyhow::bail!("--post-hoc sidak requires --contrast-list");
            }
            if args
                .post_hoc_factor
                .as_deref()
                .or(args.omnibus_factor.as_deref())
                .unwrap_or("")
                .is_empty()
            {
                anyhow::bail!("--post-hoc sidak requires --post-hoc-factor or --omnibus-factor");
            }
            if !(0.0..=1.0).contains(&args.alpha) {
                anyhow::bail!("--alpha must be in [0, 1]");
            }
            if args.groups.is_some() {
                anyhow::bail!(
                    "--post-hoc sidak drives comparisons from --contrast-list; \
                     --groups is not allowed in this mode"
                );
            }
        }
        "tukey" => {
            if args.test != "ols" {
                anyhow::bail!(
                    "--post-hoc tukey currently requires --test ols; got --test {:?}",
                    args.test
                );
            }
            if args.design.is_none() {
                anyhow::bail!("--post-hoc tukey requires --design");
            }
            if args
                .post_hoc_factor
                .as_deref()
                .or(args.omnibus_factor.as_deref())
                .unwrap_or("")
                .is_empty()
            {
                anyhow::bail!("--post-hoc tukey requires --post-hoc-factor or --omnibus-factor");
            }
            if !(0.0..=1.0).contains(&args.alpha) {
                anyhow::bail!("--alpha must be in [0, 1]");
            }
            if args.groups.is_some() {
                anyhow::bail!(
                    "--post-hoc tukey compares all pairs within the factor; \
                     --groups is not allowed in this mode"
                );
            }
            // `--contrast-list` is allowed (caller can restrict to a
            // subset of pairs) but not required — tukey generates all
            // ordered pairs from observed factor levels when absent.
        }
        "dunnett" => {
            if args.test != "ols" {
                anyhow::bail!(
                    "--post-hoc dunnett currently requires --test ols; got --test {:?}",
                    args.test
                );
            }
            if args.design.is_none() {
                anyhow::bail!("--post-hoc dunnett requires --design");
            }
            if args
                .post_hoc_factor
                .as_deref()
                .or(args.omnibus_factor.as_deref())
                .unwrap_or("")
                .is_empty()
            {
                anyhow::bail!("--post-hoc dunnett requires --post-hoc-factor or --omnibus-factor");
            }
            if !(0.0..=1.0).contains(&args.alpha) {
                anyhow::bail!("--alpha must be in [0, 1]");
            }
            if args.groups.is_some() {
                anyhow::bail!(
                    "--post-hoc dunnett compares each non-control level to the control; \
                     --groups is not allowed in this mode"
                );
            }
        }
        other => anyhow::bail!("unknown --post-hoc {other:?}; supported: sidak, tukey, dunnett"),
    }
    if args.test == "ensemble" {
        for frac in [
            args.ensemble_provisional_fraction,
            args.ensemble_sign_fraction,
        ] {
            if !frac.is_finite() || !(0.0..=1.0).contains(&frac) {
                anyhow::bail!("ensemble fraction thresholds must be in [0, 1]");
            }
        }
        if args.ensemble_provisional_fraction > args.ensemble_sign_fraction {
            anyhow::bail!(
                "--ensemble-provisional-fraction ({}) must be <= --ensemble-sign-fraction ({})",
                args.ensemble_provisional_fraction,
                args.ensemble_sign_fraction,
            );
        }
        if !args.ensemble_q_threshold.is_finite()
            || !(0.0..=1.0).contains(&args.ensemble_q_threshold)
        {
            anyhow::bail!("--ensemble-q-threshold must be in [0, 1]");
        }
    }
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

    // Ensemble: fan out to each requested method via sub-dispatch,
    // then aggregate. Each sub-call writes its own sidecar into a
    // tempdir that is discarded after aggregation; only the ensemble
    // output directory receives the final files.
    if args.test == "ensemble" {
        return run_ensemble(args, started_at);
    }

    // Post-hoc Sidak: owns its own pipeline (OLS once per protein,
    // then N contrasts via contrast_inference, Sidak-adjusted within
    // the contrast list per protein). Bypasses the classical
    // --groups loop entirely.
    if args.post_hoc == "sidak" {
        return run_posthoc_sidak(args, started_at);
    }
    if args.post_hoc == "tukey" {
        return run_posthoc_tukey(args, started_at);
    }
    if args.post_hoc == "dunnett" {
        return run_posthoc_dunnett(args, started_at);
    }

    // Parse comparisons after samples are available below; formula OLS can
    // infer the reference condition from --contrast when exactly two
    // non-control conditions are present.

    // Read inputs.
    let measurements = read_measurements_long(&args.input_dir.join("measurements.tsv"))?;
    let mut samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    if args.paired_by != "subject_id" {
        override_subject_id(
            &mut samples,
            &args.input_dir.join("samples.tsv"),
            &args.paired_by,
        )?;
    }
    let samples_path = args.input_dir.join("samples.tsv");
    if args.condition_col != "condition" {
        override_condition(&mut samples, &samples_path, &args.condition_col)?;
    }
    let subset_preds = if args.subset.is_empty() {
        Vec::new()
    } else {
        crate::design::parse_predicates(&args.subset.join(";"))?
    };
    let n_subset_dropped = if subset_preds.is_empty() {
        0
    } else {
        apply_subset(&mut samples, &samples_path, &subset_preds)?
    };
    if args.include_controls {
        for s in samples.iter_mut() {
            s.is_control = false;
        }
    }
    if let Some(label) = &args.collapse_others {
        let groups = args
            .groups
            .as_deref()
            .ok_or_else(|| anyhow::anyhow!("--collapse-others requires --groups"))?;
        let named: std::collections::BTreeSet<String> = parse_comparisons(groups)?
            .into_iter()
            .flat_map(|(a, b)| [a, b])
            .collect();
        if !named.contains(label) {
            anyhow::bail!(
                "--collapse-others {:?} must be one side of a --groups comparison",
                label
            );
        }
        for s in samples.iter_mut() {
            if let Some(c) = &s.condition {
                if !named.contains(c) {
                    s.condition = Some(label.clone());
                }
            }
        }
    }
    if !(0.0..=1.0).contains(&args.max_missing_fraction) {
        anyhow::bail!("--max-missing-fraction must be in [0, 1]");
    }
    if args.collapse_genes != crate::design::CollapseGenes::None
        && (args.test == "limma" || args.test == "msqrob")
    {
        anyhow::bail!("--collapse-genes is not supported with --test limma or msqrob");
    }
    let proteins = read_proteins(&args.input_dir.join("proteins.tsv"))?;
    let comparisons = resolve_comparisons(
        args.groups.as_deref(),
        &samples,
        args.test.as_str(),
        args.design.as_deref().or(args.fixed.as_deref()),
        args.contrast.as_deref(),
    )?;
    if args.test == "mixed" {
        validate_random_intercept_design(&samples, &comparisons, args.min_pairs)?;
    }
    let has_continuous = comparisons.iter().any(|(_, b)| b.is_empty());
    if has_continuous && !args.adjust_for.is_empty() {
        anyhow::bail!("--adjust-for is not supported with a continuous --contrast");
    }

    // Lookup: sample_id → (subject, condition, is_control)
    let sample_by_id: HashMap<&str, &atman_core::Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    // Lookup: (panel, gene_symbol) → (first_assay_id, uniprot) for reporting.
    // Gene symbols are unique within a panel for Olink Explore; we pick the
    // first OlinkID encountered as the canonical AssayId for the row.
    let mut gene_meta: BTreeMap<(String, String), (String, Vec<String>)> = BTreeMap::new();
    for p in &proteins {
        if let Some(gene) = &p.gene_symbol {
            if let Some(panel) = &p.panel {
                gene_meta
                    .entry((panel.clone(), gene.clone()))
                    .or_insert_with(|| (p.assay_id.0.clone(), p.uniprot.clone()));
            }
        }
    }

    // Build per-(panel, gene_symbol, condition) → Vec<(subject_or_sample, value)>.
    // For paired-t / moderated the first element is subject_id (for the paired
    // join); for welch-t it is unused. For `ols` we build a separate
    // sample-level structure below and use `cells_by_sample` instead.
    type CellsByCondition = BTreeMap<String, Vec<(String, f64)>>;
    let mut cells: BTreeMap<(String, String), CellsByCondition> = BTreeMap::new();
    // (panel, gene, sample_id) → abundance — only built for `--test ols` so
    // per-protein fits can look up abundance per sample and drop QC-masked
    // cells without needing to pair by subject.
    let mut cells_by_sample: HashMap<(String, String, String), f64> = HashMap::new();
    // (panel, gene) → sample_id → [(assay_id, abundance)], collapsed below
    // according to --collapse-genes.
    type RawCells = BTreeMap<(String, String), BTreeMap<String, Vec<(String, f64)>>>;
    let mut raw_cells: RawCells = BTreeMap::new();
    let mut sample_info: HashMap<String, (String, String)> = HashMap::new();

    let n_total_measurements = measurements.len();
    let mut n_skip_qc_masked = 0usize;
    let mut n_skip_unknown_sample = 0usize;
    let mut n_skip_control_sample = 0usize;
    let mut n_skip_no_subject_id = 0usize;
    let mut n_skip_no_condition = 0usize;
    let mut n_skip_no_panel = 0usize;
    let mut n_skip_no_gene_symbol = 0usize;
    let mut n_kept_measurements = 0usize;

    for m in &measurements {
        let abundance = match m.effective_abundance() {
            Some(a) => a,
            None => {
                n_skip_qc_masked += 1;
                continue;
            }
        };
        let s = match sample_by_id.get(m.sample_id.as_str()) {
            Some(s) => *s,
            None => {
                n_skip_unknown_sample += 1;
                continue;
            }
        };
        if s.is_control {
            n_skip_control_sample += 1;
            continue;
        }
        let subject = match &s.subject_id {
            Some(id) => id.clone(),
            None => {
                n_skip_no_subject_id += 1;
                continue;
            }
        };
        let condition = match &s.condition {
            Some(c) => c.clone(),
            None => {
                n_skip_no_condition += 1;
                continue;
            }
        };
        let panel = match &m.panel {
            Some(p) => p.clone(),
            None => {
                n_skip_no_panel += 1;
                continue;
            }
        };
        let gene = match &m.gene_symbol {
            Some(g) => g.clone(),
            None => {
                n_skip_no_gene_symbol += 1;
                continue;
            }
        };
        n_kept_measurements += 1;
        sample_info
            .entry(m.sample_id.clone())
            .or_insert((subject, condition));
        raw_cells
            .entry((panel, gene))
            .or_default()
            .entry(m.sample_id.clone())
            .or_default()
            .push((m.assay_id.0.clone(), abundance));
    }

    // Collapse protein groups that share a gene symbol to one value per
    // sample (--collapse-genes), then fill the per-condition cells.
    // cells_by_sample is used by OLS, mixed, and any path routing to OLS
    // (welch-t + adjust-for). paired-t + adjust-for also needs per-sample
    // lookup to compute paired differences.
    let needs_by_sample = args.test == "ols"
        || args.test == "mixed"
        || (!args.adjust_for.is_empty() && (args.test == "welch-t" || args.test == "paired-t"));
    let collapse = args.collapse_genes;
    let mut n_genes_multi_assay = 0usize;
    let mut n_extra_assays = 0usize;
    let mut n_assays_total = 0usize;
    for ((panel, gene), by_sample) in raw_cells {
        let mut assay_counts: BTreeMap<String, usize> = BTreeMap::new();
        for values in by_sample.values() {
            for (assay, _) in values {
                *assay_counts.entry(assay.clone()).or_default() += 1;
            }
        }
        n_assays_total += assay_counts.len();
        if assay_counts.len() > 1 {
            n_genes_multi_assay += 1;
            n_extra_assays += assay_counts.len() - 1;
        }
        let representative = collapse
            .representative(&assay_counts)
            .cloned()
            .unwrap_or_default();
        if assay_counts.len() > 1 {
            if let Some(meta) = gene_meta.get_mut(&(panel.clone(), gene.clone())) {
                meta.0 = representative.clone();
            }
        }
        for (sample_id, values) in by_sample {
            let Some(value) = collapse.collapse(&values, &representative) else {
                continue;
            };
            let (subject, condition) = sample_info[&sample_id].clone();
            if needs_by_sample {
                cells_by_sample.insert((panel.clone(), gene.clone(), sample_id.clone()), value);
            }
            cells
                .entry((panel.clone(), gene.clone()))
                .or_default()
                .entry(condition)
                .or_default()
                .push((subject, value));
        }
    }
    let n_genes_after_collapse = cells.len();
    if n_genes_multi_assay > 0 {
        eprintln!(
            "de: {} gene symbols are carried by more than one assay ({} extra assays); --collapse-genes {}",
            n_genes_multi_assay,
            n_extra_assays,
            collapse.as_str()
        );
    }

    let mixed_fixed_formula = args.fixed.as_ref().map(|fixed| format!("~ {fixed}"));
    // welch-t + adjust-for routes to OLS internally (plain OLS SE, not HC3).
    // paired-t + adjust-for uses a dedicated paired-diff OLS path (see below).
    let is_ols_routed = args.test == "ols"
        || args.test == "mixed"
        || (!args.adjust_for.is_empty() && args.test == "welch-t");
    let model_setup = if is_ols_routed {
        Some(build_ols_setup(
            &args.input_dir.join("samples.tsv"),
            if args.test == "mixed" {
                mixed_fixed_formula.as_deref()
            } else {
                args.design.as_deref()
            },
            if args.test == "ols" {
                args.covariates.as_deref()
            } else {
                None
            },
            args.contrast.as_deref(),
        )?)
    } else {
        None
    };

    // For each comparison, run paired_t on every (panel, gene) and collect.
    // Results are accumulated in a flat Vec for the final TSV + a per-family
    // Vec for BH-FDR correction. A family = one comparison.
    let mut all_rows: Vec<DeResultRow> = Vec::new();
    // For the report sidecar.
    let mut report_rows: Vec<DeReportRow> = Vec::new();

    // Per-covariate estimates across all proteins × comparisons, emitted to
    // de_covariates.tsv when --test ols. One row per (comparison, panel,
    // gene, covariate) with beta, se, t, p for the non-intercept and
    // non-group columns.
    let mut covariate_rows: Vec<CovariateRow> = Vec::new();
    let mut omnibus_rows: Vec<OmnibusRow> = Vec::new();
    let mut design_rows: Vec<DesignReportRow> = Vec::new();

    // Load and merge all --adjust-for external covariate files.
    // Each file is wide-format: sample_id | cov1 | cov2 | …
    // Multiple files are joined; overlapping covariate names are an error.
    // Supported for --test limma and --test msqrob (gate enforced above).
    let external_covariates: Option<
        std::collections::BTreeMap<String, std::collections::BTreeMap<String, f64>>,
    > = if args.adjust_for.is_empty() {
        None
    } else {
        use atman_core::de::{join_external_covariates, read_external_covariates};
        let mut merged: std::collections::BTreeMap<
            String,
            std::collections::BTreeMap<String, f64>,
        > = std::collections::BTreeMap::new();
        for path in &args.adjust_for {
            let ext = read_external_covariates(path)
                .map_err(|e| anyhow::anyhow!("--adjust-for {:?}: {}", path, e))?;
            // Merge into `merged`: for each sample, add its covariate columns.
            for (sample_id, covs) in ext {
                let entry = merged.entry(sample_id).or_default();
                for (cov_name, cov_val) in covs {
                    if entry.contains_key(&cov_name) {
                        anyhow::bail!(
                            "--adjust-for: duplicate covariate name {:?} found in {:?}",
                            cov_name,
                            path
                        );
                    }
                    entry.insert(cov_name, cov_val);
                }
            }
        }
        // Validate: every non-control sample with a condition must be in the merged map.
        let relevant_samples: Vec<&atman_core::Sample> = samples
            .iter()
            .filter(|s| !s.is_control && s.condition.is_some())
            .collect();
        let ref_samples: Vec<atman_core::Sample> =
            relevant_samples.iter().map(|s| (*s).clone()).collect();
        join_external_covariates(&ref_samples, &merged)
            .map_err(|e| anyhow::anyhow!("--adjust-for join failed: {}", e))?;
        Some(merged)
    };

    if args.test == "limma" {
        let (limma_rows, limma_reports) = run_limma(
            &args,
            &samples,
            &proteins,
            &measurements,
            &comparisons,
            external_covariates.as_ref(),
        )?;
        all_rows.extend(limma_rows);
        report_rows.extend(limma_reports);
    }

    if args.test == "msqrob" {
        let (msqrob_rows, msqrob_reports) = run_msqrob(
            &args,
            &samples,
            &proteins,
            &comparisons,
            external_covariates.as_ref(),
        )?;
        all_rows.extend(msqrob_rows);
        report_rows.extend(msqrob_reports);
    }

    let mut missingness_by_comparison = JsonMap::new();
    if args.test != "limma" && args.test != "msqrob" {
        for (comp_a, comp_b) in &comparisons {
            let continuous = comp_b.is_empty();
            let comparison_label = if continuous {
                comp_a.clone()
            } else {
                format!("{}-{}", comp_a, comp_b)
            };
            let n_comparison_samples = samples
                .iter()
                .filter(|s| {
                    !s.is_control
                        && match s.condition.as_deref() {
                            Some(c) => continuous || c == comp_a || c == comp_b,
                            None => false,
                        }
                })
                .count();
            let mut n_missingness_dropped = 0usize;
            let mut n_missingness_retained = 0usize;
            // Per-family p-value vector aligned with `family_rows` order.
            let mut family_p: Vec<Option<f64>> = Vec::new();
            let mut family_rows: Vec<DeResultRow> = Vec::new();

            // OLS mode: precompute the encoded design matrix (excluding y) for
            // this comparison so every per-protein fit reuses it, only swapping
            // the y vector (with complete-case filtering on non-finite abundance).
            // welch-t + adjust-for also routes here (plain OLS SE).
            let model_design: Option<OlsDesign> = if is_ols_routed {
                let setup = model_setup
                    .as_ref()
                    .expect("model_setup populated when is_ols_routed");
                let mut design = if continuous {
                    build_ols_design_continuous(&samples, setup)?
                } else {
                    build_ols_design(&samples, comp_a, comp_b, setup)?
                };
                // Augment with --adjust-for external covariate columns.
                if let Some(ext) = external_covariates.as_ref() {
                    augment_ols_design_with_external(&mut design, ext, comp_a, comp_b)?;
                }
                design_rows.extend(design.report_rows.clone());
                Some(design)
            } else {
                None
            };

            for ((panel, gene), by_condition) in &cells {
                if args.max_missing_fraction < 1.0 {
                    let present = if continuous {
                        by_condition.values().map(Vec::len).sum()
                    } else {
                        by_condition.get(comp_a).map(Vec::len).unwrap_or(0)
                            + by_condition.get(comp_b).map(Vec::len).unwrap_or(0)
                    };
                    let missing = 1.0 - present as f64 / n_comparison_samples.max(1) as f64;
                    if missing > args.max_missing_fraction + 1e-12 {
                        n_missingness_dropped += 1;
                        continue;
                    }
                }
                n_missingness_retained += 1;
                let (assay_id, uniprot) = match gene_meta.get(&(panel.clone(), gene.clone())) {
                    Some(meta) => meta.clone(),
                    None => (String::new(), vec![]),
                };

                // Build per-condition cell vectors. In paired mode we intersect
                // by subject_id; in unpaired (welch-t) mode we treat each
                // condition's samples as an independent group and drop the
                // subject linkage.
                let empty: Vec<(String, f64)> = Vec::new();
                let va = by_condition.get(comp_a).unwrap_or(&empty);
                let vb = by_condition.get(comp_b).unwrap_or(&empty);

                let (result, robust) = if args.test == "mixed" {
                    let design = model_design
                        .as_ref()
                        .expect("model_design populated when test=mixed");
                    let mut y: Vec<f64> = Vec::with_capacity(design.rows.len());
                    let mut rows: Vec<Vec<f64>> = Vec::with_capacity(design.rows.len());
                    let mut groups: Vec<String> = Vec::with_capacity(design.rows.len());
                    for (sid, row) in design.rows.iter() {
                        if let Some(abundance) =
                            cells_by_sample.get(&(panel.clone(), gene.clone(), sid.clone()))
                        {
                            let sample = sample_by_id
                                .get(sid.as_str())
                                .expect("design sample exists in sample_by_id");
                            y.push(*abundance);
                            rows.push(row.clone());
                            groups.push(sample.subject_id.as_ref().unwrap_or(sid).clone());
                        }
                    }
                    let fit = mixed_random_intercept(&rows, &y, &groups, args.min_pairs);
                    let (mean_a_raw, mean_b_raw) = raw_group_means(&y, &rows, design.group_col);
                    (
                        ols_to_paired_t_result(
                            fit,
                            design,
                            mean_a_raw,
                            mean_b_raw,
                            panel,
                            gene,
                            &comparison_label,
                            &mut covariate_rows,
                        ),
                        RobustStats::default(),
                    )
                } else if args.test == "ols" || (is_ols_routed && args.test == "welch-t") {
                    // OLS path: also handles welch-t + --adjust-for (routed to
                    // plain OLS; HC3 robust SE is not available without new
                    // dependencies — documented in DONE_WITH_CONCERNS K6).
                    let design = model_design
                        .as_ref()
                        .expect("ols_design populated when is_ols_routed");
                    // Look up per-sample abundance for this protein; drop
                    // samples where the cell is QC-masked (missing from
                    // cells_by_sample).
                    let mut y: Vec<f64> = Vec::with_capacity(design.rows.len());
                    let mut rows: Vec<Vec<f64>> = Vec::with_capacity(design.rows.len());
                    for (sid, row) in design.rows.iter() {
                        if let Some(abundance) =
                            cells_by_sample.get(&(panel.clone(), gene.clone(), sid.clone()))
                        {
                            y.push(*abundance);
                            rows.push(row.clone());
                        }
                    }
                    let fit = ols(&rows, &y, args.min_pairs);
                    // Pooled-SD Cohen d over the fitted subjects when the
                    // group column is the 0/1 condition indicator.
                    let robust = if design.group_is_indicator {
                        let mut ga = Vec::new();
                        let mut gb = Vec::new();
                        for (row, value) in rows.iter().zip(y.iter()) {
                            if row[design.group_col] == 1.0 {
                                ga.push(*value);
                            } else {
                                gb.push(*value);
                            }
                        }
                        let d = atman_core::contrast::cohen_d_pooled(&ga, &gb);
                        let ci =
                            d.and_then(|d| atman_core::contrast::cohen_d_ci(d, ga.len(), gb.len()));
                        RobustStats {
                            effect_size: d,
                            effect_size_method: if d.is_some() {
                                "cohen_d".to_string()
                            } else {
                                String::new()
                            },
                            ci_low: ci.map(|c| c.0),
                            ci_high: ci.map(|c| c.1),
                            n_a: Some(ga.len()),
                            n_b: Some(gb.len()),
                            ..RobustStats::default()
                        }
                    } else {
                        RobustStats::default()
                    };
                    // Omnibus F-test on the --omnibus-factor's columns,
                    // captured before the fit is consumed downstream. The
                    // factor is identified by matching design_labels
                    // starting with its name (mirrors how categorical
                    // one-hot columns are labelled in OlsDesign).
                    if let (OlsOutcome::Computed(f), Some(factor_name)) =
                        (&fit, args.omnibus_factor.as_deref())
                    {
                        let factor_cols: Vec<usize> = design
                            .design_labels
                            .iter()
                            .enumerate()
                            .filter(|(idx, l)| {
                                *idx != 0 && *idx != design.group_col && l.starts_with(factor_name)
                            })
                            .map(|(i, _)| i)
                            .collect();
                        if factor_cols.len() >= 2 {
                            if let Some(om) = atman_core::de::omnibus_f_test(
                                &rows,
                                &f.beta,
                                &factor_cols,
                                f.sigma2,
                                f.df,
                            ) {
                                omnibus_rows.push(OmnibusRow {
                                    panel: panel.clone(),
                                    assay_id: assay_id.clone(),
                                    gene_symbol: gene.clone(),
                                    factor: factor_name.to_string(),
                                    comparison: comparison_label.clone(),
                                    f_statistic: om.f_statistic,
                                    df_num: om.df_num,
                                    df_den: om.df_den,
                                    p_value: om.p_value,
                                    bh_q: None,
                                });
                            }
                        }
                    }
                    // Compute raw group means for schema compatibility with
                    // paired-t / welch-t (undefined for a continuous contrast).
                    let (mean_a_raw, mean_b_raw) = if design.group_is_indicator {
                        raw_group_means(&y, &rows, design.group_col)
                    } else {
                        (None, None)
                    };
                    (
                        ols_to_paired_t_result(
                            fit,
                            design,
                            mean_a_raw,
                            mean_b_raw,
                            panel,
                            gene,
                            &comparison_label,
                            &mut covariate_rows,
                        ),
                        robust,
                    )
                } else if args.test == "paired-t" && !args.adjust_for.is_empty() {
                    // paired-t + --adjust-for: paired-difference ANCOVA.
                    // Computes d_k = abundance_b_k - abundance_a_k per subject,
                    // covariate_diff_k = cov_b_k - cov_a_k, then fits
                    // OLS: d ~ 1 + cov_diff_1 + cov_diff_2 + ...
                    // The intercept is the adjusted mean paired difference.
                    let ext = external_covariates
                        .as_ref()
                        .expect("external_covariates set");
                    let result = paired_t_covariate_adjusted(
                        va,
                        vb,
                        ext,
                        &sample_by_id,
                        comp_a,
                        comp_b,
                        args.min_pairs,
                    );
                    (result, RobustStats::default())
                } else if is_unpaired {
                    let a_vals: Vec<f64> = va.iter().map(|(_, v)| *v).collect();
                    let b_vals: Vec<f64> = vb.iter().map(|(_, v)| *v).collect();
                    (
                        welch_t(&a_vals, &b_vals, args.min_pairs),
                        robust_unpaired(&a_vals, &b_vals, args.min_pairs),
                    )
                } else {
                    // Subject → value, for fast paired join.
                    let a_by_subj: HashMap<&str, f64> =
                        va.iter().map(|(s, v)| (s.as_str(), *v)).collect();
                    let b_by_subj: HashMap<&str, f64> =
                        vb.iter().map(|(s, v)| (s.as_str(), *v)).collect();
                    let pairs: Vec<(f64, f64)> = a_by_subj
                        .iter()
                        .filter_map(|(subj, a)| b_by_subj.get(subj).map(|b| (*a, *b)))
                        .collect();
                    (
                        paired_t(&pairs, args.min_pairs),
                        robust_paired(&pairs, args.min_pairs),
                    )
                };
                let row = match &result {
                    PairedTResult::Computed {
                        n_pairs,
                        mean_a,
                        mean_b,
                        mean_diff,
                        t,
                        df,
                        p_value,
                    } => DeResultRow {
                        panel: panel.clone(),
                        assay_id: assay_id.clone(),
                        gene_symbol: gene.clone(),
                        uniprot: uniprot.join(","),
                        comparison: comparison_label.clone(),
                        n_pairs: *n_pairs,
                        mean_a: Some(*mean_a).filter(|v| v.is_finite()),
                        mean_b: Some(*mean_b).filter(|v| v.is_finite()),
                        mean_diff: Some(*mean_diff),
                        t: Some(*t),
                        df: Some(*df),
                        p_value: Some(*p_value),
                        bh_q: None, // filled after BH sweep
                        effect_size: robust.effect_size,
                        effect_size_method: robust.effect_size_method.clone(),
                        ci_low: robust.ci_low,
                        ci_high: robust.ci_high,
                        wilcoxon_p: robust.wilcoxon_p,
                        wilcoxon_method: robust.wilcoxon_method.clone(),
                        median_diff: robust.median_diff,
                        trimmed_mean_diff: robust.trimmed_mean_diff,
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
                        method: args.test.clone(),
                        posthoc_method: String::new(),
                        posthoc_p: None,
                        posthoc_adj_p: None,
                        n_a: robust.n_a,
                        n_b: robust.n_b,
                    },
                    PairedTResult::Skipped { reason, n_pairs } => DeResultRow {
                        panel: panel.clone(),
                        assay_id: assay_id.clone(),
                        gene_symbol: gene.clone(),
                        uniprot: uniprot.join(","),
                        comparison: comparison_label.clone(),
                        n_pairs: *n_pairs,
                        mean_a: None,
                        mean_b: None,
                        mean_diff: None,
                        t: None,
                        df: None,
                        p_value: None,
                        bh_q: None,
                        effect_size: robust.effect_size,
                        effect_size_method: robust.effect_size_method.clone(),
                        ci_low: robust.ci_low,
                        ci_high: robust.ci_high,
                        wilcoxon_p: robust.wilcoxon_p,
                        wilcoxon_method: robust.wilcoxon_method.clone(),
                        median_diff: robust.median_diff,
                        trimmed_mean_diff: robust.trimmed_mean_diff,
                        skip_reason: match reason {
                            SkipReason::InsufficientPairs => "insufficient_pairs".to_string(),
                            SkipReason::ZeroVariance => "zero_variance".to_string(),
                            SkipReason::NonFiniteInput => "non_finite_input".to_string(),
                        },
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
                        method: args.test.clone(),
                        posthoc_method: String::new(),
                        posthoc_p: None,
                        posthoc_adj_p: None,
                        n_a: robust.n_a,
                        n_b: robust.n_b,
                    },
                };
                family_p.push(row.p_value);
                family_rows.push(row);
            }

            if args.test == "moderated" {
                apply_moderated_shrinkage(&mut family_rows, args.moderation_prior_df)?;
                family_p = family_rows.iter().map(|r| r.p_value).collect();
            }

            // BH-FDR within this comparison family.
            let qs = bh_fdr(&family_p);
            for (row, q) in family_rows.iter_mut().zip(qs.into_iter()) {
                row.bh_q = q;
            }

            // Per-panel report summary.
            let mut per_panel: BTreeMap<String, ReportAccumulator> = BTreeMap::new();
            for row in &family_rows {
                let acc = per_panel.entry(row.panel.clone()).or_default();
                acc.n_tests += 1;
                if let Some(q) = row.bh_q {
                    if q < args.report_q_strict {
                        acc.n_q_strict += 1;
                    }
                    if q < args.report_q_relaxed {
                        acc.n_q_relaxed += 1;
                    }
                    if acc.min_q.map(|m| q < m).unwrap_or(true) {
                        acc.min_q = Some(q);
                    }
                } else {
                    acc.n_skipped += 1;
                }
                if let Some(d) = row.mean_diff {
                    let abs_d = d.abs();
                    if acc.max_abs_effect.map(|m| abs_d > m).unwrap_or(true) {
                        acc.max_abs_effect = Some(abs_d);
                    }
                }
            }
            for (panel, acc) in per_panel {
                report_rows.push(DeReportRow {
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

            all_rows.extend(family_rows);
            missingness_by_comparison.insert(
                comparison_label.clone(),
                json!({
                    "n_samples": n_comparison_samples,
                    "n_retained": n_missingness_retained,
                    "n_dropped": n_missingness_dropped,
                }),
            );
        }
    } // end of gating the paired/ols/welch/mixed dispatch

    // Stable sort for inspection: (comparison, bh_q asc with None last, |mean_diff| desc).
    all_rows.sort_by(|a, b| {
        a.comparison
            .cmp(&b.comparison)
            .then_with(|| match (a.bh_q, b.bh_q) {
                (Some(qa), Some(qb)) => qa.partial_cmp(&qb).unwrap_or(std::cmp::Ordering::Equal),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            })
            .then_with(|| match (a.mean_diff, b.mean_diff) {
                (Some(da), Some(db)) => db
                    .abs()
                    .partial_cmp(&da.abs())
                    .unwrap_or(std::cmp::Ordering::Equal),
                _ => std::cmp::Ordering::Equal,
            })
            .then_with(|| a.gene_symbol.cmp(&b.gene_symbol))
    });

    let results_path = args.output_dir.join("de_results.tsv");
    let report_path = args.output_dir.join("de_report.tsv");
    let covariates_path = args.output_dir.join("de_covariates.tsv");
    let proxy_path = args.output_dir.join("de_proxy_summary.tsv");
    let design_path = args.output_dir.join("de_design.tsv");
    write_de_results(&results_path, &all_rows)?;
    write_de_report(&report_path, &report_rows)?;
    let mut outputs: Vec<PathBuf> = vec![results_path.clone(), report_path.clone()];
    if (args.test == "ols" || args.test == "mixed") && !covariate_rows.is_empty() {
        write_covariate_rows(&covariates_path, &covariate_rows)?;
        outputs.push(covariates_path.clone());
    }
    if let Some(proxy) = args.per_subject_proxy.as_deref() {
        write_proxy_summary(&proxy_path, proxy, &covariate_rows, args.proxy_p_threshold)?;
        outputs.push(proxy_path.clone());
    }
    if args.test == "ols" || args.test == "mixed" {
        write_design_rows(&design_path, &design_rows)?;
        outputs.push(design_path.clone());
    }
    if !omnibus_rows.is_empty() {
        // BH-adjust omnibus p-values within each (comparison, panel).
        use std::collections::BTreeMap;
        let mut indices_by_group: BTreeMap<(String, String), Vec<usize>> = BTreeMap::new();
        for (i, r) in omnibus_rows.iter().enumerate() {
            indices_by_group
                .entry((r.comparison.clone(), r.panel.clone()))
                .or_default()
                .push(i);
        }
        for indices in indices_by_group.values() {
            let ps: Vec<Option<f64>> = indices
                .iter()
                .map(|&i| {
                    let v = omnibus_rows[i].p_value;
                    if v.is_finite() {
                        Some(v)
                    } else {
                        None
                    }
                })
                .collect();
            let qs = atman_core::bh_fdr(&ps);
            for (j, &i) in indices.iter().enumerate() {
                omnibus_rows[i].bh_q = qs[j];
            }
        }
        let omnibus_path = args.output_dir.join("de_omnibus.tsv");
        write_omnibus_rows(&omnibus_path, &omnibus_rows)?;
        outputs.push(omnibus_path.clone());
        eprintln!(
            "de: omnibus factor={:?} rows={}",
            args.omnibus_factor.as_deref().unwrap_or(""),
            omnibus_rows.len()
        );
    }

    let computed = all_rows.iter().filter(|r| r.p_value.is_some()).count();
    let skipped = all_rows.len() - computed;
    eprintln!(
        "de: test={} comparisons={} rows={} computed={} skipped={} min_pairs={} prior_df={} design={}",
        args.test,
        comparisons.len(),
        all_rows.len(),
        computed,
        skipped,
        args.min_pairs,
        args.moderation_prior_df,
        model_setup
            .as_ref()
            .map(|s| s.label.as_str())
            .unwrap_or(""),
    );
    eprintln!(
        "de: measurement_filter total={} kept={} qc_masked={} unknown_sample={} control_sample={} no_subject_id={} no_condition={} no_panel={} no_gene_symbol={}",
        n_total_measurements,
        n_kept_measurements,
        n_skip_qc_masked,
        n_skip_unknown_sample,
        n_skip_control_sample,
        n_skip_no_subject_id,
        n_skip_no_condition,
        n_skip_no_panel,
        n_skip_no_gene_symbol,
    );

    let finished_at = SystemTime::now();
    let mut inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &["measurements.tsv", "samples.tsv", "proteins.tsv"],
    )?;
    // Hash each --adjust-for file and add to inputs so the sidecar
    // captures full chained-input provenance.
    for path in &args.adjust_for {
        let label = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("adjust_for");
        let labeled_entries: Vec<(&str, &std::path::Path)> = vec![(label, path.as_path())];
        let extra = hash_labeled_inputs(&labeled_entries)?;
        inputs_sha256.extend(extra);
    }
    let sidecar = sidecar_path_for(&results_path);
    let mut extras = JsonMap::new();
    extras.insert(
        "measurement_filter".into(),
        json!({
            "n_total": n_total_measurements,
            "n_kept": n_kept_measurements,
            "n_skip": {
                "qc_masked": n_skip_qc_masked,
                "unknown_sample": n_skip_unknown_sample,
                "control_sample": n_skip_control_sample,
                "no_subject_id": n_skip_no_subject_id,
                "no_condition": n_skip_no_condition,
                "no_panel": n_skip_no_panel,
                "no_gene_symbol": n_skip_no_gene_symbol,
            },
        }),
    );
    extras.insert(
        "missingness_filter".into(),
        serde_json::Value::Object(missingness_by_comparison),
    );
    extras.insert("n_subset_dropped".into(), json!(n_subset_dropped));
    extras.insert(
        "gene_symbol_collapse".into(),
        json!({
            "rule": collapse.as_str(),
            "n_assays": n_assays_total,
            "n_genes_with_multiple_assays": n_genes_multi_assay,
            "n_extra_assays": n_extra_assays,
            "n_genes_after_collapse": n_genes_after_collapse,
        }),
    );
    extras.insert(
        "report_thresholds".into(),
        json!({
            "q_strict": args.report_q_strict,
            "q_relaxed": args.report_q_relaxed,
            "proxy_p": args.proxy_p_threshold,
        }),
    );
    write_run_sidecar(
        &sidecar,
        "de",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output-dir": args.output_dir.display().to_string(),
            "test": args.test,
            "paired-by": args.paired_by,
            "covariates": args.covariates,
            "design": args.design,
            "contrast": args.contrast,
            "per-subject-proxy": args.per_subject_proxy,
            "fixed": args.fixed,
            "random": args.random,
            "groups": args.groups,
            "min-pairs": args.min_pairs,
            "moderation-prior-df": args.moderation_prior_df,
            "lfc-threshold": args.lfc_threshold,
            "trend": args.trend,
            "robust": args.robust,
            "limma-winsor-lower": args.limma_winsor_lower,
            "limma-winsor-upper": args.limma_winsor_upper,
            "peptide-measurements": args.peptide_measurements.as_ref().map(|p| p.display().to_string()),
            "peptide-metadata": args.peptide_metadata.as_ref().map(|p| p.display().to_string()),
            "ridge-lambda": args.ridge_lambda,
            "min-peptides": args.min_peptides,
            "omnibus-factor": args.omnibus_factor,
            "report-q-strict": args.report_q_strict,
            "report-q-relaxed": args.report_q_relaxed,
            "proxy-p-threshold": args.proxy_p_threshold,
            "adjust-for": args.adjust_for.iter().map(|p| p.display().to_string()).collect::<Vec<_>>(),
            "include-controls": args.include_controls,
            "condition-col": args.condition_col,
            "subset": args.subset,
            "collapse-others": args.collapse_others,
            "max-missing-fraction": args.max_missing_fraction,
            "collapse-genes": collapse.as_str(),
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        Some(extras),
    )?;
    eprintln!("de: sidecar={}", sidecar.display());
    Ok(())
}
