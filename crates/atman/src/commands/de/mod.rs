//! Differential abundance command.

use anyhow::{anyhow, Context, Result};
use atman_core::de::{
    bh_fdr, mixed_random_intercept, ols, paired_t, welch_t, OlsOutcome, PairedTResult, SkipReason,
};
use atman_core::stats::mean;
use atman_core::Sample;
use clap::Args as ClapArgs;
use serde_json::json;
use statrs::distribution::{ContinuousCDF, Normal, StudentsT};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use super::parse_comparisons;
use crate::io::{
    atomic_write, hash_canonical_inputs, read_measurements_long, read_proteins, read_samples,
    sidecar_path_for, write_de_report, write_de_results, write_run_sidecar, DeReportRow,
    DeResultRow,
};

#[derive(ClapArgs, Debug, Clone)]
pub struct Args {
    /// Directory containing qc_measurements.tsv, samples.tsv, proteins.tsv.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output directory for de_results.tsv + de_report.tsv.
    #[arg(long)]
    output_dir: PathBuf,

    /// Test kind: `paired-t` (default), `moderated` (paired variance
    /// shrinkage), `welch-t` (unpaired two-sample with Welch-Satterthwaite
    /// df), or `ols` (covariate-adjusted unpaired linear model).
    #[arg(long, default_value = "paired-t")]
    test: String,

    /// Biological-replicate key column on samples.tsv. Required for
    /// `paired-t` and `moderated`; ignored for `welch-t` and `ols`.
    #[arg(long, default_value = "participant")]
    paired_by: String,

    /// Comma-separated covariate column names from samples.tsv for
    /// `--test ols`. Numeric columns (parseable as f64 in every non-missing
    /// row) are used as continuous predictors; everything else is treated
    /// as categorical and one-hot encoded, dropping the alphabetically
    /// first level as reference. Samples with any missing covariate value
    /// are dropped before fitting. Ignored for other tests.
    #[arg(long)]
    covariates: Option<String>,

    /// Formula-style design for `--test ols`, for example
    /// `~ condition + age + sex + batch`. The intercept is included
    /// automatically. `condition` is the Atman sample condition; other terms
    /// are read from samples.tsv.
    #[arg(long)]
    design: Option<String>,

    /// Coefficient to test for `--test ols --design`, for example
    /// `conditionCase`. If omitted with --groups, Atman uses the comparison
    /// coefficient for each `A-B` group.
    #[arg(long)]
    contrast: Option<String>,

    /// Physiological subject-level proxy from samples.tsv to include as a continuous OLS regressor.
    #[arg(long)]
    per_subject_proxy: Option<String>,

    /// Fixed-effect terms for `--test mixed`, for example
    /// `condition + age + sex`. The intercept is included automatically.
    #[arg(long)]
    fixed: Option<String>,

    /// Random-effect structure for `--test mixed`. Initial support is
    /// random intercept by subject only: `1|subject_id`.
    #[arg(long)]
    random: Option<String>,

    /// Comma-separated comparisons in `A-B` form. Each is a separate
    /// hypothesis family for FDR. Example: "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2".
    #[arg(long)]
    groups: Option<String>,

    /// Minimum number of samples required per group for a test to run.
    /// For paired tests this is the number of matched pairs; for `welch-t`
    /// it is the minimum of `|group_a|` and `|group_b|`. Below this the
    /// protein is emitted as Skipped(InsufficientPairs) with NaN p/q.
    #[arg(long, default_value_t = 5)]
    min_pairs: usize,

    /// Prior degrees of freedom for `--test moderated` variance shrinkage.
    #[arg(long, default_value_t = 4.0)]
    moderation_prior_df: f64,

    /// Minimum log-fold-change threshold for TREAT testing under `--test limma`.
    /// `0.0` (default) reduces to the standard moderated-t p-value.
    #[arg(long, default_value_t = 0.0)]
    lfc_threshold: f64,

    /// Fit the mean-variance trend before eBayes shrinkage (limma only).
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set, num_args = 1)]
    trend: bool,

    /// Use the robust (Winsorized) prior fit (limma only).
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set, num_args = 1)]
    robust: bool,

    /// Peptide-level long TSV (one row per sample × peptide).
    /// Required for `--test msqrob`. Rejected for other tests.
    /// Schema: `sample_id, peptide_id, abundance, abundance_unit,
    /// dropped_by_qc, below_lod`.
    #[arg(long)]
    peptide_measurements: Option<PathBuf>,

    /// Peptide catalog TSV mapping `peptide_id → assay_id` (parent
    /// protein) plus optional `sequence, charge, modifications,
    /// missed_cleavages`. Required for `--test msqrob`.
    #[arg(long)]
    peptide_metadata: Option<PathBuf>,

    /// L2 ridge penalty on non-intercept fixed-effect coefficients
    /// for `--test msqrob`. Numeric value ≥ 0 or `auto` (currently
    /// equivalent to `0.0`; data-driven selection is a follow-on).
    #[arg(long, default_value = "auto")]
    ridge_lambda: String,

    /// Minimum distinct peptides observed per protein required for
    /// `--test msqrob`. Below this, the protein emits
    /// `Skipped(insufficient_peptides)`.
    #[arg(long, default_value_t = 2)]
    min_peptides: usize,

    /// Comma-separated list of methods to run under `--test ensemble`.
    /// Each of `paired-t`, `welch-t`, `ols`, `mixed`, `limma`, `msqrob`
    /// is valid. Methods whose required inputs are missing
    /// (`--peptide-measurements` for msqrob, a paired-subject layout
    /// for paired-t, etc.) are auto-skipped with a note in
    /// `methods_skipped` rather than aborting the run.
    #[arg(long, default_value = "welch-t,ols,limma,msqrob")]
    ensemble_methods: String,

    /// Per-method BH-q significance threshold for ensemble grading.
    #[arg(long, default_value_t = 0.05)]
    ensemble_q_threshold: f64,

    /// Fraction of applicable methods that must pass the q-threshold
    /// for a VALIDATED grade.
    #[arg(long, default_value_t = 0.80)]
    ensemble_validated_fraction: f64,

    /// Fraction for PROVISIONAL (below VALIDATED).
    #[arg(long, default_value_t = 0.50)]
    ensemble_provisional_fraction: f64,

    /// Fraction of methods whose mean_diff sign must match the
    /// majority sign for either VALIDATED or PROVISIONAL. Default
    /// 1.00 — any directional disagreement drops the grade.
    #[arg(long, default_value_t = 1.00)]
    ensemble_sign_fraction: f64,

    /// Name of the categorical factor (a sample metadata column
    /// referenced by `--design`) whose omnibus F-test is emitted to
    /// `de_omnibus.tsv`. Requires `--test ols` and a `--design` that
    /// includes the factor. At least 3 levels required.
    #[arg(long)]
    omnibus_factor: Option<String>,

    /// Post-hoc contrast adjustment method. `""` (default) disables
    /// post-hoc; `sidak` uses `p_adj = 1 − (1 − p)^m` across the
    /// `--contrast-list`. `tukey` and `dunnett` are rejected until
    /// DEBT-5 and DEBT-6 ship.
    #[arg(long, default_value = "")]
    post_hoc: String,

    /// Comma-separated `level_a-level_b` contrast list applied to
    /// `--post-hoc-factor` (or `--omnibus-factor` as fallback).
    /// Example: `"MCI-CN,AD-CN,AD-MCI"` for `stage` ∈ {CN, MCI, AD}.
    #[arg(long)]
    contrast_list: Option<String>,

    /// Factor whose levels `--contrast-list` references. Falls back
    /// to `--omnibus-factor` when unset.
    #[arg(long)]
    post_hoc_factor: Option<String>,

    /// Family-wise significance threshold for the post-hoc
    /// `decision` column in `de_results.tsv`.
    #[arg(long, default_value_t = 0.05)]
    alpha: f64,
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    let mut args = args;
    apply_per_subject_proxy(&mut args)?;
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
    if !is_unpaired && args.paired_by != "participant" {
        anyhow::bail!(
            "paired-by {:?} not supported for paired tests (expected `participant`)",
            args.paired_by
        );
    }
    if args.min_pairs < 2 {
        anyhow::bail!("min-pairs must be >= 2");
    }
    if args.test == "moderated"
        && (!args.moderation_prior_df.is_finite() || args.moderation_prior_df <= 0.0)
    {
        anyhow::bail!("moderation-prior-df must be > 0 for moderated test");
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
    } else if args.peptide_metadata.is_some()
        && args.test != "limma"
        && args.test != "ensemble"
    {
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
                anyhow::bail!(
                    "--post-hoc sidak requires --post-hoc-factor or --omnibus-factor"
                );
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
                anyhow::bail!(
                    "--post-hoc tukey requires --post-hoc-factor or --omnibus-factor"
                );
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
                anyhow::bail!(
                    "--post-hoc dunnett requires --post-hoc-factor or --omnibus-factor"
                );
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
        other => anyhow::bail!(
            "unknown --post-hoc {other:?}; supported: sidak (tukey / dunnett pending)"
        ),
    }
    if args.test == "ensemble" {
        for frac in [
            args.ensemble_validated_fraction,
            args.ensemble_provisional_fraction,
            args.ensemble_sign_fraction,
        ] {
            if !frac.is_finite() || !(0.0..=1.0).contains(&frac) {
                anyhow::bail!("ensemble fraction thresholds must be in [0, 1]");
            }
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
    let measurements = read_measurements_long(&args.input_dir.join("qc_measurements.tsv"))?;
    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
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

    for m in &measurements {
        let abundance = match m.effective_abundance() {
            Some(a) => a,
            None => continue, // QC-masked
        };
        let s = match sample_by_id.get(m.sample_id.as_str()) {
            Some(s) => *s,
            None => continue,
        };
        if s.is_control {
            continue;
        }
        let subject = match &s.subject_id {
            Some(id) => id.clone(),
            None => continue,
        };
        let condition = match &s.condition {
            Some(c) => c.clone(),
            None => continue,
        };
        let panel = match &m.panel {
            Some(p) => p.clone(),
            None => continue,
        };
        let gene = match &m.gene_symbol {
            Some(g) => g.clone(),
            None => continue,
        };
        if args.test == "ols" || args.test == "mixed" {
            cells_by_sample.insert(
                (panel.clone(), gene.clone(), m.sample_id.clone()),
                abundance,
            );
        }
        cells
            .entry((panel, gene))
            .or_default()
            .entry(condition)
            .or_default()
            .push((subject, abundance));
    }

    let mixed_fixed_formula = args.fixed.as_ref().map(|fixed| format!("~ {fixed}"));
    let model_setup = if args.test == "ols" || args.test == "mixed" {
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

    if args.test == "limma" {
        let (limma_rows, limma_reports) =
            run_limma(&args, &samples, &proteins, &measurements, &comparisons)?;
        all_rows.extend(limma_rows);
        report_rows.extend(limma_reports);
    }

    if args.test == "msqrob" {
        let (msqrob_rows, msqrob_reports) =
            run_msqrob(&args, &samples, &proteins, &comparisons)?;
        all_rows.extend(msqrob_rows);
        report_rows.extend(msqrob_reports);
    }

    if args.test != "limma" && args.test != "msqrob" {
    for (comp_a, comp_b) in &comparisons {
        let comparison_label = format!("{}-{}", comp_a, comp_b);
        // Per-family p-value vector aligned with `family_rows` order.
        let mut family_p: Vec<Option<f64>> = Vec::new();
        let mut family_rows: Vec<DeResultRow> = Vec::new();

        // OLS mode: precompute the encoded design matrix (excluding y) for
        // this comparison so every per-protein fit reuses it, only swapping
        // the y vector (with complete-case filtering on non-finite abundance).
        let model_design: Option<OlsDesign> = if args.test == "ols" || args.test == "mixed" {
            let setup = model_setup
                .as_ref()
                .expect("model_setup populated when test=ols or mixed");
            let design = build_ols_design(&samples, comp_a, comp_b, setup)?;
            design_rows.extend(design.report_rows.clone());
            Some(design)
        } else {
            None
        };

        for ((panel, gene), by_condition) in &cells {
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
            } else if args.test == "ols" {
                let design = model_design
                    .as_ref()
                    .expect("ols_design populated when test=ols");
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
                // paired-t / welch-t.
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
                    mean_a: Some(*mean_a),
                    mean_b: Some(*mean_b),
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
                if q < 0.05 {
                    acc.n_q_lt_05 += 1;
                }
                if q < 0.10 {
                    acc.n_q_lt_10 += 1;
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
                n_q_lt_05: acc.n_q_lt_05,
                n_q_lt_10: acc.n_q_lt_10,
                min_q: acc.min_q,
                max_abs_effect: acc.max_abs_effect,
                limma_trend_fallback_used: None,
            });
        }

        all_rows.extend(family_rows);
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
        write_proxy_summary(&proxy_path, proxy, &covariate_rows)?;
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

    let finished_at = SystemTime::now();
    let input_dir_sha256 = hash_canonical_inputs(
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
            "peptide-measurements": args.peptide_measurements.as_ref().map(|p| p.display().to_string()),
            "peptide-metadata": args.peptide_metadata.as_ref().map(|p| p.display().to_string()),
            "ridge-lambda": args.ridge_lambda,
            "min-peptides": args.min_peptides,
            "omnibus-factor": args.omnibus_factor,
        }),
        &input_dir_sha256,
        &outputs,
        started_at,
        finished_at,
    )?;
    eprintln!("de: sidecar={}", sidecar.display());
    Ok(())
}

fn apply_per_subject_proxy(args: &mut Args) -> Result<()> {
    let Some(proxy) = args
        .per_subject_proxy
        .as_ref()
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
    else {
        return Ok(());
    };
    if args.test == "mixed" {
        anyhow::bail!("--per-subject-proxy currently supports OLS DE; use --test ols");
    }
    if args.fixed.is_some() || args.random.is_some() {
        anyhow::bail!("--per-subject-proxy cannot be combined with --fixed/--random");
    }
    args.test = "ols".to_string();
    let proxy = canonical_proxy_name(proxy);
    if let Some(design) = args.design.as_mut() {
        let terms = parse_design_terms(design)?;
        if !terms.iter().any(|term| match term {
            DesignTerm::Covariate(name) => name == &proxy,
            DesignTerm::Condition => false,
        }) {
            design.push_str(" + ");
            design.push_str(&proxy);
        }
    } else if let Some(covariates) = args.covariates.as_mut() {
        let mut values: Vec<String> = covariates
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect();
        if !values.iter().any(|value| value == &proxy) {
            values.push(proxy.clone());
        }
        *covariates = values.join(",");
    } else {
        args.covariates = Some(proxy.clone());
    }
    args.per_subject_proxy = Some(proxy);
    Ok(())
}

fn canonical_proxy_name(proxy: &str) -> String {
    match proxy.to_ascii_lowercase().as_str() {
        "qalb" | "q_alb" | "q-alb" => "QAlb".to_string(),
        "qigg" | "q_igg" | "q-igg" | "q_{igg}" => "QIgG".to_string(),
        "evans" | "evans_index" | "evans-index" => "Evans_index".to_string(),
        "ventricular_volume" | "ventricular-volume" => "ventricular_volume".to_string(),
        _ => proxy.to_string(),
    }
}

#[derive(Default)]
struct ReportAccumulator {
    n_tests: usize,
    n_skipped: usize,
    n_q_lt_05: usize,
    n_q_lt_10: usize,
    min_q: Option<f64>,
    max_abs_effect: Option<f64>,
}

#[derive(Debug, Clone, Default)]
struct RobustStats {
    effect_size: Option<f64>,
    effect_size_method: String,
    ci_low: Option<f64>,
    ci_high: Option<f64>,
    wilcoxon_p: Option<f64>,
    wilcoxon_method: String,
    median_diff: Option<f64>,
    trimmed_mean_diff: Option<f64>,
}

fn robust_paired(pairs: &[(f64, f64)], min_pairs: usize) -> RobustStats {
    let diffs: Vec<f64> = pairs.iter().map(|(a, b)| a - b).collect();
    if diffs.len() < min_pairs || diffs.len() < 2 || diffs.iter().any(|v| !v.is_finite()) {
        return RobustStats::default();
    }
    let mean_diff = mean(&diffs);
    let sd = sample_sd(&diffs, mean_diff);
    let (ci_low, ci_high) = mean_ci(mean_diff, sd, diffs.len());
    RobustStats {
        effect_size: if sd > 0.0 { Some(mean_diff / sd) } else { None },
        effect_size_method: "cohen_dz".to_string(),
        ci_low,
        ci_high,
        wilcoxon_p: wilcoxon_signed_rank_p(&diffs),
        wilcoxon_method: "signed_rank".to_string(),
        median_diff: median(diffs.clone()),
        trimmed_mean_diff: trimmed_mean(diffs, 0.2),
    }
}

fn robust_unpaired(a: &[f64], b: &[f64], min_pairs: usize) -> RobustStats {
    if a.len() < min_pairs
        || b.len() < min_pairs
        || a.len() < 2
        || b.len() < 2
        || a.iter().chain(b.iter()).any(|v| !v.is_finite())
    {
        return RobustStats::default();
    }
    let mean_a = mean(a);
    let mean_b = mean(b);
    let var_a = sample_var(a, mean_a);
    let var_b = sample_var(b, mean_b);
    let mean_diff = mean_a - mean_b;
    let df_pooled = (a.len() + b.len() - 2) as f64;
    let pooled_var = ((a.len() - 1) as f64 * var_a + (b.len() - 1) as f64 * var_b) / df_pooled;
    let hedges_j = 1.0 - 3.0 / (4.0 * df_pooled - 1.0);
    let effect_size = if pooled_var > 0.0 {
        Some(hedges_j * mean_diff / pooled_var.sqrt())
    } else {
        None
    };
    let se = (var_a / a.len() as f64 + var_b / b.len() as f64).sqrt();
    let df = welch_df(var_a, var_b, a.len(), b.len());
    let (ci_low, ci_high) = if se > 0.0 && df.is_finite() && df > 0.0 {
        let dist = StudentsT::new(0.0, 1.0, df).ok();
        let crit = dist.map(|d| d.inverse_cdf(0.975));
        (
            crit.map(|c| mean_diff - c * se),
            crit.map(|c| mean_diff + c * se),
        )
    } else {
        (None, None)
    };
    RobustStats {
        effect_size,
        effect_size_method: "hedges_g".to_string(),
        ci_low,
        ci_high,
        wilcoxon_p: wilcoxon_rank_sum_p(a, b),
        wilcoxon_method: "rank_sum".to_string(),
        median_diff: median(a.to_vec())
            .zip(median(b.to_vec()))
            .map(|(ma, mb)| ma - mb),
        trimmed_mean_diff: trimmed_mean(a.to_vec(), 0.2)
            .zip(trimmed_mean(b.to_vec(), 0.2))
            .map(|(ma, mb)| ma - mb),
    }
}


fn sample_var(values: &[f64], mean: f64) -> f64 {
    values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (values.len() - 1) as f64
}

fn sample_sd(values: &[f64], mean: f64) -> f64 {
    sample_var(values, mean).sqrt()
}

fn mean_ci(mean: f64, sd: f64, n: usize) -> (Option<f64>, Option<f64>) {
    if n < 2 || sd == 0.0 {
        return (None, None);
    }
    let df = (n - 1) as f64;
    let Some(crit) = StudentsT::new(0.0, 1.0, df)
        .ok()
        .map(|d| d.inverse_cdf(0.975))
    else {
        return (None, None);
    };
    let se = sd / (n as f64).sqrt();
    (Some(mean - crit * se), Some(mean + crit * se))
}

fn welch_df(var_a: f64, var_b: f64, n_a: usize, n_b: usize) -> f64 {
    let va = var_a / n_a as f64;
    let vb = var_b / n_b as f64;
    let den = va.powi(2) / (n_a - 1) as f64 + vb.powi(2) / (n_b - 1) as f64;
    if den == 0.0 {
        f64::NAN
    } else {
        (va + vb).powi(2) / den
    }
}

fn median(mut values: Vec<f64>) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = values.len() / 2;
    if values.len() % 2 == 1 {
        Some(values[mid])
    } else {
        Some((values[mid - 1] + values[mid]) / 2.0)
    }
}

fn trimmed_mean(mut values: Vec<f64>, proportion: f64) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let trim = ((values.len() as f64) * proportion).floor() as usize;
    if trim * 2 >= values.len() {
        return None;
    }
    Some(mean(&values[trim..(values.len() - trim)]))
}

fn wilcoxon_signed_rank_p(diffs: &[f64]) -> Option<f64> {
    let nonzero: Vec<f64> = diffs.iter().copied().filter(|d| *d != 0.0).collect();
    let n = nonzero.len();
    if n < 2 {
        return None;
    }
    let abs_values: Vec<f64> = nonzero.iter().map(|d| d.abs()).collect();
    let ranks = average_ranks(&abs_values);
    let w_plus: f64 = nonzero
        .iter()
        .zip(ranks.iter())
        .filter_map(|(d, r)| if *d > 0.0 { Some(*r) } else { None })
        .sum();
    let n_f = n as f64;
    let mean_w = n_f * (n_f + 1.0) / 4.0;
    let var_w = n_f * (n_f + 1.0) * (2.0 * n_f + 1.0) / 24.0;
    normal_two_sided_p(w_plus, mean_w, var_w)
}

fn wilcoxon_rank_sum_p(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.len() < 2 || b.len() < 2 {
        return None;
    }
    let mut values = Vec::with_capacity(a.len() + b.len());
    values.extend_from_slice(a);
    values.extend_from_slice(b);
    let ranks = average_ranks(&values);
    let rank_sum_a: f64 = ranks.iter().take(a.len()).sum();
    let n_a = a.len() as f64;
    let n_b = b.len() as f64;
    let u_a = rank_sum_a - n_a * (n_a + 1.0) / 2.0;
    let mean_u = n_a * n_b / 2.0;
    let var_u = n_a * n_b * (n_a + n_b + 1.0) / 12.0;
    normal_two_sided_p(u_a, mean_u, var_u)
}

fn normal_two_sided_p(stat: f64, mean: f64, var: f64) -> Option<f64> {
    if !(var.is_finite() && var > 0.0) {
        return None;
    }
    let z = (stat - mean).abs() / var.sqrt();
    let dist = Normal::new(0.0, 1.0).ok()?;
    Some(2.0 * (1.0 - dist.cdf(z)))
}

fn average_ranks(values: &[f64]) -> Vec<f64> {
    let mut indexed: Vec<(usize, f64)> = values.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let mut ranks = vec![0.0; values.len()];
    let mut i = 0;
    while i < indexed.len() {
        let mut j = i + 1;
        while j < indexed.len() && indexed[j].1 == indexed[i].1 {
            j += 1;
        }
        let avg_rank = (i + 1 + j) as f64 / 2.0;
        for k in i..j {
            ranks[indexed[k].0] = avg_rank;
        }
        i = j;
    }
    ranks
}

fn apply_moderated_shrinkage(rows: &mut [DeResultRow], prior_df: f64) -> Result<()> {
    let mut vars: Vec<f64> = rows
        .iter()
        .filter_map(paired_variance_from_row)
        .filter(|v| v.is_finite() && *v > 0.0)
        .collect();
    if vars.is_empty() {
        return Ok(());
    }
    vars.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let prior_var = if vars.len() % 2 == 1 {
        vars[vars.len() / 2]
    } else {
        (vars[vars.len() / 2 - 1] + vars[vars.len() / 2]) / 2.0
    };

    for row in rows.iter_mut() {
        let mean_diff = match row.mean_diff {
            Some(v) => v,
            None => continue,
        };
        let df_i = match row.df {
            Some(v) if v > 0.0 => v,
            _ => continue,
        };
        let n = row.n_pairs as f64;
        if n < 2.0 {
            continue;
        }
        let var_i = match paired_variance_from_row(row) {
            Some(v) if v.is_finite() && v >= 0.0 => v,
            _ => continue,
        };
        let post_var = (prior_df * prior_var + df_i * var_i) / (prior_df + df_i);
        if !(post_var.is_finite() && post_var > 0.0) {
            continue;
        }
        let t_mod = mean_diff / (post_var / n).sqrt();
        let df_mod = prior_df + df_i;
        let dist = StudentsT::new(0.0, 1.0, df_mod)
            .with_context(|| format!("building t distribution with df={}", df_mod))?;
        let p_mod = 2.0 * (1.0 - dist.cdf(t_mod.abs()));
        if p_mod.is_finite() {
            row.t = Some(t_mod);
            row.df = Some(df_mod);
            row.p_value = Some(p_mod);
        }
    }
    Ok(())
}

fn paired_variance_from_row(row: &DeResultRow) -> Option<f64> {
    let t = row.t?;
    let mean = row.mean_diff?;
    let n = row.n_pairs as f64;
    if n < 2.0 {
        return None;
    }
    if t == 0.0 {
        return Some(0.0);
    }
    let se = mean / t;
    Some((se * n.sqrt()).powi(2))
}

// -----------------------------------------------------------------
// OLS (covariate-adjusted) helpers
// -----------------------------------------------------------------

#[derive(Debug, Clone)]
enum DesignTerm {
    Condition,
    Covariate(String),
}

struct OlsSetup {
    terms: Vec<DesignTerm>,
    cov_names: Vec<String>,
    cov_raw: HashMap<String, Vec<Option<String>>>,
    contrast: Option<String>,
    label: String,
}

fn build_ols_setup(
    samples_path: &Path,
    design: Option<&str>,
    covariates: Option<&str>,
    contrast: Option<&str>,
) -> Result<OlsSetup> {
    if let Some(formula) = design {
        let terms = parse_design_terms(formula)?;
        if !terms.iter().any(|t| matches!(t, DesignTerm::Condition)) {
            anyhow::bail!("--design must include `condition` for DE contrasts");
        }
        let cov_names = covariate_names_from_terms(&terms);
        let cov_raw = if cov_names.is_empty() {
            HashMap::new()
        } else {
            read_covariate_columns(samples_path, &cov_names)?
        };
        Ok(OlsSetup {
            terms,
            cov_names,
            cov_raw,
            contrast: contrast.map(str::to_string),
            label: formula.to_string(),
        })
    } else {
        let cov_names: Vec<String> = covariates
            .unwrap_or("")
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect();
        let cov_raw = if cov_names.is_empty() {
            HashMap::new()
        } else {
            read_covariate_columns(samples_path, &cov_names)?
        };
        let mut terms = vec![DesignTerm::Condition];
        terms.extend(cov_names.iter().cloned().map(DesignTerm::Covariate));
        Ok(OlsSetup {
            terms,
            cov_names,
            cov_raw,
            contrast: contrast.map(str::to_string),
            label: if covariates.unwrap_or("").trim().is_empty() {
                "~ condition".to_string()
            } else {
                format!("~ condition + {}", covariates.unwrap_or("").trim())
            },
        })
    }
}

fn parse_design_terms(formula: &str) -> Result<Vec<DesignTerm>> {
    let formula = formula.trim();
    let rhs = formula
        .strip_prefix('~')
        .ok_or_else(|| anyhow!("--design must start with `~`"))?
        .trim();
    if rhs.is_empty() {
        anyhow::bail!("--design must contain at least `condition`");
    }
    let mut terms = Vec::new();
    let mut seen = BTreeSet::new();
    for raw in rhs.split('+') {
        let term = raw.trim();
        if term.is_empty() {
            anyhow::bail!("empty term in --design");
        }
        if term == "1" {
            continue;
        }
        if term == "0" || term == "-1" {
            anyhow::bail!("intercept removal is not supported; Atman includes an intercept");
        }
        if !seen.insert(term.to_string()) {
            anyhow::bail!("duplicate term {:?} in --design", term);
        }
        if term == "condition" {
            terms.push(DesignTerm::Condition);
        } else {
            terms.push(DesignTerm::Covariate(term.to_string()));
        }
    }
    if terms.is_empty() {
        anyhow::bail!("--design must contain at least `condition`");
    }
    Ok(terms)
}

fn covariate_names_from_terms(terms: &[DesignTerm]) -> Vec<String> {
    terms
        .iter()
        .filter_map(|t| match t {
            DesignTerm::Condition => None,
            DesignTerm::Covariate(name) => Some(name.clone()),
        })
        .collect()
}

fn resolve_comparisons(
    groups: Option<&str>,
    samples: &[Sample],
    test: &str,
    design: Option<&str>,
    contrast: Option<&str>,
) -> Result<Vec<(String, String)>> {
    if let Some(groups) = groups {
        return parse_comparisons(groups);
    }
    if (test != "ols" && test != "mixed") || design.is_none() {
        anyhow::bail!("--groups is required unless using --test ols --design with --contrast");
    }
    let contrast =
        contrast.ok_or_else(|| anyhow!("--contrast is required when --groups is omitted"))?;
    let comp_a = contrast
        .strip_prefix("condition")
        .filter(|s| !s.is_empty())
        .ok_or_else(|| {
            anyhow!("can only infer --groups from condition contrasts like conditionCase")
        })?;
    let conditions: BTreeSet<String> = samples
        .iter()
        .filter(|s| !s.is_control)
        .filter_map(|s| s.condition.clone())
        .collect();
    if !conditions.contains(comp_a) {
        anyhow::bail!(
            "contrast condition {:?} is not present in samples.tsv",
            comp_a
        );
    }
    let others: Vec<String> = conditions.into_iter().filter(|c| c != comp_a).collect();
    if others.len() != 1 {
        anyhow::bail!(
            "cannot infer --groups for contrast {}; provide --groups A-B when more than two conditions are present",
            contrast
        );
    }
    Ok(vec![(comp_a.to_string(), others[0].clone())])
}

fn validate_random_intercept_design(
    samples: &[Sample],
    comparisons: &[(String, String)],
    min_samples: usize,
) -> Result<()> {
    for (a, b) in comparisons {
        let mut by_subject: BTreeMap<String, usize> = BTreeMap::new();
        for sample in samples {
            if sample.is_control {
                continue;
            }
            let Some(condition) = &sample.condition else {
                continue;
            };
            if condition != a && condition != b {
                continue;
            }
            let subject = sample.subject_id.as_ref().unwrap_or(&sample.sample_id);
            *by_subject.entry(subject.clone()).or_default() += 1;
        }
        let repeated_subjects = by_subject.values().filter(|n| **n >= 2).count();
        let n_samples: usize = by_subject.values().sum();
        if n_samples < min_samples || repeated_subjects < 2 {
            anyhow::bail!(
                "mixed model for {a}-{b} requires at least {} samples and at least two repeated subjects for --random '1|subject_id'",
                min_samples
            );
        }
    }
    Ok(())
}

/// Covariate kind determined at column-classification time. `Numeric`
/// columns are passed through as f64. `Categorical` columns are one-hot
/// encoded using the supplied levels; the first element of `levels` is the
/// reference and is dropped.
#[derive(Debug, Clone)]
enum CovKind {
    Numeric,
    Categorical { levels: Vec<String> },
}

impl CovKind {
    /// Human label for each design column this covariate contributes.
    /// Used for the header of de_covariates.tsv.
    fn col_labels(&self, base: &str) -> Vec<String> {
        match self {
            CovKind::Numeric => vec![base.to_string()],
            CovKind::Categorical { levels } => levels
                .iter()
                .skip(1)
                .map(|lvl| format!("{}{}", base, lvl))
                .collect(),
        }
    }
}

/// Encoded design matrix for one comparison. `rows` lists every sample
/// that (a) is in `comp_a` or `comp_b`, (b) is non-control, and (c) has a
/// non-empty value for every covariate; each entry is
/// `(sample_id, [intercept, group, cov_col_0, cov_col_1, …])`. `group_col`
/// is the index of the group column (always 1 under the current design)
/// so the caller knows which β to extract as the effect of interest.
/// `design_labels[i]` names the i-th design column for de_covariates.tsv.
struct OlsDesign {
    rows: Vec<(String, Vec<f64>)>,
    group_col: usize,
    design_labels: Vec<String>,
    report_rows: Vec<DesignReportRow>,
}

/// Row shape for de_omnibus.tsv — per (comparison × protein) omnibus
/// F-test of the user-specified factor's columns.
struct OmnibusRow {
    panel: String,
    assay_id: String,
    gene_symbol: String,
    factor: String,
    comparison: String,
    f_statistic: f64,
    df_num: usize,
    df_den: f64,
    p_value: f64,
    bh_q: Option<f64>,
}

fn write_omnibus_rows(path: &Path, rows: &[OmnibusRow]) -> Result<()> {
    let mut buf = String::from(
        "panel\tassay_id\tgene_symbol\tfactor\tcomparison\tf_statistic\tdf_num\tdf_den\tp_value\tbh_q\n",
    );
    for r in rows {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.panel,
            r.assay_id,
            r.gene_symbol,
            r.factor,
            r.comparison,
            format_opt(r.f_statistic),
            r.df_num,
            format_opt(r.df_den),
            format_opt(r.p_value),
            r.bh_q.map(format_opt).unwrap_or_default(),
        ));
    }
    atomic_write(path, buf.as_bytes())
}

fn format_opt(v: f64) -> String {
    if !v.is_finite() {
        if v.is_nan() {
            "NaN".into()
        } else if v > 0.0 {
            "Inf".into()
        } else {
            "-Inf".into()
        }
    } else {
        format!("{v}")
    }
}

/// Row shape for de_covariates.tsv — per (comparison × protein × covariate
/// column) estimate, skipping the intercept and the group indicator.
struct CovariateRow {
    comparison: String,
    panel: String,
    gene_symbol: String,
    covariate: String,
    beta: f64,
    se: f64,
    t: f64,
    p_value: f64,
    df: f64,
    n: usize,
}

#[derive(Clone)]
struct DesignReportRow {
    comparison: String,
    sample_id: String,
    condition: String,
    included: bool,
    drop_reason: String,
    columns: String,
    values: String,
}

/// Raw-read `samples.tsv` and return `sample_id → [cov_value]` for the
/// requested covariate names. Values that are missing (empty cell) are
/// returned as `None`; the caller drops them before fitting.
fn read_covariate_columns(
    path: &Path,
    names: &[String],
) -> Result<HashMap<String, Vec<Option<String>>>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let sample_col = headers
        .iter()
        .position(|h| h == "sample_id")
        .ok_or_else(|| anyhow!("samples.tsv missing 'sample_id' column"))?;
    let cov_cols: Vec<usize> = names
        .iter()
        .map(|n| {
            headers
                .iter()
                .position(|h| h == n.as_str())
                .ok_or_else(|| anyhow!("covariate {:?} not found in samples.tsv", n))
        })
        .collect::<Result<_>>()?;
    let mut out: HashMap<String, Vec<Option<String>>> = HashMap::new();
    for row in reader.records() {
        let row = row?;
        let sid = row
            .get(sample_col)
            .ok_or_else(|| anyhow!("samples.tsv row missing sample_id"))?
            .to_string();
        let vals: Vec<Option<String>> = cov_cols
            .iter()
            .map(|&i| {
                let v = row.get(i).unwrap_or("");
                if v.is_empty() {
                    None
                } else {
                    Some(v.to_string())
                }
            })
            .collect();
        out.insert(sid, vals);
    }
    Ok(out)
}

/// Classify each covariate as Numeric or Categorical based on whether
/// every non-missing value parses as f64.
fn classify_covariates(
    names: &[String],
    values_per_sample: &[Vec<Option<String>>],
) -> Vec<CovKind> {
    let mut out = Vec::with_capacity(names.len());
    for (i, _) in names.iter().enumerate() {
        let mut all_numeric = true;
        let mut levels: BTreeSet<String> = BTreeSet::new();
        for row in values_per_sample {
            if let Some(v) = &row[i] {
                if v.parse::<f64>().is_err() {
                    all_numeric = false;
                }
                levels.insert(v.clone());
            }
        }
        if all_numeric {
            out.push(CovKind::Numeric);
        } else {
            out.push(CovKind::Categorical {
                levels: levels.into_iter().collect(),
            });
        }
    }
    out
}

/// Build the design matrix for one comparison. Samples not in `comp_a` or
/// `comp_b`, control samples, and samples with any missing covariate are
/// dropped before encoding.
fn build_ols_design(
    samples: &[Sample],
    comp_a: &str,
    comp_b: &str,
    setup: &OlsSetup,
) -> Result<OlsDesign> {
    let comparison = format!("{comp_a}-{comp_b}");
    let condition_label = format!("condition{comp_a}");

    // Pick samples in either group.
    let in_group: Vec<&Sample> = samples
        .iter()
        .filter(|s| {
            !s.is_control
                && match &s.condition {
                    Some(c) => c == comp_a || c == comp_b,
                    None => false,
                }
        })
        .collect();

    // Pull covariate values aligned with `in_group`; drop samples missing
    // any formula covariate before classifying/encoding.
    let mut kept: Vec<(&Sample, Vec<Option<String>>)> = Vec::new();
    let mut report_rows = Vec::new();
    for s in &in_group {
        let vals = if setup.cov_names.is_empty() {
            Vec::new()
        } else {
            match setup.cov_raw.get(&s.sample_id) {
                Some(v) => v.clone(),
                None => {
                    report_rows.push(DesignReportRow {
                        comparison: comparison.clone(),
                        sample_id: s.sample_id.clone(),
                        condition: s.condition.clone().unwrap_or_default(),
                        included: false,
                        drop_reason: "missing_covariate_row".to_string(),
                        columns: String::new(),
                        values: String::new(),
                    });
                    continue;
                }
            }
        };
        if vals.iter().any(|v| v.is_none()) {
            report_rows.push(DesignReportRow {
                comparison: comparison.clone(),
                sample_id: s.sample_id.clone(),
                condition: s.condition.clone().unwrap_or_default(),
                included: false,
                drop_reason: "missing_covariate".to_string(),
                columns: String::new(),
                values: String::new(),
            });
            continue;
        }
        kept.push((s, vals));
    }

    // Classify covariate kinds from the retained samples.
    let cov_kinds = classify_covariates(
        &setup.cov_names,
        &kept.iter().map(|(_, v)| v.clone()).collect::<Vec<_>>(),
    );

    let mut design_labels = vec!["(Intercept)".to_string()];
    for term in &setup.terms {
        match term {
            DesignTerm::Condition => design_labels.push(condition_label.clone()),
            DesignTerm::Covariate(name) => {
                let idx = cov_index(&setup.cov_names, name)?;
                design_labels.extend(cov_kinds[idx].col_labels(name));
            }
        }
    }
    let contrast = setup
        .contrast
        .as_deref()
        .unwrap_or(condition_label.as_str());
    let group_col = design_labels
        .iter()
        .position(|label| label == contrast)
        .ok_or_else(|| {
            anyhow!(
                "contrast {:?} not found in design columns: {}",
                contrast,
                design_labels.join(",")
            )
        })?;

    let mut rows: Vec<(String, Vec<f64>)> = Vec::with_capacity(kept.len());
    for (s, vals) in kept {
        let mut row: Vec<f64> = Vec::with_capacity(design_labels.len());
        row.push(1.0);
        for term in &setup.terms {
            match term {
                DesignTerm::Condition => {
                    row.push(if s.condition.as_deref() == Some(comp_a) {
                        1.0
                    } else {
                        0.0
                    });
                }
                DesignTerm::Covariate(name) => {
                    let idx = cov_index(&setup.cov_names, name)?;
                    push_encoded_covariate(&mut row, name, &cov_kinds[idx], &vals[idx])?;
                }
            }
        }
        report_rows.push(DesignReportRow {
            comparison: comparison.clone(),
            sample_id: s.sample_id.clone(),
            condition: s.condition.clone().unwrap_or_default(),
            included: true,
            drop_reason: String::new(),
            columns: design_labels.join(","),
            values: row
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(","),
        });
        rows.push((s.sample_id.clone(), row));
    }
    ensure_full_rank(&rows, &design_labels, &comparison)?;

    Ok(OlsDesign {
        rows,
        group_col,
        design_labels,
        report_rows,
    })
}

fn cov_index(names: &[String], name: &str) -> Result<usize> {
    names
        .iter()
        .position(|n| n == name)
        .ok_or_else(|| anyhow!("internal error: missing covariate {:?}", name))
}

fn push_encoded_covariate(
    row: &mut Vec<f64>,
    name: &str,
    kind: &CovKind,
    value: &Option<String>,
) -> Result<()> {
    let v = value
        .as_ref()
        .expect("missing covariate filtered out above");
    match kind {
        CovKind::Numeric => {
            let parsed = v
                .parse::<f64>()
                .with_context(|| format!("numeric covariate {:?} contains {:?}", name, v))?;
            row.push(parsed);
        }
        CovKind::Categorical { levels } => {
            for lvl in levels.iter().skip(1) {
                row.push(if v == lvl { 1.0 } else { 0.0 });
            }
        }
    }
    Ok(())
}

fn ensure_full_rank(
    rows: &[(String, Vec<f64>)],
    labels: &[String],
    comparison: &str,
) -> Result<()> {
    let rank = matrix_rank(&rows.iter().map(|(_, row)| row.clone()).collect::<Vec<_>>());
    if rank < labels.len() {
        anyhow::bail!(
            "singular design matrix for {comparison}: rank {} < {} columns ({})",
            rank,
            labels.len(),
            labels.join(",")
        );
    }
    Ok(())
}

fn matrix_rank(rows: &[Vec<f64>]) -> usize {
    if rows.is_empty() {
        return 0;
    }
    let mut a = rows.to_vec();
    let n_rows = a.len();
    let n_cols = a[0].len();
    let mut rank = 0;
    let eps = 1e-10;
    for col in 0..n_cols {
        let pivot = (rank..n_rows).max_by(|&i, &j| {
            a[i][col]
                .abs()
                .partial_cmp(&a[j][col].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let Some(pivot) = pivot else {
            continue;
        };
        if a[pivot][col].abs() <= eps {
            continue;
        }
        a.swap(rank, pivot);
        let pivot_value = a[rank][col];
        for value in &mut a[rank][col..n_cols] {
            *value /= pivot_value;
        }
        let pivot_tail = a[rank][col..n_cols].to_vec();
        for (r, row) in a.iter_mut().enumerate().take(n_rows) {
            if r == rank {
                continue;
            }
            let factor = row[col];
            for (value, pivot_value) in row[col..n_cols].iter_mut().zip(pivot_tail.iter()) {
                *value -= factor * pivot_value;
            }
        }
        rank += 1;
        if rank == n_rows {
            break;
        }
    }
    rank
}

/// Raw per-group means of the response within the fitted rows. Used to
/// populate `mean_a` / `mean_b` on the `DeResultRow` for schema
/// compatibility with paired-t / welch-t output.
fn raw_group_means(y: &[f64], rows: &[Vec<f64>], group_col: usize) -> (Option<f64>, Option<f64>) {
    let mut sum_a = 0.0;
    let mut n_a = 0;
    let mut sum_b = 0.0;
    let mut n_b = 0;
    for (yi, row) in y.iter().zip(rows.iter()) {
        if row[group_col] > 0.5 {
            sum_a += yi;
            n_a += 1;
        } else {
            sum_b += yi;
            n_b += 1;
        }
    }
    let ma = if n_a > 0 {
        Some(sum_a / n_a as f64)
    } else {
        None
    };
    let mb = if n_b > 0 {
        Some(sum_b / n_b as f64)
    } else {
        None
    };
    (ma, mb)
}

/// Convert an `OlsOutcome` into the `PairedTResult` enum used by the
/// downstream writer. The group coefficient becomes `mean_diff`; residual
/// df becomes `df`; raw group means (supplied) become `mean_a` / `mean_b`.
/// Non-intercept, non-group coefficients are recorded into
/// `covariate_rows` for emission to de_covariates.tsv.
#[allow(clippy::too_many_arguments)]
fn ols_to_paired_t_result(
    outcome: OlsOutcome,
    design: &OlsDesign,
    mean_a: Option<f64>,
    mean_b: Option<f64>,
    panel: &str,
    gene: &str,
    comparison: &str,
    covariate_rows: &mut Vec<CovariateRow>,
) -> PairedTResult {
    match outcome {
        OlsOutcome::Computed(fit) => {
            let g = design.group_col;
            // Record per-covariate estimates (everything except intercept and group).
            for (col, label) in design.design_labels.iter().enumerate() {
                if col == 0 || col == g {
                    continue;
                }
                covariate_rows.push(CovariateRow {
                    comparison: comparison.to_string(),
                    panel: panel.to_string(),
                    gene_symbol: gene.to_string(),
                    covariate: label.clone(),
                    beta: fit.beta[col],
                    se: fit.se[col],
                    t: fit.t[col],
                    p_value: fit.p_value[col],
                    df: fit.df,
                    n: fit.n,
                });
            }
            PairedTResult::Computed {
                n_pairs: fit.n,
                mean_a: mean_a.unwrap_or(f64::NAN),
                mean_b: mean_b.unwrap_or(f64::NAN),
                mean_diff: fit.beta[g],
                t: fit.t[g],
                df: fit.df,
                p_value: fit.p_value[g],
            }
        }
        OlsOutcome::Skipped { reason, n } => PairedTResult::Skipped { reason, n_pairs: n },
    }
}

/// Writer for `de_covariates.tsv`. Tab-separated, one row per
/// (comparison × panel × gene × covariate column).
fn write_covariate_rows(path: &Path, rows: &[CovariateRow]) -> Result<()> {
    let mut buf =
        String::from("comparison\tpanel\tgene_symbol\tcovariate\tbeta\tse\tt\tp_value\tdf\tn\n");
    for r in rows {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.comparison,
            r.panel,
            r.gene_symbol,
            r.covariate,
            r.beta,
            r.se,
            r.t,
            r.p_value,
            r.df,
            r.n,
        ));
    }
    atomic_write(path, buf.as_bytes())
}

fn write_proxy_summary(path: &Path, proxy: &str, rows: &[CovariateRow]) -> Result<()> {
    let mut by_comparison: BTreeMap<&str, Vec<&CovariateRow>> = BTreeMap::new();
    for row in rows {
        if row.covariate == proxy {
            by_comparison
                .entry(row.comparison.as_str())
                .or_default()
                .push(row);
        }
    }
    let mut buf = String::from(
        "proxy\tcomparison\tn_tests\tn_p_lt_05\tmedian_abs_beta\tmedian_p_value\tinterpretation\n",
    );
    for (comparison, rows) in by_comparison {
        let n_tests = rows.len();
        let n_p_lt_05 = rows.iter().filter(|row| row.p_value < 0.05).count();
        let median_abs_beta = median(rows.iter().map(|row| row.beta.abs()).collect());
        let median_p_value = median(rows.iter().map(|row| row.p_value).collect());
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            proxy,
            comparison,
            n_tests,
            n_p_lt_05,
            fmt_opt(median_abs_beta),
            fmt_opt(median_p_value),
            proxy_interpretation(proxy)
        ));
    }
    atomic_write(path, buf.as_bytes())
}

fn proxy_interpretation(proxy: &str) -> &'static str {
    match proxy {
        "QAlb" => "albumin quotient barrier-clearance adjustment",
        "QIgG" => "immunoglobulin quotient barrier-clearance adjustment",
        "Evans_index" => "ventricular size barrier-clearance adjustment",
        "ventricular_volume" => "ventricular volume barrier-clearance adjustment",
        _ => "physiological proxy adjustment",
    }
}

fn fmt_opt(value: Option<f64>) -> String {
    value
        .filter(|v| v.is_finite())
        .map(|v| format!("{v:.6}"))
        .unwrap_or_else(|| "NA".to_string())
}

/// Build the two-group design matrix for a limma `(a, b)` comparison.
///
/// Rows correspond to non-control samples whose `condition` is `a` or `b`,
/// in the order they appear in the input `samples` slice. Columns are
/// `[intercept, group_b_indicator]` — `group_b_indicator = 1.0` when the
/// sample's condition is `b`, `0.0` when it is `a`. The returned
/// `sample_ids` vector is in the same row order and is used to align
/// abundance lookups in the `y` matrix.
fn build_two_group_design(
    samples: &[Sample],
    a: &str,
    b: &str,
) -> (Vec<Vec<f64>>, Vec<String>) {
    let mut design = Vec::new();
    let mut sample_ids = Vec::new();
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
        design.push(vec![1.0, group_b]);
        sample_ids.push(s.sample_id.clone());
    }
    (design, sample_ids)
}

/// Dispatch for `--test limma`. Builds a two-group design per comparison,
/// constructs a features × samples `y` matrix (missing cells filled with
/// NaN so limma's per-feature ols skips them), calls
/// [`atman_core::limma::limma_fit`], then maps each `LimmaRow` into a
/// `DeResultRow` and emits one `DeReportRow` per (comparison, panel)
/// family. BH-FDR is applied per panel within each comparison.
fn run_limma(
    args: &Args,
    samples: &[Sample],
    proteins: &[atman_core::ProteinIdentity],
    measurements: &[atman_core::MeasurementRecord],
    comparisons: &[(String, String)],
) -> Result<(Vec<DeResultRow>, Vec<DeReportRow>)> {
    use atman_core::limma::{limma_fit, LimmaOptions};
    use crate::io::read_peptides;

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
        let (design, sample_ids) = build_two_group_design(samples, a, b);
        if design.is_empty() {
            // No samples in either group — skip this comparison silently;
            // no rows, no report.
            continue;
        }

        // Feature set = (panel, gene_symbol) pairs observed in
        // measurements that have at least one retained-sample value in
        // this comparison. `gene_meta` lookups fall back to empty
        // (assay_id, uniprot) when a measured gene isn't listed in
        // proteins.tsv — matching the welch-t / ols behaviour. Iteration
        // is deterministic via BTreeSet.
        let mut features: Vec<(String, String, String, Vec<String>)> = Vec::new();
        for (panel, gene) in &measured_features {
            let has_any = sample_ids.iter().any(|sid| {
                cells_by_sample
                    .contains_key(&(panel.clone(), gene.clone(), sid.clone()))
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
                n_q_lt_05: 0,
                n_q_lt_10: 0,
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
        let contrast_matrix: Vec<Vec<f64>> = vec![vec![0.0], vec![-1.0]];

        // (iv) Limma fit. When --peptide-metadata is supplied on the
        // limma path, count peptides per feature and switch the trend
        // covariate from mean-log2-abundance to `log(peptide_count +
        // 1)` (DEqMS, Zhu et al. 2020). The peptide→parent mapping
        // comes from the ingested `peptides.tsv`; features without
        // peptide records get count 0.
        let peptide_counts: Option<Vec<u32>> = if let Some(counts_map) = peptides_per_assay.as_ref() {
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
                        n_q_lt_05: 0,
                        n_q_lt_10: 0,
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
        let panel_by_feature: Vec<String> =
            features.iter().map(|(p, _, _, _)| p.clone()).collect();
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
        let t_crit = StudentsT::new(0.0, 1.0, output.df_total)
            .ok()
            .map(|d| d.inverse_cdf(0.975));

        #[derive(Default)]
        struct PanelAcc {
            n_tests: usize,
            n_skipped: usize,
            n_q_lt_05: usize,
            n_q_lt_10: usize,
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
                    Some(c) if s.is_finite() && c.is_finite() => {
                        (Some(e - c * s), Some(e + c * s))
                    }
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
                    if qv < 0.05 {
                        acc.n_q_lt_05 += 1;
                    }
                    if qv < 0.10 {
                        acc.n_q_lt_10 += 1;
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
                n_q_lt_05: acc.n_q_lt_05,
                n_q_lt_10: acc.n_q_lt_10,
                min_q: acc.min_q,
                max_abs_effect: acc.max_abs_effect,
                limma_trend_fallback_used: Some(output.trend_fallback_used),
            });
        }
    }

    Ok((all_result_rows, all_report_rows))
}

/// Dispatch for `--test msqrob`. One ridge-regularized linear mixed
/// model per protein, fit over its peptide-level observations in the
/// two-group comparison. The design matrix is `[intercept, group_b]`;
/// the reported effect is `mean_a − mean_b`, matching the sign
/// convention used by the paired-t, welch-t, and ols dispatches
/// (contrast = −β_group_b).
fn run_msqrob(
    args: &Args,
    samples: &[Sample],
    proteins: &[atman_core::ProteinIdentity],
    comparisons: &[(String, String)],
) -> Result<(Vec<DeResultRow>, Vec<DeReportRow>)> {
    use atman_core::msqrob::{fit_msqrob, squeeze_variance, MsqrobFit, MsqrobOutcome};
    use crate::io::{read_peptide_measurements, read_peptides};

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

    let ridge_lambda: f64 = match args.ridge_lambda.trim() {
        "auto" => 0.0,
        other => other
            .parse::<f64>()
            .with_context(|| format!("--ridge-lambda must be a number or `auto`, got {other:?}"))?,
    };
    if !ridge_lambda.is_finite() || ridge_lambda < 0.0 {
        anyhow::bail!("--ridge-lambda must be finite and >= 0 (got {ridge_lambda})");
    }

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
                if sample_by_id.get(sample_id.as_str()).is_none() {
                    continue;
                }
                for pid in peptide_ids {
                    let key = (sample_id.clone(), pid.clone());
                    if let Some(&abund) = pep_cells.get(&key) {
                        design.push(vec![1.0, *group_b]);
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
            let ci_low = effect - 1.96 * se;
            let ci_high = effect + 1.96 * se;
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
        for (_panel, indices) in &by_panel {
            let ps: Vec<Option<f64>> = indices.iter().map(|&i| fits[i].1.de.p_value).collect();
            let qs = atman_core::bh_fdr(&ps);
            for (local, &global) in indices.iter().enumerate() {
                fits[global].1.de.bh_q = qs[local];
                if let Some(q) = qs[local] {
                    let panel = fits[global].0.clone();
                    let acc = per_panel.entry(panel).or_default();
                    if q < 0.05 {
                        acc.n_q_lt_05 += 1;
                    }
                    if q < 0.10 {
                        acc.n_q_lt_10 += 1;
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
                n_q_lt_05: acc.n_q_lt_05,
                n_q_lt_10: acc.n_q_lt_10,
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
    n_q_lt_05: usize,
    n_q_lt_10: usize,
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
    let mean_a = if n_a > 0 { sum_a / n_a as f64 } else { f64::NAN };
    let mean_b = if n_b > 0 { sum_b / n_b as f64 } else { f64::NAN };
    (mean_a, mean_b)
}

/// Dispatch for `--test ensemble`: fan out to every requested method
/// in a tempdir, then aggregate per (comparison, protein) into a
/// Stouffer-combined ensemble p / q and a VALIDATED / PROVISIONAL /
/// INSUFFICIENT grade. Each sub-method runs via a cloned `Args` with
/// conflicting flags stripped; errors are collected into
/// `methods_skipped` rather than aborting the whole run.
fn run_ensemble(args: Args, started_at: SystemTime) -> Result<()> {
    use atman_core::ensemble::{
        aggregate_per_protein, assign_grade, EnsembleInput, GradeThresholds,
    };
    use crate::io::{read_de_results, write_de_ensemble, EnsembleRow};
    use tempfile::TempDir;

    let requested: Vec<String> = args
        .ensemble_methods
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    if requested.is_empty() {
        anyhow::bail!("--ensemble-methods resolved to zero methods");
    }
    for m in &requested {
        if !matches!(
            m.as_str(),
            "paired-t" | "moderated" | "welch-t" | "ols" | "mixed" | "limma" | "msqrob"
        ) {
            anyhow::bail!("unsupported ensemble method {:?}", m);
        }
    }

    let thresholds = GradeThresholds {
        q_threshold: args.ensemble_q_threshold,
        validated_sign_fraction: args.ensemble_sign_fraction,
        provisional_sign_fraction: args.ensemble_provisional_fraction,
    };
    // `--ensemble-validated-fraction` is kept in the CLI surface for
    // future per-method significance gates; the current grading uses
    // ensemble_q + sign consistency only.
    let _ = args.ensemble_validated_fraction;

    let tmp = TempDir::new().with_context(|| "creating ensemble tempdir")?;
    let mut merged_rows: Vec<DeResultRow> = Vec::new();
    let mut skipped: Vec<String> = Vec::new();
    let mut applied: Vec<String> = Vec::new();

    for method in &requested {
        let mut sub = args.clone();
        sub.test = method.clone();
        sub.output_dir = tmp.path().join(method);
        // Strip peptide flags for methods that don't consume them.
        match method.as_str() {
            "msqrob" => { /* keeps both */ }
            "limma" => {
                sub.peptide_measurements = None;
            }
            _ => {
                sub.peptide_measurements = None;
                sub.peptide_metadata = None;
            }
        }
        // Strip formula / mixed flags for methods that can't use them.
        match method.as_str() {
            "ols" => {
                sub.fixed = None;
                sub.random = None;
            }
            "mixed" => {
                sub.covariates = None;
                sub.design = None;
                sub.contrast = None;
                sub.per_subject_proxy = None;
            }
            _ => {
                sub.covariates = None;
                sub.design = None;
                sub.contrast = None;
                sub.per_subject_proxy = None;
                sub.fixed = None;
                sub.random = None;
            }
        }
        match run(sub.clone()) {
            Ok(()) => {
                let results_path = sub.output_dir.join("de_results.tsv");
                let rows = read_de_results(&results_path).with_context(|| {
                    format!("reading sub-method de_results for {method:?}")
                })?;
                merged_rows.extend(rows);
                applied.push(method.clone());
            }
            Err(e) => {
                let short = e.to_string().lines().next().unwrap_or("").to_string();
                eprintln!("ensemble: method {method:?} skipped — {short}");
                skipped.push(format!("{method}:{short}"));
            }
        }
    }

    // Aggregate per (comparison, assay_id).
    use std::collections::BTreeMap;
    type Key = (String, String);
    let mut grouped: BTreeMap<
        Key,
        (Vec<EnsembleInput>, String, String, String), // (inputs, panel, gene, uniprot)
    > = BTreeMap::new();
    for r in &merged_rows {
        let key = (r.comparison.clone(), r.assay_id.clone());
        let entry = grouped.entry(key).or_insert_with(|| {
            (
                Vec::new(),
                r.panel.clone(),
                r.gene_symbol.clone(),
                r.uniprot.clone(),
            )
        });
        entry.0.push(EnsembleInput {
            method: r.method.clone(),
            mean_diff: r.mean_diff,
            p_value: r.p_value,
            bh_q: r.bh_q,
        });
    }

    let mut ensemble_rows: Vec<EnsembleRow> = Vec::new();
    let mut by_comparison: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    let methods_skipped_str = if skipped.is_empty() {
        String::new()
    } else {
        skipped.join("|")
    };
    for ((comparison, assay_id), (inputs, panel, gene, uniprot)) in &grouped {
        let agg = aggregate_per_protein(inputs, thresholds);
        let idx = ensemble_rows.len();
        by_comparison
            .entry(comparison.clone())
            .or_default()
            .push(idx);
        ensemble_rows.push(EnsembleRow {
            comparison: comparison.clone(),
            panel: panel.clone(),
            assay_id: assay_id.clone(),
            gene_symbol: gene.clone(),
            uniprot: uniprot.clone(),
            n_applied: agg.n_applied,
            n_significant: agg.n_significant,
            n_sign_consistent: agg.n_sign_consistent,
            majority_sign: agg.majority_sign,
            ensemble_p: agg.ensemble_p,
            ensemble_q: None,
            grade: String::new(), // filled after BH
            methods_applied: agg.methods_applied.join(","),
            methods_skipped: methods_skipped_str.clone(),
        });
    }
    // BH on ensemble_p within each comparison, then assign grade per
    // protein from (ensemble_q, sign consistency).
    for indices in by_comparison.values() {
        let ps: Vec<Option<f64>> = indices.iter().map(|&i| ensemble_rows[i].ensemble_p).collect();
        let qs = atman_core::bh_fdr(&ps);
        for (j, &i) in indices.iter().enumerate() {
            ensemble_rows[i].ensemble_q = qs[j];
            let row = &ensemble_rows[i];
            let grade = assign_grade(
                row.ensemble_q,
                row.n_applied,
                row.n_sign_consistent,
                thresholds,
            );
            ensemble_rows[i].grade = grade.as_str().to_string();
        }
    }

    // Stable output ordering: (comparison, grade-priority, ensemble_q asc).
    ensemble_rows.sort_by(|a, b| {
        a.comparison
            .cmp(&b.comparison)
            .then_with(|| grade_rank(&a.grade).cmp(&grade_rank(&b.grade)))
            .then_with(|| match (a.ensemble_q, b.ensemble_q) {
                (Some(qa), Some(qb)) => qa.partial_cmp(&qb).unwrap_or(std::cmp::Ordering::Equal),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            })
            .then_with(|| a.assay_id.cmp(&b.assay_id))
    });
    merged_rows.sort_by(|a, b| {
        a.comparison
            .cmp(&b.comparison)
            .then_with(|| a.assay_id.cmp(&b.assay_id))
            .then_with(|| a.method.cmp(&b.method))
    });

    let results_path = args.output_dir.join("de_results.tsv");
    let ensemble_path = args.output_dir.join("de_ensemble.tsv");
    write_de_results(&results_path, &merged_rows)?;
    write_de_ensemble(&ensemble_path, &ensemble_rows)?;

    let finished_at = SystemTime::now();
    let input_dir_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &[
            "qc_measurements.tsv",
            "measurements.tsv",
            "samples.tsv",
            "proteins.tsv",
        ],
    )?;
    let outputs: Vec<PathBuf> = vec![results_path.clone(), ensemble_path.clone()];
    let sidecar = sidecar_path_for(&results_path);
    write_run_sidecar(
        &sidecar,
        "de",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output-dir": args.output_dir.display().to_string(),
            "test": "ensemble",
            "groups": args.groups,
            "min-pairs": args.min_pairs,
            "peptide-measurements": args.peptide_measurements.as_ref().map(|p| p.display().to_string()),
            "peptide-metadata": args.peptide_metadata.as_ref().map(|p| p.display().to_string()),
            "ensemble-methods": args.ensemble_methods,
            "ensemble-applied-methods": applied.join(","),
            "ensemble-skipped-methods": skipped.join("|"),
            "ensemble-q-threshold": args.ensemble_q_threshold,
            "ensemble-validated-fraction": args.ensemble_validated_fraction,
            "ensemble-provisional-fraction": args.ensemble_provisional_fraction,
            "ensemble-sign-fraction": args.ensemble_sign_fraction,
        }),
        &input_dir_sha256,
        &outputs,
        started_at,
        finished_at,
    )?;
    eprintln!(
        "de ensemble: {} rows, {} applied: [{}], {} skipped: [{}]",
        ensemble_rows.len(),
        applied.len(),
        applied.join(", "),
        skipped.len(),
        skipped.join(", ")
    );
    eprintln!("de: sidecar={}", sidecar.display());
    Ok(())
}

fn grade_rank(grade: &str) -> u8 {
    match grade {
        "VALIDATED" => 0,
        "PROVISIONAL" => 1,
        _ => 2,
    }
}

/// Parse a `--contrast-list` like `"MCI-CN,AD-CN,AD-MCI"` into
/// owned `(level_a, level_b)` pairs. Rejects empty tokens and
/// malformed entries.
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
fn run_posthoc_sidak(args: Args, started_at: SystemTime) -> Result<()> {
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
    let (design_rows, design_labels, cov_kinds) =
        build_posthoc_full_design(&samples, &setup, &factor_name)?;
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
    let input_dir_sha256 = hash_canonical_inputs(
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
        &input_dir_sha256,
        &outputs,
        started_at,
        finished_at,
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
fn run_posthoc_tukey(args: Args, started_at: SystemTime) -> Result<()> {
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

    let (design_rows, design_labels, cov_kinds) =
        build_posthoc_full_design(&samples, &setup, &factor_name)?;
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
    let input_dir_sha256 = hash_canonical_inputs(
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
        &input_dir_sha256,
        &outputs,
        started_at,
        finished_at,
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
fn run_posthoc_dunnett(args: Args, started_at: SystemTime) -> Result<()> {
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

    let (design_rows, design_labels, cov_kinds) =
        build_posthoc_full_design(&samples, &setup, &factor_name)?;
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
    let input_dir_sha256 = hash_canonical_inputs(
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
        &input_dir_sha256,
        &outputs,
        started_at,
        finished_at,
    )?;
    eprintln!("de posthoc dunnett: sidecar={}", sidecar.display());
    Ok(())
}

fn build_posthoc_full_design(
    samples: &[Sample],
    setup: &OlsSetup,
    _factor_name: &str,
) -> Result<(Vec<(String, Vec<f64>)>, Vec<String>, Vec<CovKind>)> {
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
    Ok((rows, design_labels, cov_kinds))
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

fn write_design_rows(path: &Path, rows: &[DesignReportRow]) -> Result<()> {
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
