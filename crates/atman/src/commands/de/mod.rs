//! Differential abundance command.

mod ensemble;
mod limma;
mod moderated;
mod msqrob;
mod posthoc;
mod preflight;
mod robust_stats;

use ensemble::run_ensemble;
use limma::run_limma;
use moderated::apply_moderated_shrinkage;
use msqrob::run_msqrob;
use posthoc::{run_posthoc_dunnett, run_posthoc_sidak, run_posthoc_tukey, write_design_rows};
use robust_stats::{median, robust_paired, robust_unpaired, RobustStats};

use anyhow::{anyhow, Context, Result};
use atman_core::de::{
    bh_fdr, mixed_random_intercept, ols, paired_t, welch_t, OlsOutcome, PairedTResult, SkipReason,
};
use atman_core::Sample;
use clap::Args as ClapArgs;
use serde_json::{json, Map as JsonMap};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use super::parse_comparisons;
use crate::io::{
    atomic_write, hash_canonical_inputs, hash_labeled_inputs, read_measurements_long,
    read_proteins, read_samples, sidecar_path_for, write_de_report, write_de_results,
    write_run_sidecar, DeReportRow, DeResultRow,
};

#[derive(ClapArgs, Debug, Clone)]
pub struct Args {
    /// Directory containing measurements.tsv, samples.tsv, proteins.tsv.
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

    /// Biological-replicate key. Currently fixed to the canonical
    /// `subject_id` column from samples.tsv; the legacy alias
    /// `participant` is also accepted for backward compatibility.
    /// Required for `paired-t` and `moderated`; ignored for the
    /// unpaired tests (`welch-t`, `ols`, `mixed`, `limma`, `msqrob`,
    /// `ensemble`).
    #[arg(long, default_value = "subject_id")]
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

    /// Lower Winsor tail fraction for the limma robust prior fit. Must be
    /// in `(0, 0.5)`. Ignored when `--robust false`. Default matches limma's
    /// `winsor.tail.p[1]`.
    #[arg(long, default_value_t = 0.05)]
    limma_winsor_lower: f64,

    /// Upper Winsor tail fraction for the limma robust prior fit. Must be
    /// in `(0, 0.5)` and `lower + upper < 1`. Ignored when `--robust false`.
    /// Default matches limma's `winsor.tail.p[2]`.
    #[arg(long, default_value_t = 0.10)]
    limma_winsor_upper: f64,

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

    /// Sign-fraction threshold for the PROVISIONAL grade (below
    /// VALIDATED, above INSUFFICIENT). At least this fraction of
    /// methods must agree with the majority sign on the protein's
    /// effect.
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

    /// Maximum allowed fraction of non-dropped measurements flagged
    /// `below_lod=1` in `measurements.tsv`. Default 0.5. Left-censored
    /// (MNAR) data above this fraction is likely to distort unpaired-t,
    /// limma, and OLS tests that treat missing-at-random by drop. Raise
    /// the threshold, pass `--allow-censored`, or impute/filter upstream.
    #[arg(long, default_value_t = 0.5)]
    max_below_lod_fraction: f64,

    /// Bypass the below-LOD safety gate. Set only when you have confirmed
    /// the chosen test handles left-censored data (e.g. an explicit
    /// MNAR-aware model upstream, or after imputation).
    #[arg(long, default_value_t = false)]
    allow_censored: bool,

    /// Strict BH-q threshold counted in the `n_q_strict` column of
    /// `de_report.tsv`. Reporting policy only — does not change which
    /// tests are run, what is written to `de_results.tsv`, or the
    /// post-hoc `decision` column (see `--alpha`). Must satisfy
    /// `0 < strict < relaxed < 1`. Resolved value is stamped in the
    /// run sidecar's `report_thresholds` block.
    #[arg(long, default_value_t = 0.05)]
    report_q_strict: f64,

    /// Relaxed BH-q threshold counted in the `n_q_relaxed` column of
    /// `de_report.tsv`. Same scope as `--report-q-strict`. Must
    /// satisfy `0 < strict < relaxed < 1`.
    #[arg(long, default_value_t = 0.10)]
    report_q_relaxed: f64,

    /// p-value threshold counted in the `n_p_strict` column of the
    /// per-subject-proxy summary TSV. Reporting policy only. Must lie
    /// in `(0, 1)`.
    #[arg(long, default_value_t = 0.05)]
    proxy_p_threshold: f64,

    /// External covariate TSV(s) for `--test limma`. Wide format:
    /// `sample_id` column followed by one column per numeric covariate.
    /// Multiple paths append all covariates (joined by sample_id; error
    /// on missing sample IDs). Covariate columns from each file must not
    /// overlap each other or the samples.tsv columns.
    /// Distinct from `--covariates` (which reads from samples.tsv).
    /// Sidecar records each path and its SHA-256.
    #[arg(long, action = clap::ArgAction::Append)]
    adjust_for: Vec<PathBuf>,
}

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
    if !is_unpaired && args.paired_by != "subject_id" {
        anyhow::bail!(
            "paired-by {:?} not supported for paired tests; pairing currently uses the \
             canonical `subject_id` column from samples.tsv.",
            args.paired_by
        );
    }
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
        other => anyhow::bail!(
            "unknown --post-hoc {other:?}; supported: sidak (tukey / dunnett pending)"
        ),
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
        // cells_by_sample is used by OLS, mixed, and any path routing to OLS
        // (welch-t + adjust-for). paired-t + adjust-for also needs per-sample
        // lookup to compute paired differences.
        if args.test == "ols"
            || args.test == "mixed"
            || (!args.adjust_for.is_empty()
                && (args.test == "welch-t" || args.test == "paired-t"))
        {
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
    let external_covariates: Option<std::collections::BTreeMap<String, std::collections::BTreeMap<String, f64>>> =
        if args.adjust_for.is_empty() {
            None
        } else {
            use atman_core::de::{join_external_covariates, read_external_covariates};
            let mut merged: std::collections::BTreeMap<String, std::collections::BTreeMap<String, f64>> =
                std::collections::BTreeMap::new();
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
        let (limma_rows, limma_reports) =
            run_limma(&args, &samples, &proteins, &measurements, &comparisons, external_covariates.as_ref())?;
        all_rows.extend(limma_rows);
        report_rows.extend(limma_reports);
    }

    if args.test == "msqrob" {
        let (msqrob_rows, msqrob_reports) =
            run_msqrob(&args, &samples, &proteins, &comparisons, external_covariates.as_ref())?;
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
            // welch-t + adjust-for also routes here (plain OLS SE).
            let model_design: Option<OlsDesign> = if is_ols_routed {
                let setup = model_setup
                    .as_ref()
                    .expect("model_setup populated when is_ols_routed");
                let mut design =
                    build_ols_design(&samples, comp_a, comp_b, setup)?;
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
                } else if args.test == "ols"
                    || (is_ols_routed && args.test == "welch-t")
                {
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
                } else if args.test == "paired-t" && !args.adjust_for.is_empty() {
                    // paired-t + --adjust-for: paired-difference ANCOVA.
                    // Computes d_k = abundance_b_k - abundance_a_k per subject,
                    // covariate_diff_k = cov_b_k - cov_a_k, then fits
                    // OLS: d ~ 1 + cov_diff_1 + cov_diff_2 + ...
                    // The intercept is the adjusted mean paired difference.
                    let ext = external_covariates.as_ref().expect("external_covariates set");
                    let result = paired_t_covariate_adjusted(
                        va, vb, ext, &sample_by_id, comp_a, comp_b, args.min_pairs,
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
        &[
            "measurements.tsv",
            "measurements.tsv",
            "samples.tsv",
            "proteins.tsv",
        ],
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
    n_q_strict: usize,
    n_q_relaxed: usize,
    min_q: Option<f64>,
    max_abs_effect: Option<f64>,
}

// -----------------------------------------------------------------
// OLS (covariate-adjusted) helpers
// -----------------------------------------------------------------

#[derive(Debug, Clone)]
pub(super) enum DesignTerm {
    Condition,
    Covariate(String),
}

pub(super) struct OlsSetup {
    terms: Vec<DesignTerm>,
    cov_names: Vec<String>,
    cov_raw: HashMap<String, Vec<Option<String>>>,
    contrast: Option<String>,
    label: String,
}

pub(super) fn build_ols_setup(
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

pub(super) fn parse_design_terms(formula: &str) -> Result<Vec<DesignTerm>> {
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

pub(super) fn covariate_names_from_terms(terms: &[DesignTerm]) -> Vec<String> {
    terms
        .iter()
        .filter_map(|t| match t {
            DesignTerm::Condition => None,
            DesignTerm::Covariate(name) => Some(name.clone()),
        })
        .collect()
}

pub(super) fn resolve_comparisons(
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
pub(super) enum CovKind {
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
pub(super) struct OmnibusRow {
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

pub(super) fn write_omnibus_rows(path: &Path, rows: &[OmnibusRow]) -> Result<()> {
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

pub(super) fn format_opt(v: f64) -> String {
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
pub(super) struct CovariateRow {
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
pub(super) struct DesignReportRow {
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
pub(super) fn read_covariate_columns(
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
pub(super) fn classify_covariates(
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

pub(super) fn push_encoded_covariate(
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
pub(super) fn write_covariate_rows(path: &Path, rows: &[CovariateRow]) -> Result<()> {
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

fn write_proxy_summary(
    path: &Path,
    proxy: &str,
    rows: &[CovariateRow],
    p_threshold: f64,
) -> Result<()> {
    let mut by_comparison: BTreeMap<&str, Vec<&CovariateRow>> = BTreeMap::new();
    for row in rows {
        if row.covariate == proxy {
            by_comparison
                .entry(row.comparison.as_str())
                .or_default()
                .push(row);
        }
    }
    let mut buf =
        String::from("proxy\tcomparison\tn_tests\tn_p_strict\tmedian_abs_beta\tmedian_p_value\n");
    for (comparison, rows) in by_comparison {
        let n_tests = rows.len();
        let n_p_strict = rows.iter().filter(|row| row.p_value < p_threshold).count();
        let median_abs_beta = median(rows.iter().map(|row| row.beta.abs()).collect());
        let median_p_value = median(rows.iter().map(|row| row.p_value).collect());
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            proxy,
            comparison,
            n_tests,
            n_p_strict,
            fmt_opt(median_abs_beta),
            fmt_opt(median_p_value),
        ));
    }
    atomic_write(path, buf.as_bytes())
}

pub(super) fn fmt_opt(value: Option<f64>) -> String {
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
pub(super) fn build_two_group_design(
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

/// Augment an already-built `OlsDesign` with external covariate columns from
/// `--adjust-for`. Each covariate is appended as a numeric column to every
/// design row; samples whose `sample_id` is absent from `ext` are dropped
/// (with a drop-reason record). The `group_col` index is unchanged because
/// external columns are appended after the existing columns.
///
/// External covariates are always treated as numeric (the caller has already
/// parsed them as `f64` via `read_external_covariates`).
fn augment_ols_design_with_external(
    design: &mut OlsDesign,
    ext: &std::collections::BTreeMap<String, std::collections::BTreeMap<String, f64>>,
    comp_a: &str,
    comp_b: &str,
) -> Result<()> {
    if ext.is_empty() {
        return Ok(());
    }

    // Collect sorted covariate names from the first sample that's in the design.
    // All samples must have the same covariate names (enforced by join_external_covariates).
    let ext_cov_names: Vec<String> = {
        let sample_with_covs = design
            .rows
            .iter()
            .find_map(|(sid, _)| ext.get(sid.as_str()));
        match sample_with_covs {
            Some(covs) => covs.keys().cloned().collect(),
            None => return Ok(()), // no samples in ext — nothing to augment
        }
    };

    if ext_cov_names.is_empty() {
        return Ok(());
    }

    // Drop rows whose sample_id is missing from ext, keep the rest.
    let comparison = format!("{comp_a}-{comp_b}");
    let mut new_rows: Vec<(String, Vec<f64>)> = Vec::new();
    for (sid, mut row) in design.rows.drain(..) {
        match ext.get(&sid) {
            Some(cov_map) => {
                for name in &ext_cov_names {
                    let val = cov_map.get(name).copied().unwrap_or(f64::NAN);
                    row.push(val);
                }
                new_rows.push((sid, row));
            }
            None => {
                // Sample missing from external covariate file: drop with report.
                design.report_rows.push(DesignReportRow {
                    comparison: comparison.clone(),
                    sample_id: sid,
                    condition: String::new(),
                    included: false,
                    drop_reason: "missing_external_covariate_row".to_string(),
                    columns: String::new(),
                    values: String::new(),
                });
            }
        }
    }
    design.rows = new_rows;

    // Append column labels for the external covariates.
    for name in &ext_cov_names {
        design.design_labels.push(name.clone());
    }

    Ok(())
}

/// Paired-t with covariates via paired-difference ANCOVA.
///
/// For each matched subject k: `d_k = abundance_b_k - abundance_a_k`.
/// For each external covariate c: `dc_k = cov_c(sample_b_k) - cov_c(sample_a_k)`.
/// Fits OLS: `d ~ 1 + dc_1 + dc_2 + ...` per protein.
/// The intercept is the adjusted mean paired difference (log_fc).
/// t and p are from the intercept coefficient with df = n_pairs - (1 + n_cov).
///
/// `va` and `vb` are `(subject_id, abundance)` pairs for conditions a and b.
/// `ext` maps sample_id → covariate_name → value.
/// `comp_a` and `comp_b` are the condition labels (used to look up sample_ids).
#[allow(clippy::too_many_arguments)]
fn paired_t_covariate_adjusted(
    va: &[(String, f64)],
    vb: &[(String, f64)],
    ext: &std::collections::BTreeMap<String, std::collections::BTreeMap<String, f64>>,
    sample_by_id: &HashMap<&str, &atman_core::Sample>,
    comp_a: &str,
    comp_b: &str,
    min_pairs: usize,
) -> PairedTResult {
    use atman_core::de::{OlsOutcome, SkipReason};

    // Build subject → abundance maps.
    let a_by_subj: HashMap<&str, f64> = va.iter().map(|(s, v)| (s.as_str(), *v)).collect();
    let b_by_subj: HashMap<&str, f64> = vb.iter().map(|(s, v)| (s.as_str(), *v)).collect();

    // Build subject → sample_id lookup for each condition from sample_by_id.
    let mut subj_to_sid_a: HashMap<&str, &str> = HashMap::new();
    let mut subj_to_sid_b: HashMap<&str, &str> = HashMap::new();
    for (sid, s) in sample_by_id.iter() {
        if s.is_control {
            continue;
        }
        let subj = match s.subject_id.as_deref() {
            Some(id) => id,
            None => continue,
        };
        match s.condition.as_deref() {
            Some(c) if c == comp_a => { subj_to_sid_a.insert(subj, sid); }
            Some(c) if c == comp_b => { subj_to_sid_b.insert(subj, sid); }
            _ => {}
        }
    }

    // Collect sorted covariate names.
    let ext_cov_names: Vec<String> = {
        let any = ext.values().next();
        match any {
            Some(m) => m.keys().cloned().collect(),
            None => return PairedTResult::Skipped {
                reason: SkipReason::InsufficientPairs,
                n_pairs: 0,
            },
        }
    };
    let n_cov = ext_cov_names.len();

    // Build per-subject: (d_k, [dc_k_1, dc_k_2, ...]).
    let mut diffs: Vec<f64> = Vec::new();
    let mut cov_diffs: Vec<Vec<f64>> = Vec::new();

    let subjects: std::collections::BTreeSet<&str> = a_by_subj
        .keys()
        .copied()
        .filter(|s| b_by_subj.contains_key(s))
        .collect();

    for subj in &subjects {
        let abund_a = match a_by_subj.get(subj) { Some(v) => *v, None => continue };
        let abund_b = match b_by_subj.get(subj) { Some(v) => *v, None => continue };
        let sid_a = match subj_to_sid_a.get(subj) { Some(s) => *s, None => continue };
        let sid_b = match subj_to_sid_b.get(subj) { Some(s) => *s, None => continue };

        let cov_a = match ext.get(sid_a) { Some(m) => m, None => continue };
        let cov_b = match ext.get(sid_b) { Some(m) => m, None => continue };

        let mut cds: Vec<f64> = Vec::with_capacity(n_cov);
        let mut valid = true;
        for name in &ext_cov_names {
            let va_cov = match cov_a.get(name) { Some(v) => *v, None => { valid = false; break; } };
            let vb_cov = match cov_b.get(name) { Some(v) => *v, None => { valid = false; break; } };
            // cov_a is the comp_a sample (e.g. condition "b"), cov_b is comp_b (e.g. condition "a").
            // Covariate diff = comp_a_cov - comp_b_cov, consistent with abundance diff sign.
            let d = va_cov - vb_cov;
            if !d.is_finite() { valid = false; break; }
            cds.push(d);
        }
        if !valid { continue; }
        if !abund_a.is_finite() || !abund_b.is_finite() { continue; }

        // diff = comp_a_value - comp_b_value, consistent with regular paired_t
        // which computes diffs as a_i - b_i where a = va[subj] (comp_a condition).
        // With --groups b-a: comp_a = "b", so diff = condition_b - condition_a
        // (positive when b > a, matching atman sign convention).
        diffs.push(abund_a - abund_b);
        cov_diffs.push(cds);
    }

    let n_pairs = diffs.len();
    if n_pairs < min_pairs || n_pairs < 2 {
        return PairedTResult::Skipped {
            reason: SkipReason::InsufficientPairs,
            n_pairs,
        };
    }

    // Build design matrix for OLS: [1, dc_1, dc_2, ...]
    // The intercept (column 0) is the adjusted mean paired difference.
    let p = 1 + n_cov;
    let design: Vec<Vec<f64>> = diffs
        .iter()
        .zip(cov_diffs.iter())
        .map(|(_, cds)| {
            let mut row = Vec::with_capacity(p);
            row.push(1.0);
            row.extend_from_slice(cds);
            row
        })
        .collect();

    // Use atman_core::de::ols() to fit d ~ 1 + cov_diffs.
    // The intercept (column 0) is the adjusted mean diff.
    let fit = match atman_core::de::ols(&design, &diffs, min_pairs) {
        OlsOutcome::Computed(f) => f,
        OlsOutcome::Skipped { reason, n } => {
            return PairedTResult::Skipped { reason, n_pairs: n };
        }
    };

    // mean_a / mean_b are not readily available from the paired-diff design
    // (we only have diffs, not the paired subset's individual abundances).
    // Report NaN; the key outputs (mean_diff, t, p) are fully determined.

    PairedTResult::Computed {
        n_pairs,
        mean_a: f64::NAN,
        mean_b: f64::NAN,
        mean_diff: fit.beta[0], // intercept = adjusted mean diff
        t: fit.t[0],
        df: fit.df,
        p_value: fit.p_value[0],
    }
}
