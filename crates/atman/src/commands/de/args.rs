//! Clap `Args` for the differential-abundance command plus the
//! `--per-subject-proxy` preprocessing that rewrites the parsed args
//! before the `run` dispatch executes.

use anyhow::Result;
use clap::Args as ClapArgs;
use std::path::PathBuf;

use super::{parse_design_terms, DesignTerm};

#[derive(ClapArgs, Debug, Clone)]
pub struct Args {
    /// Directory containing measurements.tsv, samples.tsv, proteins.tsv.
    #[arg(long)]
    pub(super) input_dir: PathBuf,

    /// Output directory for de_results.tsv + de_report.tsv.
    #[arg(long)]
    pub(super) output_dir: PathBuf,

    /// Test kind: `paired-t` (default), `moderated` (paired variance
    /// shrinkage), `welch-t` (unpaired two-sample with Welch-Satterthwaite
    /// df), or `ols` (covariate-adjusted unpaired linear model).
    #[arg(long, default_value = "paired-t")]
    pub(super) test: String,

    /// Biological-replicate key: the samples.tsv column used to pair
    /// observations. Defaults to the canonical `subject_id` column, but
    /// may be overridden to any column present in samples.tsv; its value
    /// then replaces the in-memory `subject_id` used for the paired join.
    /// The legacy alias `participant` is also accepted for backward
    /// compatibility. Required for `paired-t` and `moderated`; ignored
    /// for the unpaired tests (`welch-t`, `ols`, `mixed`, `limma`,
    /// `msqrob`, `ensemble`).
    #[arg(long, default_value = "subject_id")]
    pub(super) paired_by: String,

    /// Comma-separated covariate column names from samples.tsv for
    /// `--test ols`. Numeric columns (parseable as f64 in every non-missing
    /// row) are used as continuous predictors; everything else is treated
    /// as categorical and one-hot encoded, dropping the alphabetically
    /// first level as reference. Samples with any missing covariate value
    /// are dropped before fitting. Ignored for other tests.
    #[arg(long)]
    pub(super) covariates: Option<String>,

    /// Formula-style design for `--test ols`, for example
    /// `~ condition + age + sex + batch`. The intercept is included
    /// automatically. `condition` is the Atman sample condition; other terms
    /// are read from samples.tsv.
    #[arg(long)]
    pub(super) design: Option<String>,

    /// Coefficient to test for `--test ols --design`, for example
    /// `conditionCase`. If omitted with --groups, Atman uses the comparison
    /// coefficient for each `A-B` group.
    #[arg(long)]
    pub(super) contrast: Option<String>,

    /// Physiological subject-level proxy from samples.tsv to include as a continuous OLS regressor.
    #[arg(long)]
    pub(super) per_subject_proxy: Option<String>,

    /// Fixed-effect terms for `--test mixed`, for example
    /// `condition + age + sex`. The intercept is included automatically.
    #[arg(long)]
    pub(super) fixed: Option<String>,

    /// Random-effect structure for `--test mixed`. Initial support is
    /// random intercept by subject only: `1|subject_id`.
    #[arg(long)]
    pub(super) random: Option<String>,

    /// Comma-separated comparisons in `A-B` form. Each is a separate
    /// hypothesis family for FDR. Example: "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2".
    #[arg(long)]
    pub(super) groups: Option<String>,

    /// Minimum number of samples required per group for a test to run.
    /// For paired tests this is the number of matched pairs; for `welch-t`
    /// it is the minimum of `|group_a|` and `|group_b|`. Below this the
    /// protein is emitted as Skipped(InsufficientPairs) with NaN p/q.
    #[arg(long, default_value_t = 5)]
    pub(super) min_pairs: usize,

    /// Prior degrees of freedom for `--test moderated` variance shrinkage.
    #[arg(long, default_value_t = 4.0)]
    pub(super) moderation_prior_df: f64,

    /// Minimum log-fold-change threshold for TREAT testing under `--test limma`.
    /// `0.0` (default) reduces to the standard moderated-t p-value.
    #[arg(long, default_value_t = 0.0)]
    pub(super) lfc_threshold: f64,

    /// Fit the mean-variance trend before eBayes shrinkage (limma only).
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set, num_args = 1)]
    pub(super) trend: bool,

    /// Use the robust (Winsorized) prior fit (limma only).
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set, num_args = 1)]
    pub(super) robust: bool,

    /// Lower Winsor tail fraction for the limma robust prior fit. Must be
    /// in `(0, 0.5)`. Ignored when `--robust false`. Default matches limma's
    /// `winsor.tail.p[1]`.
    #[arg(long, default_value_t = 0.05)]
    pub(super) limma_winsor_lower: f64,

    /// Upper Winsor tail fraction for the limma robust prior fit. Must be
    /// in `(0, 0.5)` and `lower + upper < 1`. Ignored when `--robust false`.
    /// Default matches limma's `winsor.tail.p[2]`.
    #[arg(long, default_value_t = 0.10)]
    pub(super) limma_winsor_upper: f64,

    /// Peptide-level long TSV (one row per sample × peptide).
    /// Required for `--test msqrob`. Rejected for other tests.
    /// Schema: `sample_id, peptide_id, abundance, abundance_unit,
    /// dropped_by_qc, below_lod`.
    #[arg(long)]
    pub(super) peptide_measurements: Option<PathBuf>,

    /// Peptide catalog TSV mapping `peptide_id → assay_id` (parent
    /// protein) plus optional `sequence, charge, modifications,
    /// missed_cleavages`. Required for `--test msqrob`.
    #[arg(long)]
    pub(super) peptide_metadata: Option<PathBuf>,

    /// L2 ridge penalty on non-intercept fixed-effect coefficients
    /// for `--test msqrob`. Numeric value ≥ 0 or `auto` (currently
    /// equivalent to `0.0`; data-driven selection is a follow-on).
    #[arg(long, default_value = "auto")]
    pub(super) ridge_lambda: String,

    /// Minimum distinct peptides observed per protein required for
    /// `--test msqrob`. Below this, the protein emits
    /// `Skipped(insufficient_peptides)`.
    #[arg(long, default_value_t = 2)]
    pub(super) min_peptides: usize,

    /// Comma-separated list of methods to run under `--test ensemble`.
    /// Each of `paired-t`, `welch-t`, `ols`, `mixed`, `limma`, `msqrob`
    /// is valid. Methods whose required inputs are missing
    /// (`--peptide-measurements` for msqrob, a paired-subject layout
    /// for paired-t, etc.) are auto-skipped with a note in
    /// `methods_skipped` rather than aborting the run.
    #[arg(long, default_value = "welch-t,ols,limma,msqrob")]
    pub(super) ensemble_methods: String,

    /// Per-method BH-q significance threshold for ensemble grading.
    #[arg(long, default_value_t = 0.05)]
    pub(super) ensemble_q_threshold: f64,

    /// Sign-fraction threshold for the PROVISIONAL grade (below
    /// VALIDATED, above INSUFFICIENT). At least this fraction of
    /// methods must agree with the majority sign on the protein's
    /// effect.
    #[arg(long, default_value_t = 0.50)]
    pub(super) ensemble_provisional_fraction: f64,

    /// Fraction of methods whose mean_diff sign must match the
    /// majority sign for either VALIDATED or PROVISIONAL. Default
    /// 1.00 — any directional disagreement drops the grade.
    #[arg(long, default_value_t = 1.00)]
    pub(super) ensemble_sign_fraction: f64,

    /// Name of the categorical factor (a sample metadata column
    /// referenced by `--design`) whose omnibus F-test is emitted to
    /// `de_omnibus.tsv`. Requires `--test ols` and a `--design` that
    /// includes the factor. At least 3 levels required.
    #[arg(long)]
    pub(super) omnibus_factor: Option<String>,

    /// Post-hoc contrast adjustment method. `""` (default) disables
    /// post-hoc; `sidak` uses `p_adj = 1 − (1 − p)^m` across the
    /// `--contrast-list`; `tukey` runs the Tukey HSD studentized-range
    /// adjustment over all pairs of the post-hoc factor's levels; and
    /// `dunnett` compares each non-control level to the control level.
    /// All three are fully implemented and dispatched. Each requires
    /// `--test ols` with a `--design`.
    #[arg(long, default_value = "")]
    pub(super) post_hoc: String,

    /// Comma-separated `level_a-level_b` contrast list applied to
    /// `--post-hoc-factor` (or `--omnibus-factor` as fallback).
    /// Example: `"MCI-CN,AD-CN,AD-MCI"` for `stage` ∈ {CN, MCI, AD}.
    #[arg(long)]
    pub(super) contrast_list: Option<String>,

    /// Factor whose levels `--contrast-list` references. Falls back
    /// to `--omnibus-factor` when unset.
    #[arg(long)]
    pub(super) post_hoc_factor: Option<String>,

    /// Family-wise significance threshold for the post-hoc
    /// `decision` column in `de_results.tsv`.
    #[arg(long, default_value_t = 0.05)]
    pub(super) alpha: f64,

    /// Maximum allowed fraction of non-dropped measurements flagged
    /// `below_lod=1` in `measurements.tsv`. Default 0.5. Left-censored
    /// (MNAR) data above this fraction is likely to distort unpaired-t,
    /// limma, and OLS tests that treat missing-at-random by drop. Raise
    /// the threshold, pass `--allow-censored`, or impute/filter upstream.
    #[arg(long, default_value_t = 0.5)]
    pub(super) max_below_lod_fraction: f64,

    /// Bypass the below-LOD safety gate. Set only when you have confirmed
    /// the chosen test handles left-censored data (e.g. an explicit
    /// MNAR-aware model upstream, or after imputation).
    #[arg(long, default_value_t = false)]
    pub(super) allow_censored: bool,

    /// Strict BH-q threshold counted in the `n_q_strict` column of
    /// `de_report.tsv`. Reporting policy only — does not change which
    /// tests are run, what is written to `de_results.tsv`, or the
    /// post-hoc `decision` column (see `--alpha`). Must satisfy
    /// `0 < strict < relaxed < 1`. Resolved value is stamped in the
    /// run sidecar's `report_thresholds` block.
    #[arg(long, default_value_t = 0.05)]
    pub(super) report_q_strict: f64,

    /// Relaxed BH-q threshold counted in the `n_q_relaxed` column of
    /// `de_report.tsv`. Same scope as `--report-q-strict`. Must
    /// satisfy `0 < strict < relaxed < 1`.
    #[arg(long, default_value_t = 0.10)]
    pub(super) report_q_relaxed: f64,

    /// p-value threshold counted in the `n_p_strict` column of the
    /// per-subject-proxy summary TSV. Reporting policy only. Must lie
    /// in `(0, 1)`.
    #[arg(long, default_value_t = 0.05)]
    pub(super) proxy_p_threshold: f64,

    /// External covariate TSV(s) for `--test limma`. Wide format:
    /// `sample_id` column followed by one column per numeric covariate.
    /// Multiple paths append all covariates (joined by sample_id; error
    /// on missing sample IDs). Covariate columns from each file must not
    /// overlap each other or the samples.tsv columns.
    /// Distinct from `--covariates` (which reads from samples.tsv).
    /// Sidecar records each path and its SHA-256.
    #[arg(long, action = clap::ArgAction::Append)]
    pub(super) adjust_for: Vec<PathBuf>,

    /// Let samples flagged `is_control=1` join a group when their condition
    /// is named in `--groups`. Default: control samples are excluded from
    /// every comparison (historical behaviour).
    #[arg(long, default_value_t = false)]
    pub(super) include_controls: bool,

    /// samples.tsv column used as the condition label.
    #[arg(long, default_value = "condition")]
    pub(super) condition_col: String,

    /// Keep only samples whose raw samples.tsv values satisfy every
    /// predicate (`col==v`, `col=v`, `col!=v`; repeatable or `;`-separated).
    #[arg(long, action = clap::ArgAction::Append)]
    pub(super) subset: Vec<String>,

    /// Relabel every condition not named in `--groups` to this label
    /// (which must itself appear in `--groups`).
    #[arg(long)]
    pub(super) collapse_others: Option<String>,

    /// Keep only samples whose samples.tsv row has a non-empty value in
    /// every listed column (comma-separated or repeatable), e.g.
    /// `--require-cols total_protein` to fit the same subjects as an
    /// absolute-scale run.
    #[arg(long, action = clap::ArgAction::Append)]
    pub(super) require_cols: Vec<String>,

    /// Like `--require-cols`, but the value must also parse as a finite
    /// number (drops placeholders such as `not measured`).
    #[arg(long, action = clap::ArgAction::Append)]
    pub(super) require_numeric: Vec<String>,

    /// Drop proteins whose fraction of missing samples within a comparison
    /// exceeds this value (1.0 = keep every protein).
    #[arg(long, default_value_t = 1.0)]
    pub(super) max_missing_fraction: f64,

    /// How protein groups sharing a gene symbol are reduced to one value per
    /// sample: `none` (lexically first assay), `mean` (mean of the observed
    /// assays), `max-observed` (assay observed in the most samples). Applies
    /// to paired-t, welch-t, ols, and mixed; the missingness filter runs on
    /// the collapsed gene.
    #[arg(long, value_enum, default_value_t = crate::design::CollapseGenes::None)]
    pub(super) collapse_genes: crate::design::CollapseGenes,
}

pub(super) fn apply_per_subject_proxy(args: &mut Args) -> Result<()> {
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
            DesignTerm::Condition | DesignTerm::Expression { .. } => false,
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
