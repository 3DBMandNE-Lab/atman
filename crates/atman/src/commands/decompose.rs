use anyhow::{bail, Context, Result};
use atman_core::bh_fdr;
use atman_core::compositional::{apply_transform, Transform};
use atman_core::decompose_unmix::{
    bootstrap_ci, ora_enrichment, select_k_auto, unmix, AbundanceCi, AbundanceMethod,
    AnnotationRow, EndmemberMethod, KSweepRow, LoadingCi, UnmixConfig, UnmixResult,
};
use atman_core::ica::{
    fast_ica, jaccard_top_n, pca_whiten, select_k_cumulative_variance, IcaResult,
};
use atman_core::nmf::{
    multi_seed_nmf, nmf as nmf_core, select_k as nmf_select_k, BetaLoss, Init, KSelection,
    MultiSeedConfig, NmfConfig, StabilityMetric as NmfStabilityMetric,
};
use atman_core::ica_null::{archetype_null, ArchetypeNullRow, NullMode, NullParams};
use atman_core::variance_decomposition::{decompose_archetype_variance, FixedFactor, VarianceRow};
use clap::{Args as ClapArgs, Subcommand, ValueEnum};
use csv::ReaderBuilder;
use serde_json::json;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{
    atomic_write, format_float, hash_canonical_inputs, read_measurements_long, read_samples,
    sidecar_path_for, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Multi-seed FastICA decomposition with seed-stability reporting.
    Ica(IcaArgs),
    /// Permutation-null calibration of archetype stability.
    Null(NullArgs),
    /// Partition archetype activation variance into fixed-effect,
    /// random-intercept, and residual components per archetype.
    Variance(VarianceArgs),
    /// Geometric compartmental unmixing (VCA endmember extraction
    /// + FCLS/UCLS abundance estimation).
    Unmix(UnmixArgs),
    /// Reconstruct the abundance matrix with one program's activation
    /// set to a target value (default 0). Produces a delta matrix
    /// showing each program's per-protein, per-sample contribution.
    Counterfactual(CounterfactualArgs),
    /// Multiplicative-updates NMF (Frobenius loss, Lee & Seung).
    Nmf(NmfArgs),
}

#[derive(ClapArgs, Debug)]
pub struct IcaArgs {
    /// Canonical Atman input directory (contains `measurements.tsv`).
    #[arg(long)]
    input_dir: PathBuf,

    /// Explicit number of components. Overrides `--k-selection` if set.
    #[arg(long)]
    k: Option<usize>,

    /// K-selection method; currently supports `cumulative-variance=<target>`.
    #[arg(long, default_value = "cumulative-variance=0.80")]
    k_selection: String,

    /// Minimum allowed K (lower clamp for `--k-selection`).
    #[arg(long, default_value_t = 3)]
    k_min: usize,

    /// Maximum allowed K (upper clamp for `--k-selection`).
    #[arg(long, default_value_t = 30)]
    k_max: usize,

    /// Number of seeds (must be >= 1; first seed is the reference).
    #[arg(long, default_value_t = 50)]
    n_seeds: usize,

    /// First seed value; subsequent seeds are sequential (seed, seed+1, ...).
    #[arg(long, default_value_t = 20260418)]
    seed: u64,

    /// Jaccard threshold for a program to be considered "recovered" in another seed.
    #[arg(long, default_value_t = 0.9)]
    seed_stability_threshold: f64,

    /// A program is flagged unstable when the fraction of alternative
    /// seeds that recover it (best-Jaccard >= `--seed-stability-threshold`)
    /// falls below this value. Default 0.9 — i.e. a program must recover
    /// in at least 90% of alternative seeds to pass.
    #[arg(long, default_value_t = 0.9)]
    min_stable_seed_fraction: f64,

    /// Stability metric; currently `jaccard-top20` (i.e. top-N |loading| overlap).
    #[arg(long, value_enum, default_value_t = StabilityMetric::JaccardTop20)]
    stability_metric: StabilityMetric,

    /// Top-N loadings used by the stability metric.
    #[arg(long, default_value_t = 20)]
    stability_top_n: usize,

    /// Maximum FastICA iterations per seed.
    #[arg(long, default_value_t = 300)]
    max_iter: usize,

    /// FastICA convergence tolerance.
    #[arg(long, default_value_t = 1e-4)]
    tol: f64,

    /// Drop assays where the fraction of missing samples exceeds this value. 0 = strict complete-case.
    #[arg(long, default_value_t = 0.0)]
    max_missing_fraction: f64,

    /// Imputation for residual missingness after the assay filter: `none` (fail) or `mean` (per-assay mean).
    #[arg(long, default_value = "none")]
    impute: String,

    /// Loadings output TSV path (`program\tassay_id\tgene_symbol\tloading`).
    #[arg(long)]
    output_loadings: PathBuf,

    /// Activations output TSV path (`cohort\tsubject_id\tsample_id\tprogram\tactivation`).
    #[arg(long)]
    output_activations: PathBuf,

    /// Stability output TSV path.
    #[arg(long)]
    output_stability: PathBuf,

    /// Cohort label written into the activations table (optional; defaults to directory name).
    #[arg(long)]
    cohort: Option<String>,

    /// Compositional transform applied to the sample × protein matrix
    /// before ICA. Atman's canonical abundance is already on a log
    /// scale, so these are linear operations on log values. `clr`
    /// (per-sample mean centering) is the recommended default for
    /// closed-sum proteomics. `alr` / `ratio-anchor` require
    /// `--alr-reference` (a gene symbol). `ilr` emits loadings in
    /// the Helmert coordinate basis rather than raw protein space.
    #[arg(long, value_enum, default_value_t = TransformArg::None)]
    transform: TransformArg,

    /// Gene symbol (matched against `samples` metadata `gene_symbol`)
    /// used as the reference for `--transform alr` or
    /// `--transform ratio-anchor`. Ignored for other transforms.
    #[arg(long)]
    alr_reference: Option<String>,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum TransformArg {
    #[value(name = "none")]
    None,
    #[value(name = "clr")]
    Clr,
    #[value(name = "alr")]
    Alr,
    #[value(name = "ilr")]
    Ilr,
    #[value(name = "ratio-anchor")]
    RatioAnchor,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum StabilityMetric {
    #[value(name = "jaccard-top20")]
    JaccardTop20,
}

#[derive(ClapArgs, Debug)]
pub struct CounterfactualArgs {
    /// Loadings TSV from `atman decompose ica` (columns: program, assay_id,
    /// gene_symbol, loading).
    #[arg(long)]
    loadings: PathBuf,

    /// Activations TSV from `atman decompose ica` (columns: cohort,
    /// subject_id, sample_id, program, activation).
    #[arg(long)]
    activations: PathBuf,

    /// Program to manipulate (e.g. `program_04`).
    #[arg(long)]
    set_program: String,

    /// Value to set the target program's activation to for all samples.
    #[arg(long, default_value_t = 0.0)]
    to: f64,

    /// Output TSV for the delta matrix (protein × sample contribution of
    /// the manipulated program). Wide format: assay_id, gene_symbol,
    /// then one column per sample_id.
    #[arg(long)]
    output: PathBuf,
}

#[derive(ClapArgs, Debug)]
pub struct NmfArgs {
    /// Canonical Atman input directory (contains `measurements.tsv`).
    #[arg(long)]
    pub input_dir: PathBuf,

    /// Number of NMF components. Required when `--k-selection fixed` (default).
    /// Ignored (with a warning) when `--k-selection` is automatic.
    #[arg(long)]
    pub k: Option<usize>,

    /// K-selection method: `fixed` (use `--k`), `cophenetic-knee`, or `rss-knee`.
    #[arg(long, default_value = "fixed")]
    pub k_selection: String,

    /// Minimum k for automatic k-selection sweep (≥ 2). Required when
    /// `--k-selection` is not `fixed`.
    #[arg(long)]
    pub k_min: Option<usize>,

    /// Maximum k for automatic k-selection sweep (≤ 50, > k_min). Required
    /// when `--k-selection` is not `fixed`.
    #[arg(long)]
    pub k_max: Option<usize>,

    /// Optional path for the k-sweep TSV when `--k-selection` is automatic.
    /// Columns: `k\tcophenetic\tmean_rss\tmean_kl\tselected`.
    #[arg(long)]
    pub output_k_sweep: Option<PathBuf>,

    /// `frobenius` or `kullback-leibler` (alias `kl`).
    #[arg(long, default_value = "frobenius")]
    pub beta_loss: String,

    /// Initialization strategy: `nndsvda` (deterministic) or `random`.
    #[arg(long, default_value = "nndsvda")]
    pub init: String,

    /// Update rule: only `mu` (multiplicative updates) is supported.
    #[arg(long, default_value = "mu")]
    pub solver: String,

    /// Maximum number of multiplicative-update iterations.
    #[arg(long, default_value_t = 400)]
    pub max_iter: usize,

    /// Convergence tolerance on the change in Frobenius error per iteration.
    #[arg(long, default_value_t = 1e-6)]
    pub tol: f64,

    /// RNG seed (used only when `--init random`).
    #[arg(long, default_value_t = 42)]
    pub seed: u64,

    /// Number of independent seeds to run. When 1 (default), single-seed NMF
    /// is used and no stability TSV is emitted.  When > 1, multi-seed NMF runs
    /// and programs are filtered by `--min-stable-seed-fraction`.
    #[arg(long, default_value_t = 1)]
    pub n_seeds: usize,

    /// Base seed for multi-seed runs; seed i = seed_base + i.
    #[arg(long, default_value_t = 42)]
    pub seed_base: u64,

    /// Minimum fraction of seeds in which a program must appear to survive
    /// filtering in multi-seed mode.
    #[arg(long, default_value_t = 0.9)]
    pub min_stable_seed_fraction: f64,

    /// Stability metric: `jaccard-top20` (Jaccard overlap of top-N loadings).
    #[arg(long, default_value = "jaccard-top20")]
    pub stability_metric: String,

    /// Top-N loadings used by the stability metric.
    #[arg(long, default_value_t = 20)]
    pub stability_top_n: usize,

    /// Optional path for the stability TSV (`program\tstable_seed_fraction\tn_seeds_present`).
    /// Only written when `--n-seeds > 1`.
    #[arg(long)]
    pub output_stability: Option<PathBuf>,

    /// Loadings output TSV path (`program\tassay_id\tgene_symbol\tloading`).
    #[arg(long)]
    pub output_loadings: PathBuf,

    /// Activations output TSV path (`sample_id\tprogram\tactivation`). Optional.
    #[arg(long)]
    pub output_activations: Option<PathBuf>,

    /// Pre-decomposition transform. Default `none` (NMF requires non-negative
    /// input; atman will reject the input loudly if any value is < 0).
    #[arg(long, default_value = "none")]
    pub transform: String,
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Ica(args) => run_ica(args),
        Command::Null(args) => run_null(args),
        Command::Variance(args) => run_variance(args),
        Command::Unmix(args) => run_unmix(args),
        Command::Counterfactual(args) => run_counterfactual(args),
        Command::Nmf(args) => nmf_run(args),
    }
}

fn nmf_run(args: NmfArgs) -> Result<()> {
    let started_at = SystemTime::now();

    // ── parse / validate discrete args ──────────────────────────────────────
    let beta_loss = match args.beta_loss.as_str() {
        "frobenius" => BetaLoss::Frobenius,
        "kullback-leibler" | "kl" => BetaLoss::KullbackLeibler,
        other => bail!(
            "--beta-loss {:?}: expected `frobenius` or `kullback-leibler`",
            other
        ),
    };

    let init = match args.init.as_str() {
        "nndsvda" => Init::Nndsvda,
        "random" => Init::Random,
        other => bail!(
            "--init {:?}: expected `nndsvda` or `random`",
            other
        ),
    };

    match args.solver.as_str() {
        "mu" => {}
        other => bail!(
            "--solver {:?}: only `mu` (multiplicative updates) is supported in this release",
            other
        ),
    }

    // ── parse k-selection ────────────────────────────────────────────────────
    let k_sel = match args.k_selection.as_str() {
        "fixed" => KSelection::Fixed,
        "cophenetic-knee" => KSelection::CopheneticKnee,
        "rss-knee" => KSelection::RssKnee,
        other => bail!(
            "--k-selection {:?}: expected `fixed`, `cophenetic-knee`, or `rss-knee`",
            other
        ),
    };

    // Validate --k / --k-min / --k-max depending on selection mode.
    let k_fixed: usize = if k_sel == KSelection::Fixed {
        let k = args.k.ok_or_else(|| {
            anyhow::anyhow!("--k is required when --k-selection fixed (the default)")
        })?;
        if k == 0 {
            bail!("--k must be at least 1");
        }
        k
    } else {
        0 // placeholder; will be determined by select_k
    };

    let (k_min, k_max) = if k_sel != KSelection::Fixed {
        let k_min = args.k_min.ok_or_else(|| {
            anyhow::anyhow!("--k-min is required when --k-selection is not fixed")
        })?;
        let k_max = args.k_max.ok_or_else(|| {
            anyhow::anyhow!("--k-max is required when --k-selection is not fixed")
        })?;
        if k_min < 2 {
            bail!("--k-min must be >= 2, got {}", k_min);
        }
        if k_max > 50 {
            bail!("--k-max must be <= 50, got {}", k_max);
        }
        if k_max <= k_min {
            bail!("--k-max ({}) must be > --k-min ({})", k_max, k_min);
        }
        (k_min, k_max)
    } else {
        (0, 0) // unused for Fixed
    };

    // ── load matrix (canonical input: measurements.tsv) ────────────────────
    let tsv = args.input_dir.join("measurements.tsv");
    let records = read_measurements_long(&tsv)?;
    if records.is_empty() {
        bail!("no measurements in {:?}", tsv);
    }

    // Pivot to sample × feature (same pattern as load_matrix for ICA).
    let mut abundance_by_key: BTreeMap<(String, String), f64> = BTreeMap::new();
    let mut sample_order: Vec<String> = Vec::new();
    let mut seen_samples: BTreeSet<String> = BTreeSet::new();
    let mut assay_meta: BTreeMap<String, String> = BTreeMap::new(); // assay_id → gene_symbol

    for r in &records {
        if r.dropped_by_qc {
            continue;
        }
        let assay = r.assay_id.0.clone();
        let sample = r.sample_id.clone();
        let abundance = r.abundance.as_f64();
        if !abundance.is_finite() {
            continue;
        }
        if seen_samples.insert(sample.clone()) {
            sample_order.push(sample.clone());
        }
        assay_meta
            .entry(assay.clone())
            .or_insert_with(|| r.gene_symbol.clone().unwrap_or_default());
        abundance_by_key.insert((assay, sample), abundance);
    }

    let n_samples = sample_order.len();
    if n_samples < 2 {
        bail!("need at least 2 samples for NMF, got {}", n_samples);
    }

    let assay_order: Vec<String> = assay_meta.keys().cloned().collect();
    let n_assays = assay_order.len();

    if k_sel == KSelection::Fixed && k_fixed > n_samples.min(n_assays) {
        bail!(
            "--k={} exceeds min(n_samples={}, n_assays={})",
            k_fixed,
            n_samples,
            n_assays
        );
    }

    // Build dense matrix: data[i][j] = sample i, assay j.
    // NMF requires no missing values; fail loudly if any are absent.
    let mut data = vec![vec![0.0_f64; n_assays]; n_samples];
    for (i, sample) in sample_order.iter().enumerate() {
        for (j, assay) in assay_order.iter().enumerate() {
            match abundance_by_key.get(&(assay.clone(), sample.clone())) {
                Some(&v) => data[i][j] = v,
                None => bail!(
                    "NMF: missing value for sample={} assay={}; NMF requires a complete matrix. \
                     Drop incomplete assays upstream or impute before running decompose nmf.",
                    sample, assay
                ),
            }
        }
    }

    // ── apply transform (only "none" supported; validate non-negativity) ────
    match args.transform.as_str() {
        "none" => {
            // Validate that all values are non-negative.
            for (i, row) in data.iter().enumerate() {
                for (j, &v) in row.iter().enumerate() {
                    if v < 0.0 {
                        bail!(
                            "NMF (--transform none) requires non-negative input, but \
                             sample={} assay={} has value {}. \
                             Apply a non-negativity-preserving transform upstream \
                             (e.g. shift to [0, ∞)) before calling decompose nmf.",
                            sample_order[i], assay_order[j], v
                        );
                    }
                }
            }
        }
        other => bail!(
            "--transform {:?}: only `none` is supported for decompose nmf in this release. \
             CLR/log are not appropriate without thought; ensure your input is already non-negative.",
            other
        ),
    }

    // ── validate multi-seed args ────────────────────────────────────────────
    if args.n_seeds == 0 {
        bail!("--n-seeds must be at least 1");
    }
    let nmf_stability_metric = match args.stability_metric.as_str() {
        "jaccard-top20" | "jaccard-topn" | "jaccard" => NmfStabilityMetric::JaccardTopN,
        other => bail!(
            "--stability-metric {:?}: expected `jaccard-top20`",
            other
        ),
    };

    // ── build NmfConfig base (k will be finalised after k-selection) ─────────
    let nmf_cfg_base = NmfConfig {
        k: if k_sel == KSelection::Fixed { k_fixed } else { k_min }, // temp; overridden below
        beta_loss,
        init,
        max_iter: args.max_iter,
        tol: args.tol,
        seed: if args.n_seeds == 1 { args.seed } else { args.seed_base },
    };

    // ── run k-selection (or Fixed pass-through) ───────────────────────────────
    let ms_cfg_for_sel = MultiSeedConfig {
        n_seeds: args.n_seeds,
        seed_base: args.seed_base,
        stability_metric: nmf_stability_metric,
        stability_top_n: args.stability_top_n,
    };

    let (effective_k, k_sweep_result) = if k_sel == KSelection::Fixed {
        (k_fixed, None)
    } else {
        eprintln!(
            "decompose nmf: k-selection={} k_min={} k_max={} n_seeds={}",
            args.k_selection, k_min, k_max, args.n_seeds,
        );
        let sel_result = nmf_select_k(&data, &nmf_cfg_base, &ms_cfg_for_sel, k_min, k_max, k_sel);
        eprintln!(
            "decompose nmf: selected_k={}",
            sel_result.selected_k,
        );
        let sk = sel_result.selected_k;
        (sk, Some(sel_result))
    };

    let nmf_cfg = NmfConfig {
        k: effective_k,
        ..nmf_cfg_base
    };

    eprintln!(
        "decompose nmf: n_samples={} n_assays={} k={} init={} max_iter={} tol={} n_seeds={}",
        n_samples, n_assays, effective_k, args.init, args.max_iter, args.tol, args.n_seeds,
    );

    // ── run NMF (single-seed or multi-seed) ─────────────────────────────────
    //
    // Outputs produced differ by mode:
    //   n_seeds == 1: run nmf_core once; H is k×p loadings, W is n×k activations.
    //   n_seeds >  1: run multi_seed_nmf; filter by min_stable_seed_fraction;
    //                 re-fit activations from surviving mean loadings via a
    //                 fixed-H NMF (1 MU update step, H held fixed).

    // final_h: Vec of (program_label, Vec<f64> of length p)
    // final_w: Vec of (sample_idx, Vec<f64> of length k_surviving)  [row order = sample_order]
    // stability_rows: Some only in multi-seed mode.

    struct SingleResult {
        h: Vec<Vec<f64>>,      // k × p
        w: Vec<Vec<f64>>,      // n × k
        n_iter: usize,
        converged: bool,
        final_error: f64,
    }

    enum RunOutput {
        Single(SingleResult),
        MultiSeed {
            /// All programs (pre-filter) with stability info.
            all_rows: Vec<atman_core::nmf::StableProgramRow>,
            /// Indices into all_rows that passed the stability filter.
            surviving: Vec<usize>,
            /// Re-fitted activations: n × k_surviving.
            w_refit: Vec<Vec<f64>>,
        },
    }

    let run_output = if args.n_seeds == 1 {
        let result = nmf_core(&data, &nmf_cfg);
        eprintln!(
            "decompose nmf: n_iter={} converged={} final_error={:.6e}",
            result.n_iter, result.converged, result.final_error
        );
        RunOutput::Single(SingleResult {
            h: result.h,
            w: result.w,
            n_iter: result.n_iter,
            converged: result.converged,
            final_error: result.final_error,
        })
    } else {
        let ms_cfg = MultiSeedConfig {
            n_seeds: args.n_seeds,
            seed_base: args.seed_base,
            stability_metric: nmf_stability_metric,
            stability_top_n: args.stability_top_n,
        };
        eprintln!(
            "decompose nmf: running {} seeds (seed_base={}) …",
            args.n_seeds, args.seed_base
        );
        let all_rows = multi_seed_nmf(&data, &nmf_cfg, &ms_cfg);
        let surviving: Vec<usize> = all_rows
            .iter()
            .enumerate()
            .filter(|(_, r)| r.stable_seed_fraction >= args.min_stable_seed_fraction)
            .map(|(i, _)| i)
            .collect();
        eprintln!(
            "decompose nmf: {} / {} programs pass stability >= {}",
            surviving.len(),
            all_rows.len(),
            args.min_stable_seed_fraction
        );

        // Re-fit W for surviving programs with H fixed (NNLS-style via MU, H locked).
        // Strategy: run NMF on data with k = k_surviving, initialised so H = mean_loadings
        // (surviving) and W = uniform.  We freeze H by using only the W-update MU rule.
        // For simplicity and reproducibility we use atman_core::nmf's existing infrastructure:
        // we run a Frobenius MU update with H held to mean_loadings by re-initialising H
        // after each W update. 50 iterations is enough for activations to converge.
        let k_surv = surviving.len();
        let w_refit = if k_surv == 0 {
            vec![]
        } else {
            // Build fixed H from mean loadings of surviving programs.
            let h_fixed: Vec<Vec<f64>> = surviving
                .iter()
                .map(|&i| all_rows[i].mean_loading.clone())
                .collect();

            // Initialise W randomly (n × k_surv), then iterate W-only MU updates.
            let mean_x: f64 = {
                let total: f64 = data.iter().flat_map(|r| r.iter()).sum();
                let cnt = (n_samples * n_assays) as f64;
                if cnt > 0.0 { total / cnt } else { 1.0 }
            };
            let scale = (mean_x / k_surv.max(1) as f64).sqrt().max(1e-9);
            let mut rng = atman_core::ica::Xoshiro256pp::new(args.seed_base.wrapping_add(9999));
            let mut w: Vec<Vec<f64>> = (0..n_samples)
                .map(|_| (0..k_surv).map(|_| rng.next_normal().abs() * scale).collect())
                .collect();

            const EPS: f64 = 1e-9;
            for _iter in 0..100 {
                // W update: W[i,a] <- W[i,a] * (X H^T)[i,a] / (W H H^T)[i,a]
                // XHt: n × k_surv
                let mut xht = vec![vec![0.0_f64; k_surv]; n_samples];
                for i in 0..n_samples {
                    for a in 0..k_surv {
                        for j in 0..n_assays {
                            xht[i][a] += data[i][j] * h_fixed[a][j];
                        }
                    }
                }
                // HHt: k_surv × k_surv
                let mut hht = vec![vec![0.0_f64; k_surv]; k_surv];
                for a in 0..k_surv {
                    for b in 0..k_surv {
                        for j in 0..n_assays {
                            hht[a][b] += h_fixed[a][j] * h_fixed[b][j];
                        }
                    }
                }
                // WHHt: n × k_surv
                let mut whht = vec![vec![0.0_f64; k_surv]; n_samples];
                for i in 0..n_samples {
                    for a in 0..k_surv {
                        for b in 0..k_surv {
                            whht[i][a] += w[i][b] * hht[b][a];
                        }
                    }
                }
                for i in 0..n_samples {
                    for a in 0..k_surv {
                        let num = xht[i][a];
                        let den = whht[i][a] + EPS;
                        w[i][a] *= num / den;
                        if w[i][a] < 0.0 { w[i][a] = 0.0; }
                    }
                }
            }
            w
        };

        RunOutput::MultiSeed { all_rows, surviving, w_refit }
    };

    // ── write loadings TSV ───────────────────────────────────────────────────
    if let Some(parent) = args.output_loadings.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output dir {:?}", parent))?;
    }
    {
        let mut out = String::from("program\tassay_id\tgene_symbol\tloading\n");
        match &run_output {
            RunOutput::Single(r) => {
                for (a, h_row) in r.h.iter().enumerate() {
                    let prog = program_name(a);
                    for (j, &load) in h_row.iter().enumerate() {
                        let assay_id = &assay_order[j];
                        let gene_symbol =
                            assay_meta.get(assay_id).map(|s| s.as_str()).unwrap_or("");
                        out.push_str(&prog);
                        out.push('\t');
                        out.push_str(assay_id);
                        out.push('\t');
                        out.push_str(gene_symbol);
                        out.push('\t');
                        out.push_str(&format_float(load));
                        out.push('\n');
                    }
                }
            }
            RunOutput::MultiSeed { all_rows, surviving, .. } => {
                for (out_idx, &surv_idx) in surviving.iter().enumerate() {
                    let row = &all_rows[surv_idx];
                    let prog = program_name(out_idx);
                    for (j, &load) in row.mean_loading.iter().enumerate() {
                        let assay_id = &assay_order[j];
                        let gene_symbol =
                            assay_meta.get(assay_id).map(|s| s.as_str()).unwrap_or("");
                        out.push_str(&prog);
                        out.push('\t');
                        out.push_str(assay_id);
                        out.push('\t');
                        out.push_str(gene_symbol);
                        out.push('\t');
                        out.push_str(&format_float(load));
                        out.push('\n');
                    }
                }
            }
        }
        atomic_write(&args.output_loadings, out.as_bytes())?;
    }
    eprintln!("decompose nmf: loadings={}", args.output_loadings.display());

    // ── write activations TSV ────────────────────────────────────────────────
    let mut output_files: Vec<PathBuf> = vec![args.output_loadings.clone()];
    if let Some(ref act_path) = args.output_activations {
        if let Some(parent) = act_path.parent() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("creating output dir {:?}", parent))?;
        }
        let mut out = String::from("sample_id\tprogram\tactivation\n");
        match &run_output {
            RunOutput::Single(r) => {
                for (i, w_row) in r.w.iter().enumerate() {
                    let sample_id = &sample_order[i];
                    for (a, &activation) in w_row.iter().enumerate() {
                        let prog = program_name(a);
                        out.push_str(sample_id);
                        out.push('\t');
                        out.push_str(&prog);
                        out.push('\t');
                        out.push_str(&format_float(activation));
                        out.push('\n');
                    }
                }
            }
            RunOutput::MultiSeed { w_refit, .. } => {
                for (i, w_row) in w_refit.iter().enumerate() {
                    let sample_id = &sample_order[i];
                    for (a, &activation) in w_row.iter().enumerate() {
                        let prog = program_name(a);
                        out.push_str(sample_id);
                        out.push('\t');
                        out.push_str(&prog);
                        out.push('\t');
                        out.push_str(&format_float(activation));
                        out.push('\n');
                    }
                }
            }
        }
        atomic_write(act_path, out.as_bytes())?;
        eprintln!("decompose nmf: activations={}", act_path.display());
        output_files.push(act_path.clone());
    }

    // ── write stability TSV (multi-seed only) ────────────────────────────────
    if let RunOutput::MultiSeed { ref all_rows, ref surviving, .. } = run_output {
        if let Some(ref stab_path) = args.output_stability {
            if let Some(parent) = stab_path.parent() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("creating stability dir {:?}", parent))?;
            }
            let mut out =
                String::from("program\tstable_seed_fraction\tn_seeds_present\n");
            for row in all_rows {
                let n_seeds_present =
                    (row.stable_seed_fraction * args.n_seeds as f64).round() as usize;
                out.push_str(&row.program_id);
                out.push('\t');
                out.push_str(&format_float(row.stable_seed_fraction));
                out.push('\t');
                out.push_str(&n_seeds_present.to_string());
                out.push('\n');
            }
            let _ = surviving; // already written in loadings
            atomic_write(stab_path, out.as_bytes())?;
            eprintln!("decompose nmf: stability={}", stab_path.display());
            output_files.push(stab_path.clone());
        }
    }

    // ── write k-sweep TSV (automatic k-selection only) ───────────────────────
    if let Some(ref sweep_result) = k_sweep_result {
        if let Some(ref sweep_path) = args.output_k_sweep {
            if let Some(parent) = sweep_path.parent() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("creating k-sweep dir {:?}", parent))?;
            }
            let mut out = String::from("k\tcophenetic\tmean_rss\tmean_kl\tselected\n");
            for row in &sweep_result.per_k {
                let cophenetic_str = if row.cophenetic.is_nan() {
                    "NA".to_string()
                } else {
                    format_float(row.cophenetic)
                };
                let mean_kl_str = row
                    .mean_kl
                    .map(|v| format_float(v))
                    .unwrap_or_else(|| "NA".to_string());
                let selected = if row.k == sweep_result.selected_k { "1" } else { "0" };
                out.push_str(&row.k.to_string());
                out.push('\t');
                out.push_str(&cophenetic_str);
                out.push('\t');
                out.push_str(&format_float(row.mean_rss));
                out.push('\t');
                out.push_str(&mean_kl_str);
                out.push('\t');
                out.push_str(selected);
                out.push('\n');
            }
            atomic_write(sweep_path, out.as_bytes())?;
            eprintln!("decompose nmf: k-sweep={}", sweep_path.display());
            output_files.push(sweep_path.clone());
        }
    }

    // ── sidecar ─────────────────────────────────────────────────────────────
    let finished_at = SystemTime::now();
    let canonical_inputs: &[&str] = &["measurements.tsv", "samples.tsv", "proteins.tsv"];
    let inputs_sha256 = hash_canonical_inputs(&args.input_dir, canonical_inputs)?;
    let sidecar = sidecar_path_for(&args.output_loadings);

    let (n_iter_val, converged_val, final_error_val) = match &run_output {
        RunOutput::Single(r) => (
            serde_json::Value::from(r.n_iter),
            serde_json::Value::from(r.converged),
            serde_json::Value::from(r.final_error),
        ),
        RunOutput::MultiSeed { .. } => (
            serde_json::Value::Null,
            serde_json::Value::Null,
            serde_json::Value::Null,
        ),
    };

    let selected_k_val = k_sweep_result
        .as_ref()
        .map(|r| serde_json::Value::from(r.selected_k))
        .unwrap_or(serde_json::Value::Null);

    write_run_sidecar(
        &sidecar,
        "decompose nmf",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "k": effective_k,
            "k-selection": args.k_selection,
            "k-min": args.k_min,
            "k-max": args.k_max,
            "selected-k": selected_k_val,
            "beta-loss": args.beta_loss,
            "init": args.init,
            "solver": args.solver,
            "max-iter": args.max_iter,
            "tol": args.tol,
            "seed": args.seed,
            "n-seeds": args.n_seeds,
            "seed-base": args.seed_base,
            "min-stable-seed-fraction": args.min_stable_seed_fraction,
            "stability-metric": args.stability_metric,
            "stability-top-n": args.stability_top_n,
            "transform": args.transform,
            "output-loadings": args.output_loadings.display().to_string(),
            "output-activations": args.output_activations.as_ref().map(|p| p.display().to_string()),
            "output-stability": args.output_stability.as_ref().map(|p| p.display().to_string()),
            "output-k-sweep": args.output_k_sweep.as_ref().map(|p| p.display().to_string()),
            "decomposition_method": "nmf",
            "n_iter": n_iter_val,
            "converged": converged_val,
            "final_error": final_error_val,
        }),
        &inputs_sha256,
        &output_files,
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("decompose nmf: sidecar={}", sidecar.display());
    Ok(())
}

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

fn run_variance(args: VarianceArgs) -> Result<()> {
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

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum NullModeArg {
    #[value(name = "sample-shuffle")]
    SampleShuffle,
    #[value(name = "protein-shuffle")]
    ProteinShuffle,
    #[value(name = "gaussian-matched")]
    GaussianMatched,
}

impl NullModeArg {
    fn to_core(self) -> NullMode {
        match self {
            Self::SampleShuffle => NullMode::SampleShuffle,
            Self::ProteinShuffle => NullMode::ProteinShuffle,
            Self::GaussianMatched => NullMode::GaussianMatched,
        }
    }
}

#[derive(ClapArgs, Debug)]
pub struct NullArgs {
    /// Canonical Atman input directory (contains `measurements.tsv`).
    #[arg(long)]
    input_dir: PathBuf,

    /// Number of ICA components to fit. Required — no k-selection here
    /// so the null run exactly mirrors the K chosen in `decompose ica`.
    #[arg(long)]
    k: usize,

    /// Number of null permutation iterations. Higher = tighter p-value
    /// resolution. Floor on p is `1/(n_perm + 1)` so use ≥ 200 for
    /// `q < 0.05` claims to mean anything.
    #[arg(long, default_value_t = 200)]
    n_perm: usize,

    /// FastICA seeds per matrix (real and each null). Seed 0 is the
    /// reference; stability is mean best-Jaccard across the remaining
    /// `n_seeds − 1` alternative seeds. Must be ≥ 2.
    #[arg(long, default_value_t = 3)]
    n_seeds: usize,

    /// Top-level seed. Each null iteration's sub-seed derives from
    /// `(seed, iter)` deterministically.
    #[arg(long, default_value_t = 20260418)]
    seed: u64,

    /// Null-matrix construction mode.
    #[arg(long, value_enum, default_value_t = NullModeArg::ProteinShuffle)]
    null_mode: NullModeArg,

    /// Top-N loadings used by the Jaccard stability metric.
    #[arg(long, default_value_t = 20)]
    top_n: usize,

    /// FastICA max iterations per seed.
    #[arg(long, default_value_t = 300)]
    max_iter: usize,

    /// FastICA convergence tolerance.
    #[arg(long, default_value_t = 1e-4)]
    tol: f64,

    /// Drop assays with more than this fraction of missing samples.
    #[arg(long, default_value_t = 0.0)]
    max_missing_fraction: f64,

    /// Residual missingness imputation: `none` fails, `mean` uses
    /// per-assay mean.
    #[arg(long, default_value = "none")]
    impute: String,

    /// BH-q threshold used by the output `decision` column
    /// (`signal` when `null_q < threshold`, else `noise`).
    #[arg(long, default_value_t = 0.05)]
    q_threshold: f64,

    /// Output TSV: one row per program with observed_stability,
    /// null distribution summary, null_p, null_q, decision.
    #[arg(long)]
    output: PathBuf,
}

fn run_null(args: NullArgs) -> Result<()> {
    let started_at = SystemTime::now();
    if args.k == 0 {
        bail!("--k must be at least 1");
    }
    if args.n_perm == 0 {
        bail!("--n-perm must be at least 1");
    }
    if args.n_seeds < 2 {
        bail!("--n-seeds must be >= 2 for stability calibration");
    }
    if args.top_n == 0 {
        bail!("--top-n must be at least 1");
    }
    if !(0.0..=1.0).contains(&args.q_threshold) {
        bail!("--q-threshold must be in [0, 1]");
    }

    // Reuse the abundance-matrix loader from `run_ica` by reconstructing
    // an equivalent IcaArgs; `load_matrix` only reads the fields below.
    let adapter = IcaArgs {
        input_dir: args.input_dir.clone(),
        k: Some(args.k),
        k_selection: "cumulative-variance=0.80".to_string(),
        k_min: 1,
        k_max: args.k.max(1),
        n_seeds: args.n_seeds,
        seed: args.seed,
        seed_stability_threshold: 0.9,
        min_stable_seed_fraction: 0.9,
        stability_metric: StabilityMetric::JaccardTop20,
        stability_top_n: args.top_n,
        max_iter: args.max_iter,
        tol: args.tol,
        max_missing_fraction: args.max_missing_fraction,
        impute: args.impute.clone(),
        output_loadings: PathBuf::new(),
        output_activations: PathBuf::new(),
        output_stability: PathBuf::new(),
        cohort: None,
        transform: TransformArg::None,
        alr_reference: None,
    };
    let matrix = load_matrix(&adapter)?;
    if args.k > matrix.samples.len() || args.k > matrix.assays.len() {
        bail!(
            "k={} exceeds min(n_samples={}, n_assays={})",
            args.k,
            matrix.samples.len(),
            matrix.assays.len()
        );
    }

    eprintln!(
        "decompose null: n_samples={} n_assays={} k={} n_perm={} n_seeds={} mode={}",
        matrix.samples.len(),
        matrix.assays.len(),
        args.k,
        args.n_perm,
        args.n_seeds,
        args.null_mode.to_core().as_str()
    );

    let params = NullParams {
        k: args.k,
        n_perm: args.n_perm,
        n_seeds: args.n_seeds,
        seed: args.seed,
        mode: args.null_mode.to_core(),
        top_n: args.top_n,
        max_iter: args.max_iter,
        tol: args.tol,
    };
    let rows: Vec<ArchetypeNullRow> =
        archetype_null(&matrix.data, &params).map_err(|e| anyhow::anyhow!(e))?;

    // BH on null_p across programs in this run.
    let ps: Vec<Option<f64>> = rows.iter().map(|r| Some(r.null_p)).collect();
    let qs = bh_fdr(&ps);

    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output dir {:?}", parent))?;
    }
    let mut out = String::from(
        "program\tobserved_stability\tnull_stability_mean\tnull_stability_p95\tnull_p\tnull_q\tdecision\n",
    );
    for (row, q) in rows.iter().zip(qs.iter()) {
        let program = program_name(row.program - 1);
        let q_val = q.unwrap_or(f64::NAN);
        let decision = if q_val.is_finite() && q_val < args.q_threshold {
            "signal"
        } else {
            "noise"
        };
        out.push_str(&format!(
            "{program}\t{}\t{}\t{}\t{}\t{}\t{decision}\n",
            format_float(row.observed_stability),
            format_float(row.null_stability_mean),
            format_float(row.null_stability_p95),
            format_float(row.null_p),
            format_float(q_val),
        ));
    }
    atomic_write(&args.output, out.as_bytes())?;

    eprintln!(
        "decompose null: wrote {} programs to {}",
        rows.len(),
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let canonical_inputs: &[&str] = &["measurements.tsv", "samples.tsv", "proteins.tsv"];
    let inputs_sha256 = hash_canonical_inputs(&args.input_dir, canonical_inputs)?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "decompose null",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output": args.output.display().to_string(),
            "k": args.k,
            "n-perm": args.n_perm,
            "n-seeds": args.n_seeds,
            "seed": args.seed,
            "null-mode": args.null_mode.to_core().as_str(),
            "top-n": args.top_n,
            "max-iter": args.max_iter,
            "tol": args.tol,
            "max-missing-fraction": args.max_missing_fraction,
            "impute": args.impute,
            "q-threshold": args.q_threshold,
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("decompose null: sidecar={}", sidecar.display());
    Ok(())
}

fn run_ica(args: IcaArgs) -> Result<()> {
    let started_at = SystemTime::now();
    if args.n_seeds == 0 {
        bail!("--n-seeds must be at least 1");
    }
    if args.k_min == 0 {
        bail!("--k-min must be at least 1");
    }
    if args.k_max < args.k_min {
        bail!("--k-max must be >= --k-min");
    }
    if !(0.0..=1.0).contains(&args.seed_stability_threshold) {
        bail!("--seed-stability-threshold must be in [0, 1]");
    }
    if !(0.0..=1.0).contains(&args.min_stable_seed_fraction) {
        bail!("--min-stable-seed-fraction must be in [0, 1]");
    }
    if args.stability_top_n == 0 {
        bail!("--stability-top-n must be at least 1");
    }

    let matrix = load_matrix(&args)?;
    let (matrix, transform_meta) = apply_compositional(matrix, &args)?;
    let k = resolve_k(&matrix, &args)?;
    if k > matrix.samples.len() || k > matrix.assays.len() {
        bail!(
            "k={k} exceeds min(n_samples={}, n_assays={})",
            matrix.samples.len(),
            matrix.assays.len()
        );
    }

    eprintln!(
        "decompose ica: n_samples={} n_assays={} k={} n_seeds={} threshold={} transform={}",
        matrix.samples.len(),
        matrix.assays.len(),
        k,
        args.n_seeds,
        args.seed_stability_threshold,
        transform_meta.name,
    );

    let mut runs = Vec::with_capacity(args.n_seeds);
    for offset in 0..args.n_seeds {
        let seed = args.seed.wrapping_add(offset as u64);
        let result = fast_ica(&matrix.data, k, seed, args.max_iter, args.tol);
        runs.push((seed, result));
    }

    let (ref_seed, ref_run) = runs.first().expect("at least one seed");
    let reference = canonicalize(ref_run);
    let cohort = args
        .cohort
        .clone()
        .unwrap_or_else(|| derive_cohort(&args.input_dir));

    write_loadings(&args.output_loadings, &matrix, &reference)?;
    write_activations(&args.output_activations, &matrix, &reference, &cohort)?;
    write_stability(
        &args.output_stability,
        &reference,
        &runs[1..],
        args.stability_top_n,
        args.seed_stability_threshold,
        args.min_stable_seed_fraction,
        *ref_seed,
    )?;

    eprintln!(
        "decompose ica: reference_seed={} loadings={} activations={} stability={}",
        ref_seed,
        args.output_loadings.display(),
        args.output_activations.display(),
        args.output_stability.display()
    );

    // Emit `transform_applied.json` next to the loadings whenever a
    // non-`none` transform ran, so downstream tools can audit which
    // coordinate system the archetypes live in.
    let transform_applied_path = if transform_meta.name != "none" {
        let p = args
            .output_loadings
            .parent()
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| PathBuf::from("."))
            .join("transform_applied.json");
        write_transform_applied(&p, &transform_meta)?;
        eprintln!("decompose ica: transform_applied={}", p.display());
        Some(p)
    } else {
        None
    };

    let finished_at = SystemTime::now();
    let canonical_inputs: &[&str] = &["measurements.tsv", "samples.tsv", "proteins.tsv"];
    let inputs_sha256 = hash_canonical_inputs(&args.input_dir, canonical_inputs)?;
    let sidecar = sidecar_path_for(&args.output_loadings);
    let stability_metric = match args.stability_metric {
        StabilityMetric::JaccardTop20 => "jaccard-top20",
    };
    let mut outputs = vec![
        args.output_loadings.clone(),
        args.output_activations.clone(),
        args.output_stability.clone(),
    ];
    if let Some(p) = &transform_applied_path {
        outputs.push(p.clone());
    }
    write_run_sidecar(
        &sidecar,
        "decompose ica",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "k": args.k,
            "k-selection": args.k_selection,
            "k-min": args.k_min,
            "k-max": args.k_max,
            "n-seeds": args.n_seeds,
            "seed": args.seed,
            "seed-stability-threshold": args.seed_stability_threshold,
            "min-stable-seed-fraction": args.min_stable_seed_fraction,
            "stability-metric": stability_metric,
            "stability-top-n": args.stability_top_n,
            "max-iter": args.max_iter,
            "tol": args.tol,
            "max-missing-fraction": args.max_missing_fraction,
            "impute": args.impute,
            "output-loadings": args.output_loadings.display().to_string(),
            "output-activations": args.output_activations.display().to_string(),
            "output-stability": args.output_stability.display().to_string(),
            "cohort": args.cohort,
            "transform": transform_meta.name,
            "alr-reference": transform_meta.alr_reference,
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        {
            // Resolved-value fields: record the integer `k` actually used
            // so a reviewer reading the sidecar doesn't have to count
            // columns in loadings.tsv to find out which `k` the
            // `--k-selection` rule picked. `k_resolution_source` names
            // the path: `cli` when --k was passed, `<rule>` otherwise.
            let mut extras = serde_json::Map::new();
            extras.insert("k_resolved".into(), serde_json::json!(k));
            let source = if args.k.is_some() {
                "cli".to_string()
            } else {
                args.k_selection.clone()
            };
            extras.insert("k_resolution_source".into(), serde_json::json!(source));
            Some(extras)
        },
    )?;
    eprintln!("decompose ica: sidecar={}", sidecar.display());
    Ok(())
}

fn resolve_k(matrix: &AbundanceMatrix, args: &IcaArgs) -> Result<usize> {
    if let Some(k) = args.k {
        if k == 0 {
            bail!("--k must be at least 1");
        }
        return Ok(k);
    }
    let target = parse_k_selection(&args.k_selection)?;
    let whitening = pca_whiten(&matrix.data, matrix.samples.len().min(matrix.assays.len()));
    let chosen = select_k_cumulative_variance(
        &whitening.full_spectrum,
        target,
        args.k_min,
        args.k_max
            .min(matrix.samples.len().saturating_sub(1))
            .max(1),
    );
    Ok(chosen)
}

fn parse_k_selection(spec: &str) -> Result<f64> {
    let trimmed = spec.trim();
    if let Some(rest) = trimmed.strip_prefix("cumulative-variance=") {
        let value: f64 = rest
            .parse()
            .with_context(|| format!("parsing k-selection target in {:?}", spec))?;
        if !(0.0..=1.0).contains(&value) {
            bail!("cumulative-variance target must be in (0, 1]");
        }
        Ok(value)
    } else {
        bail!(
            "unsupported --k-selection {:?}; use cumulative-variance=<target>",
            spec
        )
    }
}

fn derive_cohort(path: &Path) -> String {
    path.file_name()
        .and_then(|os| os.to_str())
        .unwrap_or("cohort")
        .to_string()
}

#[derive(Debug)]
struct AbundanceMatrix {
    samples: Vec<SampleInfo>,
    assays: Vec<AssayInfo>,
    data: Vec<Vec<f64>>,
}

#[derive(Debug, Clone)]
struct SampleInfo {
    sample_id: String,
    subject_id: String,
}

#[derive(Debug, Clone)]
struct AssayInfo {
    assay_id: String,
    gene_symbol: String,
}

fn load_matrix(args: &IcaArgs) -> Result<AbundanceMatrix> {
    let tsv = args.input_dir.join("measurements.tsv");
    let records = read_measurements_long(&tsv)?;
    if records.is_empty() {
        bail!("no measurements in {:?}", tsv);
    }

    let mut abundance_by_key: BTreeMap<(String, String), (f64, Option<String>)> = BTreeMap::new();
    let mut sample_order: Vec<String> = Vec::new();
    let mut seen_samples: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    let mut assay_meta: BTreeMap<String, Option<String>> = BTreeMap::new();
    for r in &records {
        if r.dropped_by_qc {
            continue;
        }
        let assay = r.assay_id.0.clone();
        let sample = r.sample_id.clone();
        let abundance = r.abundance.as_f64();
        if !abundance.is_finite() {
            continue;
        }
        if seen_samples.insert(sample.clone()) {
            sample_order.push(sample.clone());
        }
        assay_meta
            .entry(assay.clone())
            .or_insert_with(|| r.gene_symbol.clone());
        abundance_by_key.insert((assay, sample), (abundance, r.gene_symbol.clone()));
    }

    // Load samples.tsv for subject IDs (if present).
    let samples_path = args.input_dir.join("samples.tsv");
    let subject_lookup: BTreeMap<String, String> = if samples_path.exists() {
        crate::io::read_samples(&samples_path)?
            .into_iter()
            .map(|s| (s.sample_id.clone(), s.subject_id.unwrap_or_default()))
            .collect()
    } else {
        BTreeMap::new()
    };

    if !(0.0..=1.0).contains(&args.max_missing_fraction) {
        bail!("--max-missing-fraction must be in [0, 1]");
    }
    let impute_mean = match args.impute.as_str() {
        "none" => false,
        "mean" => true,
        other => bail!("--impute {:?}; expected none or mean", other),
    };

    // Drop assays whose missing-sample fraction exceeds the configured threshold.
    let n_samples = sample_order.len();
    if n_samples < 4 {
        bail!("need at least 4 samples for ICA, got {}", n_samples);
    }
    let mut kept_assays: Vec<String> = Vec::new();
    let mut dropped_assays = 0usize;
    for assay in assay_meta.keys() {
        let present = sample_order
            .iter()
            .filter(|s| abundance_by_key.contains_key(&(assay.clone(), (*s).clone())))
            .count();
        let missing_fraction = 1.0 - present as f64 / n_samples as f64;
        if missing_fraction <= args.max_missing_fraction + 1e-12 {
            kept_assays.push(assay.clone());
        } else {
            dropped_assays += 1;
        }
    }
    if kept_assays.is_empty() {
        bail!(
            "no assays retained with --max-missing-fraction={} over {} samples in {:?}",
            args.max_missing_fraction,
            n_samples,
            tsv
        );
    }
    eprintln!(
        "decompose ica: retained {} assays, dropped {} for exceeding --max-missing-fraction={}",
        kept_assays.len(),
        dropped_assays,
        args.max_missing_fraction
    );

    let assays: Vec<AssayInfo> = kept_assays
        .iter()
        .map(|a| AssayInfo {
            assay_id: a.clone(),
            gene_symbol: assay_meta.get(a).cloned().flatten().unwrap_or_default(),
        })
        .collect();
    let samples: Vec<SampleInfo> = sample_order
        .iter()
        .map(|s| SampleInfo {
            sample_id: s.clone(),
            subject_id: subject_lookup.get(s).cloned().unwrap_or_default(),
        })
        .collect();

    let mut data = vec![vec![0.0_f64; assays.len()]; samples.len()];
    let mut missing_cells = 0usize;
    for (j, assay) in assays.iter().enumerate() {
        let mut values: Vec<Option<f64>> = Vec::with_capacity(samples.len());
        let mut sum = 0.0_f64;
        let mut count = 0usize;
        for sample in &samples {
            match abundance_by_key.get(&(assay.assay_id.clone(), sample.sample_id.clone())) {
                Some((v, _)) => {
                    values.push(Some(*v));
                    sum += *v;
                    count += 1;
                }
                None => values.push(None),
            }
        }
        let mean = if count > 0 { sum / count as f64 } else { 0.0 };
        for (i, v) in values.iter().enumerate() {
            match v {
                Some(x) => data[i][j] = *x,
                None => {
                    if impute_mean {
                        data[i][j] = mean;
                        missing_cells += 1;
                    } else {
                        bail!(
                            "assay {} is missing sample {} and --impute=none; rerun with --impute mean or lower --max-missing-fraction",
                            assay.assay_id,
                            samples[i].sample_id
                        );
                    }
                }
            }
        }
    }
    if missing_cells > 0 {
        eprintln!(
            "decompose ica: mean-imputed {} missing cells across {} assays",
            missing_cells,
            assays.len()
        );
    }
    Ok(AbundanceMatrix {
        samples,
        assays,
        data,
    })
}

/// Canonicalize an ICA result: sign-flip each program so its max |loading| is positive,
/// then sort programs by descending max |loading|. This makes output deterministic
/// under sign ambiguity and yields stable program ordering for reporting.
/// Metadata written into `transform_applied.json` alongside
/// `loadings.tsv` whenever a non-`none` transform runs. Lets reviewers
/// audit which coordinate system the archetypes live in.
struct TransformMeta {
    name: &'static str,
    alr_reference: Option<String>,
    n_input_proteins: usize,
    n_output_coords: usize,
}

fn apply_compositional(
    matrix: AbundanceMatrix,
    args: &IcaArgs,
) -> Result<(AbundanceMatrix, TransformMeta)> {
    let transform_name = match args.transform {
        TransformArg::None => "none",
        TransformArg::Clr => "clr",
        TransformArg::Alr => "alr",
        TransformArg::Ilr => "ilr",
        TransformArg::RatioAnchor => "ratio-anchor",
    };
    if args.transform == TransformArg::None {
        let n = matrix.assays.len();
        return Ok((
            matrix,
            TransformMeta {
                name: transform_name,
                alr_reference: None,
                n_input_proteins: n,
                n_output_coords: n,
            },
        ));
    }
    if matches!(args.transform, TransformArg::Ilr) && args.alr_reference.is_some() {
        eprintln!("decompose ica: --alr-reference is ignored for --transform ilr");
    }
    let needs_reference = matches!(
        args.transform,
        TransformArg::Alr | TransformArg::RatioAnchor
    );
    let resolved_reference: Option<(usize, String)> = if needs_reference {
        let gene = match args.alr_reference.as_deref() {
            Some(g) if !g.is_empty() => g,
            _ => bail!("--transform {transform_name} requires --alr-reference <gene_symbol>"),
        };
        let idx = matrix
            .assays
            .iter()
            .position(|a| a.gene_symbol == gene)
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "--alr-reference {gene:?} not found among {} assays",
                    matrix.assays.len()
                )
            })?;
        Some((idx, gene.to_string()))
    } else {
        None
    };
    let core_transform = match args.transform {
        TransformArg::None => Transform::None,
        TransformArg::Clr => Transform::Clr,
        TransformArg::Ilr => Transform::Ilr,
        TransformArg::Alr => Transform::Alr {
            reference_index: resolved_reference.as_ref().unwrap().0,
        },
        TransformArg::RatioAnchor => Transform::RatioAnchor {
            reference_index: resolved_reference.as_ref().unwrap().0,
        },
    };
    let n_input = matrix.assays.len();
    let transformed = apply_transform(&matrix.data, core_transform)
        .map_err(|e| anyhow::anyhow!("transform {transform_name} failed: {e}"))?;
    let n_output = if transformed.is_empty() {
        0
    } else {
        transformed[0].len()
    };
    let assays = if matches!(args.transform, TransformArg::Ilr) {
        // ILR coordinates are in the Helmert basis, not raw
        // protein space. Emit synthetic assay IDs so downstream
        // writers don't try to match them back to a protein.
        (0..n_output)
            .map(|i| AssayInfo {
                assay_id: format!("ilr_coord_{:03}", i + 1),
                gene_symbol: format!("ILR{:03}", i + 1),
            })
            .collect()
    } else {
        matrix.assays
    };
    Ok((
        AbundanceMatrix {
            samples: matrix.samples,
            assays,
            data: transformed,
        },
        TransformMeta {
            name: transform_name,
            alr_reference: resolved_reference.map(|(_, g)| g),
            n_input_proteins: n_input,
            n_output_coords: n_output,
        },
    ))
}

fn write_transform_applied(path: &Path, meta: &TransformMeta) -> Result<()> {
    let mut obj = serde_json::Map::new();
    obj.insert("transform".into(), json!(meta.name));
    obj.insert(
        "alr_reference".into(),
        match &meta.alr_reference {
            Some(g) => json!(g),
            None => serde_json::Value::Null,
        },
    );
    obj.insert("n_input_proteins".into(), json!(meta.n_input_proteins));
    obj.insert("n_output_coords".into(), json!(meta.n_output_coords));
    atomic_write(
        path,
        serde_json::to_string_pretty(&serde_json::Value::Object(obj))?.as_bytes(),
    )?;
    Ok(())
}

fn canonicalize(result: &IcaResult) -> CanonicalIca {
    let k = result.mixing[0].len();
    let p = result.mixing.len();

    let mut signs = vec![1.0_f64; k];
    let mut loadings: Vec<Vec<f64>> = vec![vec![0.0; p]; k];
    for (c, sign) in signs.iter_mut().enumerate() {
        let mut max_abs = 0.0_f64;
        let mut argmax = 0usize;
        for (j, row) in result.mixing.iter().enumerate() {
            if row[c].abs() > max_abs {
                max_abs = row[c].abs();
                argmax = j;
            }
        }
        if result.mixing[argmax][c] < 0.0 {
            *sign = -1.0;
        }
        for (j, slot) in loadings[c].iter_mut().enumerate() {
            *slot = *sign * result.mixing[j][c];
        }
    }
    let activations: Vec<Vec<f64>> = signs
        .iter()
        .enumerate()
        .map(|(c, &sign)| {
            result
                .sources
                .iter()
                .map(|src_row| sign * src_row[c])
                .collect()
        })
        .collect();
    let mut order: Vec<usize> = (0..k).collect();
    order.sort_by(|&a, &b| {
        let ma = loadings[a].iter().fold(0.0_f64, |acc, v| acc.max(v.abs()));
        let mb = loadings[b].iter().fold(0.0_f64, |acc, v| acc.max(v.abs()));
        mb.partial_cmp(&ma).unwrap_or(Ordering::Equal)
    });
    let ordered_loadings: Vec<Vec<f64>> = order.iter().map(|&i| loadings[i].clone()).collect();
    let ordered_activations: Vec<Vec<f64>> =
        order.iter().map(|&i| activations[i].clone()).collect();
    CanonicalIca {
        loadings: ordered_loadings,
        activations: ordered_activations,
    }
}

struct CanonicalIca {
    /// `k` rows, each of length `p`.
    loadings: Vec<Vec<f64>>,
    /// `k` rows, each of length `n`.
    activations: Vec<Vec<f64>>,
}

fn program_name(index: usize) -> String {
    format!("program_{:02}", index + 1)
}

fn write_loadings(path: &Path, matrix: &AbundanceMatrix, run: &CanonicalIca) -> Result<()> {
    let mut out = String::from("program\tassay_id\tgene_symbol\tloading\n");
    for (c, row) in run.loadings.iter().enumerate() {
        let program = program_name(c);
        for (j, load) in row.iter().enumerate() {
            let assay = &matrix.assays[j];
            out.push_str(&program);
            out.push('\t');
            out.push_str(&assay.assay_id);
            out.push('\t');
            out.push_str(&assay.gene_symbol);
            out.push('\t');
            out.push_str(&format_float(*load));
            out.push('\n');
        }
    }
    atomic_write(path, out.as_bytes())
}

fn write_activations(
    path: &Path,
    matrix: &AbundanceMatrix,
    run: &CanonicalIca,
    cohort: &str,
) -> Result<()> {
    let mut out = String::from("cohort\tsubject_id\tsample_id\tprogram\tactivation\n");
    for (c, row) in run.activations.iter().enumerate() {
        let program = program_name(c);
        for (i, v) in row.iter().enumerate() {
            let sample = &matrix.samples[i];
            out.push_str(cohort);
            out.push('\t');
            out.push_str(&sample.subject_id);
            out.push('\t');
            out.push_str(&sample.sample_id);
            out.push('\t');
            out.push_str(&program);
            out.push('\t');
            out.push_str(&format_float(*v));
            out.push('\n');
        }
    }
    atomic_write(path, out.as_bytes())
}

#[allow(clippy::too_many_arguments)]
fn write_stability(
    path: &Path,
    reference: &CanonicalIca,
    others: &[(u64, IcaResult)],
    top_n: usize,
    jaccard_threshold: f64,
    min_stable_seed_fraction: f64,
    ref_seed: u64,
) -> Result<()> {
    let mut out = String::from(
        "program\treference_seed\tn_alt_seeds\tn_stable_runs\tseed_stability_fraction\tmean_best_jaccard\tflag_below_threshold\n",
    );
    let canon_alts: Vec<CanonicalIca> = others.iter().map(|(_, r)| canonicalize(r)).collect();
    for (c, row) in reference.loadings.iter().enumerate() {
        let program = program_name(c);
        let mut n_stable = 0usize;
        let mut jaccards: Vec<f64> = Vec::with_capacity(canon_alts.len());
        for alt in &canon_alts {
            let best = alt
                .loadings
                .iter()
                .map(|alt_row| jaccard_top_n(row, alt_row, top_n))
                .fold(0.0_f64, f64::max);
            jaccards.push(best);
            if best >= jaccard_threshold {
                n_stable += 1;
            }
        }
        let n_alt = canon_alts.len();
        let fraction = if n_alt == 0 {
            1.0
        } else {
            n_stable as f64 / n_alt as f64
        };
        let mean = if jaccards.is_empty() {
            f64::NAN
        } else {
            jaccards.iter().sum::<f64>() / jaccards.len() as f64
        };
        let flag = if n_alt > 0 && fraction < min_stable_seed_fraction {
            "1"
        } else {
            "0"
        };
        out.push_str(&format!(
            "{program}\t{ref_seed}\t{n_alt}\t{n_stable}\t{}\t{}\t{flag}\n",
            format_float(fraction),
            format_float(mean),
        ));
    }
    atomic_write(path, out.as_bytes())
}

#[derive(ClapArgs, Debug)]
pub struct UnmixArgs {
    /// Canonical Atman input directory (expects `measurements.tsv`,
    /// `samples.tsv`, `proteins.tsv`).
    #[arg(long)]
    input_dir: PathBuf,

    /// Number of endmembers. Pass an integer for a fixed `k`, or
    /// `auto` to sweep `[--k-min, --k-max]` and pick the smallest
    /// `k` whose marginal reconstruction-residual improvement drops
    /// below `--k-elbow-threshold` times the sweep's peak
    /// improvement. Must satisfy `2 ≤ k ≤ n_samples / 2`.
    #[arg(long)]
    k: String,

    /// Lower bound for `--k auto` sweep.
    #[arg(long, default_value_t = 2)]
    k_min: usize,

    /// Upper bound for `--k auto` sweep.
    #[arg(long, default_value_t = 8)]
    k_max: usize,

    /// Elbow threshold for `--k auto`: marginal improvement / peak
    /// improvement ratio at which the sweep stops.
    #[arg(long, default_value_t = 0.10)]
    k_elbow_threshold: f64,

    /// Endmember extraction method. `vca` (default): deterministic
    /// projection-onto-complement. `nfindr`: iterative simplex-volume
    /// maximization initialized from VCA.
    #[arg(long, default_value = "vca")]
    method: String,

    /// Max N-FINDR passes (ignored for `--method vca`). A full pass
    /// tries every (endmember slot × candidate sample) swap.
    #[arg(long, default_value_t = 20)]
    nfindr_max_passes: usize,

    /// Abundance estimator. `fcls` (default) enforces simplex
    /// constraints `α ≥ 0 ∧ Σα = 1`. `ucls` drops them.
    #[arg(long, default_value = "fcls")]
    abundance: String,

    /// Pre-transform applied to the subject × protein matrix before
    /// unmixing. `none` passes atman's already-log-scale canonical
    /// data through; `clr` centers each subject row on its own mean
    /// (valid only with `--abundance ucls`, or with
    /// `--allow-unconstrained-simplex` for FCLS); `ilr` is refused
    /// because it changes the feature ordering.
    #[arg(long, default_value = "none")]
    transform: String,

    /// Reference gene symbol for `--transform alr` / `ratio-anchor`.
    #[arg(long)]
    alr_reference: Option<String>,

    /// Escape hatch: run FCLS even with a compositional transform
    /// whose simplex interpretation is debatable (CLR, ILR). Off by
    /// default — FCLS's sum-to-one constraint has no natural
    /// meaning on log-ratio coordinates.
    #[arg(long, default_value_t = false)]
    allow_unconstrained_simplex: bool,

    /// Seed for VCA's initial projection direction. Sub-seeds for
    /// each VCA step derive from this via SplitMix64.
    #[arg(long, default_value_t = 20260420)]
    seed: u64,

    /// FCLS maximum projected-gradient iterations.
    #[arg(long, default_value_t = 1000)]
    fcls_max_iter: usize,

    /// FCLS convergence tolerance (`max |Δα|` between iterations).
    #[arg(long, default_value_t = 1e-9)]
    fcls_tol: f64,

    /// Subject-level bootstrap iterations for loading + abundance CI.
    /// 0 (default) disables bootstrap — endmembers.tsv and
    /// abundances.tsv carry `NA` in the CI columns.
    #[arg(long, default_value_t = 0)]
    n_boot: usize,

    /// Path to a marker-set TSV with columns (`set_name`,
    /// `gene_symbol`) for hypergeometric ORA on each endmember's
    /// top-N loadings. Emits `endmember_annotations.tsv` alongside
    /// the other outputs.
    #[arg(long)]
    annotate_markers: Option<PathBuf>,

    /// Top-N threshold used when ranking an endmember's proteins for
    /// ORA against `--annotate-markers`.
    #[arg(long, default_value_t = 20)]
    annotate_top_n: usize,

    /// Drop proteins with more than this fraction of missing samples.
    #[arg(long, default_value_t = 0.0)]
    max_missing_fraction: f64,

    /// `--impute none|mean` for residual missingness.
    #[arg(long, default_value = "none")]
    impute: String,

    /// Output directory for `endmembers.tsv`, `abundances.tsv`,
    /// `unmix_diagnostics.tsv`.
    #[arg(long)]
    output_dir: PathBuf,
}

fn run_unmix(args: UnmixArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let endmember_method = match args.method.as_str() {
        "vca" => EndmemberMethod::Vca,
        "nfindr" => EndmemberMethod::Nfindr {
            max_passes: args.nfindr_max_passes,
        },
        other => bail!("--method {other:?}; expected `vca` or `nfindr`"),
    };
    let abundance_method = match args.abundance.as_str() {
        "fcls" => AbundanceMethod::Fcls,
        "ucls" => AbundanceMethod::Ucls,
        other => bail!("--abundance {other:?}; expected fcls or ucls"),
    };
    match args.transform.as_str() {
        "none" | "log" => {}
        "clr" | "alr" | "ratio-anchor" => {
            if matches!(abundance_method, AbundanceMethod::Fcls)
                && !args.allow_unconstrained_simplex
            {
                bail!(
                    "--transform {:?} combined with --abundance fcls is refused: \
                     log-ratio coordinates do not admit a `Σα = 1` interpretation. \
                     Pass --allow-unconstrained-simplex to override, or switch to \
                     --abundance ucls.",
                    args.transform
                );
            }
        }
        "ilr" => {
            bail!(
                "--transform ilr is refused in `decompose unmix`: ILR changes \
                 the feature ordering, so recovered endmember loadings would \
                 not line up with the input protein labels"
            );
        }
        other => bail!("--transform {other:?}; expected none|log|clr|alr|ratio-anchor"),
    }
    let impute_mean = match args.impute.as_str() {
        "none" => false,
        "mean" => true,
        other => bail!("--impute {other:?}; expected none or mean"),
    };
    if args.fcls_max_iter < 1 {
        bail!("--fcls-max-iter must be >= 1");
    }

    // Load subject × protein matrix, keyed by gene_symbol.
    let SubjectProteinMatrix {
        sample_ids: subject_ids,
        protein_labels,
        data,
    } = load_subject_protein_matrix(
        &args.input_dir.join("measurements.tsv"),
        args.max_missing_fraction,
        impute_mean,
    )?;
    // Resolve --k: either an integer or "auto" (sweep + elbow).
    let (effective_k, k_sweep): (usize, Vec<KSweepRow>) = if args.k == "auto" {
        if args.k_min < 2 || args.k_max < args.k_min {
            bail!(
                "invalid --k auto sweep: k-min={} k-max={}",
                args.k_min,
                args.k_max
            );
        }
        if args.k_max * 2 > subject_ids.len() {
            bail!(
                "--k-max={} requires >= {} samples; got {}",
                args.k_max,
                args.k_max * 2,
                subject_ids.len()
            );
        }
        eprintln!(
            "decompose unmix: --k auto sweeping [{}, {}]",
            args.k_min, args.k_max
        );
        // Apply transform later; here we need it already for the
        // sweep's VCA+FCLS runs. Rebuild the transformed matrix
        // early.
        // (The transform below is still the canonical one — we
        // just need it before the sweep runs.)
        // Sweep uses VCA; caller can re-run with nfindr post-hoc.
        let sweep_cfg = UnmixConfig {
            seed: args.seed,
            endmember_method: EndmemberMethod::Vca,
            abundance_method,
            fcls_max_iter: args.fcls_max_iter,
            fcls_tol: args.fcls_tol,
        };
        select_k_auto(
            &data,
            args.k_min,
            args.k_max,
            sweep_cfg,
            args.k_elbow_threshold,
        )
        .map_err(|e| anyhow::anyhow!(e))?
    } else {
        let k: usize = args
            .k
            .parse()
            .with_context(|| format!("--k {:?}: expected an integer or `auto`", args.k))?;
        (k, Vec::new())
    };
    if subject_ids.len() < effective_k * 2 {
        bail!(
            "k={} requires >= {} samples (k > n/2 is under-determined); got {}",
            effective_k,
            effective_k * 2,
            subject_ids.len()
        );
    }
    // Apply transform (if any) — "log" is a no-op on atman canonical
    // data since it's already on a log scale.
    let transform = match args.transform.as_str() {
        "none" | "log" => Transform::None,
        "clr" => Transform::Clr,
        "alr" => {
            let reference = args
                .alr_reference
                .as_deref()
                .context("--transform alr requires --alr-reference <gene>")?;
            let idx = protein_labels
                .iter()
                .position(|l| l == reference)
                .with_context(|| format!("--alr-reference {reference:?} not in protein labels"))?;
            Transform::Alr {
                reference_index: idx,
            }
        }
        "ratio-anchor" => {
            let reference = args
                .alr_reference
                .as_deref()
                .context("--transform ratio-anchor requires --alr-reference <gene>")?;
            let idx = protein_labels
                .iter()
                .position(|l| l == reference)
                .with_context(|| format!("--alr-reference {reference:?} not in protein labels"))?;
            Transform::RatioAnchor {
                reference_index: idx,
            }
        }
        _ => unreachable!(),
    };
    let transformed =
        apply_transform(&data, transform).map_err(|e| anyhow::anyhow!("transform failed: {e}"))?;

    let cfg = UnmixConfig {
        seed: args.seed,
        endmember_method,
        abundance_method,
        fcls_max_iter: args.fcls_max_iter,
        fcls_tol: args.fcls_tol,
    };
    let result = unmix(&transformed, effective_k, cfg).map_err(|e| anyhow::anyhow!(e))?;
    let (loading_ci, abundance_ci) =
        bootstrap_ci(&result, &transformed, cfg, args.n_boot).map_err(|e| anyhow::anyhow!(e))?;

    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating {:?}", args.output_dir))?;
    let endmembers_path = args.output_dir.join("endmembers.tsv");
    write_unmix_endmembers(
        &endmembers_path,
        &result,
        &protein_labels,
        &subject_ids,
        loading_ci.as_ref(),
    )?;
    let abundances_path = args.output_dir.join("abundances.tsv");
    write_unmix_abundances(
        &abundances_path,
        &result,
        &subject_ids,
        abundance_ci.as_ref(),
    )?;
    let diag_path = args.output_dir.join("unmix_diagnostics.tsv");
    write_unmix_diagnostics(&diag_path, &result, &subject_ids)?;
    let mut outputs = vec![
        endmembers_path.clone(),
        abundances_path.clone(),
        diag_path.clone(),
    ];
    if !k_sweep.is_empty() {
        let k_path = args.output_dir.join("k_selection.tsv");
        write_k_selection(&k_path, &k_sweep, effective_k)?;
        outputs.push(k_path);
    }
    if let Some(markers_path) = &args.annotate_markers {
        let marker_sets = read_marker_sets(markers_path)?;
        let ann_rows = ora_enrichment(
            &result.endmember_loadings,
            &protein_labels,
            &marker_sets,
            args.annotate_top_n,
        );
        let ann_path = args.output_dir.join("endmember_annotations.tsv");
        write_endmember_annotations(&ann_path, &ann_rows)?;
        outputs.push(ann_path);
    }

    eprintln!(
        "decompose unmix: k={} n={} p={} method={} abundance={} transform={}",
        effective_k,
        subject_ids.len(),
        protein_labels.len(),
        args.method,
        args.abundance,
        args.transform,
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &[
            "measurements.tsv",
            "measurements.tsv",
            "samples.tsv",
            "proteins.tsv",
        ],
    )?;
    let sidecar = sidecar_path_for(&endmembers_path);
    write_run_sidecar(
        &sidecar,
        "decompose unmix",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "k": args.k,
            "effective-k": effective_k,
            "k-min": args.k_min,
            "k-max": args.k_max,
            "k-elbow-threshold": args.k_elbow_threshold,
            "method": args.method,
            "nfindr-max-passes": args.nfindr_max_passes,
            "abundance": args.abundance,
            "transform": args.transform,
            "alr-reference": args.alr_reference,
            "allow-unconstrained-simplex": args.allow_unconstrained_simplex,
            "seed": args.seed,
            "fcls-max-iter": args.fcls_max_iter,
            "fcls-tol": args.fcls_tol,
            "n-boot": args.n_boot,
            "annotate-markers": args.annotate_markers.as_ref().map(|p| p.display().to_string()),
            "annotate-top-n": args.annotate_top_n,
            "max-missing-fraction": args.max_missing_fraction,
            "impute": args.impute,
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        {
            // Resolved-value fields: when `--k auto` swept for the
            // elbow, `k_resolved` carries the integer actually used
            // and `k_resolution_source` names the rule. Otherwise
            // `k_resolution_source` is `cli`.
            let mut extras = serde_json::Map::new();
            extras.insert("k_resolved".into(), serde_json::json!(effective_k));
            let source = if args.k == "auto" {
                format!(
                    "auto:elbow(k_min={},k_max={},threshold={})",
                    args.k_min, args.k_max, args.k_elbow_threshold
                )
            } else {
                "cli".to_string()
            };
            extras.insert("k_resolution_source".into(), serde_json::json!(source));
            Some(extras)
        },
    )?;
    eprintln!("decompose unmix: sidecar={}", sidecar.display());
    Ok(())
}

fn read_marker_sets(path: &Path) -> Result<BTreeMap<String, Vec<String>>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let name_idx = headers
        .iter()
        .position(|h| h == "set_name")
        .context("marker TSV missing set_name column")?;
    let gene_idx = headers
        .iter()
        .position(|h| h == "gene_symbol")
        .context("marker TSV missing gene_symbol column")?;
    let mut out: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let name = row[name_idx].to_string();
        let gene = row[gene_idx].trim().to_string();
        if !name.is_empty() && !gene.is_empty() {
            out.entry(name).or_default().push(gene);
        }
    }
    if out.is_empty() {
        bail!("marker TSV {:?} produced zero sets", path);
    }
    Ok(out)
}

fn write_endmember_annotations(path: &Path, rows: &[AnnotationRow]) -> Result<()> {
    let mut buf = String::from(
        "endmember_id\tset_name\tuniverse_size\tset_size_in_universe\t\
         top_n\toverlap\tp_value\n",
    );
    for r in rows {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.6e}\n",
            r.endmember_id,
            r.set_name,
            r.universe_size,
            r.set_size_in_universe,
            r.top_n,
            r.overlap,
            r.p_value,
        ));
    }
    atomic_write(path, buf.as_bytes())
}

fn write_k_selection(path: &Path, rows: &[KSweepRow], chosen_k: usize) -> Result<()> {
    let mut buf = String::from("k\tmean_residual_norm\tmarginal_improvement\tchosen\n");
    for r in rows {
        buf.push_str(&format!(
            "{}\t{:.6}\t{:.6}\t{}\n",
            r.k,
            r.mean_residual_norm,
            r.marginal_improvement,
            if r.k == chosen_k { 1 } else { 0 },
        ));
    }
    atomic_write(path, buf.as_bytes())
}

struct SubjectProteinMatrix {
    /// Sample IDs in row order.
    sample_ids: Vec<String>,
    /// Retained protein labels in column order.
    protein_labels: Vec<String>,
    /// `[sample][protein]` abundance after missingness filter / imputation.
    data: Vec<Vec<f64>>,
}

fn load_subject_protein_matrix(
    path: &Path,
    max_missing_fraction: f64,
    impute_mean: bool,
) -> Result<SubjectProteinMatrix> {
    let records = read_measurements_long(path)?;
    let mut samples: BTreeSet<String> = BTreeSet::new();
    let mut genes: BTreeSet<String> = BTreeSet::new();
    let mut cells: BTreeMap<(String, String), f64> = BTreeMap::new();
    for r in &records {
        if r.dropped_by_qc {
            continue;
        }
        let v = r.abundance.as_f64();
        if !v.is_finite() {
            continue;
        }
        let gene = match &r.gene_symbol {
            Some(g) => g.clone(),
            None => continue,
        };
        samples.insert(r.sample_id.clone());
        genes.insert(gene.clone());
        cells.insert((gene, r.sample_id.clone()), v);
    }
    let sample_list: Vec<String> = samples.into_iter().collect();
    let gene_list: Vec<String> = genes.into_iter().collect();
    if sample_list.is_empty() || gene_list.is_empty() {
        bail!("no usable measurements in {:?}", path);
    }
    // Build with missingness map.
    let mut raw = vec![vec![f64::NAN; gene_list.len()]; sample_list.len()];
    for (si, sid) in sample_list.iter().enumerate() {
        for (gi, gene) in gene_list.iter().enumerate() {
            if let Some(&v) = cells.get(&(gene.clone(), sid.clone())) {
                raw[si][gi] = v;
            }
        }
    }
    // Drop features exceeding the missing-fraction budget.
    let n = sample_list.len() as f64;
    let keep: Vec<bool> = (0..gene_list.len())
        .map(|gi| {
            let present = (0..sample_list.len())
                .filter(|&si| raw[si][gi].is_finite())
                .count();
            let missing = 1.0 - (present as f64 / n);
            missing <= max_missing_fraction + 1e-12
        })
        .collect();
    let kept_genes: Vec<String> = gene_list
        .iter()
        .zip(keep.iter())
        .filter_map(|(g, k)| if *k { Some(g.clone()) } else { None })
        .collect();
    if kept_genes.is_empty() {
        bail!(
            "no proteins retained after --max-missing-fraction {} filter",
            max_missing_fraction
        );
    }
    let mut data = vec![vec![0.0_f64; kept_genes.len()]; sample_list.len()];
    for (new_gi, gi) in (0..gene_list.len()).filter(|i| keep[*i]).enumerate() {
        // Compute column mean for imputation.
        let column: Vec<f64> = raw.iter().map(|row| row[gi]).collect();
        let (sum, count) = column
            .iter()
            .filter(|v| v.is_finite())
            .fold((0.0_f64, 0_usize), |(s, n), v| (s + v, n + 1));
        let mean = if count > 0 { sum / count as f64 } else { 0.0 };
        for (si, data_row) in data.iter_mut().enumerate() {
            let v = column[si];
            if v.is_finite() {
                data_row[new_gi] = v;
            } else if impute_mean {
                data_row[new_gi] = mean;
            } else {
                bail!(
                    "missing value at sample {} protein {}; rerun with --impute mean",
                    sample_list[si],
                    kept_genes[new_gi]
                );
            }
        }
    }
    Ok(SubjectProteinMatrix {
        sample_ids: sample_list,
        protein_labels: kept_genes,
        data,
    })
}

fn write_unmix_endmembers(
    path: &Path,
    result: &UnmixResult,
    protein_labels: &[String],
    subject_ids: &[String],
    loading_ci: Option<&LoadingCi>,
) -> Result<()> {
    let mut buf = String::from(
        "endmember_id\tsource_sample_id\tprotein\tloading\trank_in_endmember\t\
         loading_ci_lower\tloading_ci_upper\n",
    );
    for (ai, loading) in result.endmember_loadings.iter().enumerate() {
        let src_idx = result.endmember_sample_indices[ai];
        let source_sid = subject_ids
            .get(src_idx)
            .cloned()
            .unwrap_or_else(|| format!("sample{src_idx}"));
        let mut ranking: Vec<(usize, f64)> = loading
            .iter()
            .enumerate()
            .map(|(i, v)| (i, v.abs()))
            .collect();
        ranking.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
        let mut rank_for = vec![0usize; loading.len()];
        for (r, (i, _)) in ranking.iter().enumerate() {
            rank_for[*i] = r + 1;
        }
        for (pi, &v) in loading.iter().enumerate() {
            let (lo, hi) = loading_ci
                .map(|ci| (ci.lower[ai][pi], ci.upper[ai][pi]))
                .map(|(a, b)| (format!("{a:.6}"), format!("{b:.6}")))
                .unwrap_or_else(|| ("NA".into(), "NA".into()));
            buf.push_str(&format!(
                "E{:03}\t{}\t{}\t{:.6}\t{}\t{}\t{}\n",
                ai + 1,
                source_sid,
                protein_labels[pi],
                v,
                rank_for[pi],
                lo,
                hi,
            ));
        }
    }
    atomic_write(path, buf.as_bytes())
}

fn write_unmix_abundances(
    path: &Path,
    result: &UnmixResult,
    subject_ids: &[String],
    abundance_ci: Option<&AbundanceCi>,
) -> Result<()> {
    // Emit per-endmember `E001, E001_ci_lower, E001_ci_upper, E002, ...`
    // so downstream code can extract either the point estimate or the
    // interval by column name.
    let mut buf = String::from("sample_id");
    for ai in 0..result.endmember_loadings.len() {
        let lbl = format!("E{:03}", ai + 1);
        buf.push('\t');
        buf.push_str(&lbl);
        buf.push('\t');
        buf.push_str(&format!("{lbl}_ci_lower"));
        buf.push('\t');
        buf.push_str(&format!("{lbl}_ci_upper"));
    }
    buf.push('\n');
    for (si, sid) in subject_ids.iter().enumerate() {
        buf.push_str(sid);
        for (ai, &v) in result.abundances[si].iter().enumerate() {
            buf.push('\t');
            buf.push_str(&format!("{v:.6}"));
            let (lo, hi) = abundance_ci
                .map(|ci| (ci.lower[si][ai], ci.upper[si][ai]))
                .map(|(a, b)| (format!("{a:.6}"), format!("{b:.6}")))
                .unwrap_or_else(|| ("NA".into(), "NA".into()));
            buf.push('\t');
            buf.push_str(&lo);
            buf.push('\t');
            buf.push_str(&hi);
        }
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

fn write_unmix_diagnostics(
    path: &Path,
    result: &UnmixResult,
    subject_ids: &[String],
) -> Result<()> {
    let mut buf = String::from("sample_id\treconstruction_residual_norm\tabundance_sum\n");
    for (si, sid) in subject_ids.iter().enumerate() {
        let sum: f64 = result.abundances[si].iter().sum();
        buf.push_str(&format!(
            "{}\t{:.6}\t{:.6}\n",
            sid, result.residual_norms[si], sum,
        ));
    }
    atomic_write(path, buf.as_bytes())
}

// ── counterfactual ──────────────────────────────────────────────────

fn run_counterfactual(args: CounterfactualArgs) -> Result<()> {
    let started_at = SystemTime::now();

    // Read loadings: program, assay_id, gene_symbol, loading
    let mut loadings_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&args.loadings)
        .with_context(|| format!("opening {:?}", args.loadings))?;
    let mut loadings_map: BTreeMap<(String, String), BTreeMap<String, f64>> = BTreeMap::new();
    let mut all_programs: BTreeSet<String> = BTreeSet::new();
    for record in loadings_reader.records() {
        let r = record?;
        let program = r.get(0).unwrap_or("").to_string();
        let assay_id = r.get(1).unwrap_or("").to_string();
        let gene_symbol = r.get(2).unwrap_or("").to_string();
        let loading: f64 = r
            .get(3)
            .unwrap_or("")
            .parse()
            .with_context(|| format!("parsing loading for {}/{}", program, assay_id))?;
        all_programs.insert(program.clone());
        loadings_map
            .entry((assay_id, gene_symbol))
            .or_default()
            .insert(program, loading);
    }

    if !all_programs.contains(&args.set_program) {
        bail!(
            "--set-program {:?} not found in loadings; available: {}",
            args.set_program,
            all_programs.iter().cloned().collect::<Vec<_>>().join(", ")
        );
    }

    // Read activations: cohort, subject_id, sample_id, program, activation
    let mut act_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&args.activations)
        .with_context(|| format!("opening {:?}", args.activations))?;
    let mut activations_map: BTreeMap<(String, String), f64> = BTreeMap::new();
    let mut sample_ids: BTreeSet<String> = BTreeSet::new();
    for record in act_reader.records() {
        let r = record?;
        let sample_id = r.get(2).unwrap_or("").to_string();
        let program = r.get(3).unwrap_or("").to_string();
        let activation: f64 = r
            .get(4)
            .unwrap_or("")
            .parse()
            .with_context(|| format!("parsing activation for {}/{}", sample_id, program))?;
        sample_ids.insert(sample_id.clone());
        activations_map.insert((sample_id, program), activation);
    }
    let sample_ids: Vec<String> = sample_ids.into_iter().collect();

    // Compute delta: for each protein × sample, the contribution of the
    // target program is loading[protein, program] × (activation_original - to).
    let mut buf = String::from("assay_id\tgene_symbol");
    for sid in &sample_ids {
        buf.push('\t');
        buf.push_str(sid);
    }
    buf.push('\n');

    let mut n_proteins = 0u64;
    for ((assay_id, gene_symbol), per_program) in &loadings_map {
        let loading = per_program.get(&args.set_program).copied().unwrap_or(0.0);
        buf.push_str(assay_id);
        buf.push('\t');
        buf.push_str(gene_symbol);
        for sid in &sample_ids {
            let original_activation = activations_map
                .get(&(sid.clone(), args.set_program.clone()))
                .copied()
                .unwrap_or(0.0);
            let delta = loading * (original_activation - args.to);
            buf.push('\t');
            buf.push_str(&format_float(delta));
        }
        buf.push('\n');
        n_proteins += 1;
    }

    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent)?;
    }
    atomic_write(&args.output, buf.as_bytes())?;

    eprintln!(
        "decompose counterfactual: {} proteins × {} samples, program={} set to {}, output={}",
        n_proteins,
        sample_ids.len(),
        args.set_program,
        args.to,
        args.output.display()
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = crate::io::hash_labeled_inputs(&[
        ("loadings", args.loadings.as_path()),
        ("activations", args.activations.as_path()),
    ])?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "decompose counterfactual",
        json!({
            "loadings": args.loadings.display().to_string(),
            "activations": args.activations.display().to_string(),
            "set-program": args.set_program,
            "to": args.to,
            "output": args.output.display().to_string(),
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("decompose counterfactual: sidecar={}", sidecar.display());
    Ok(())
}
