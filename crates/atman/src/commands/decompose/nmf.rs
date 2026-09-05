use anyhow::{bail, Context, Result};
use atman_core::nmf::{
    multi_seed_nmf, nmf as nmf_core, select_k as nmf_select_k, BetaLoss, Init, KSelection,
    MultiSeedConfig, NmfConfig, StabilityMetric as NmfStabilityMetric,
};
use clap::Args as ClapArgs;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;
use std::time::SystemTime;

use super::program_name;
use crate::io::{
    atomic_write, format_float, hash_canonical_inputs, read_measurements_long, sidecar_path_for,
    write_run_sidecar,
};

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

    /// Drop assays where the fraction of missing samples exceeds this value. 0 = strict complete-case.
    #[arg(long, default_value_t = 0.0)]
    pub max_missing_fraction: f64,
}

pub(super) fn nmf_run(args: NmfArgs) -> Result<()> {
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
        other => bail!("--init {:?}: expected `nndsvda` or `random`", other),
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

    if !(0.0..=1.0).contains(&args.max_missing_fraction) {
        bail!("--max-missing-fraction must be in [0, 1]");
    }

    // Drop assays whose missing-sample fraction exceeds the configured
    // threshold (same rule as `decompose ica`). The retained matrix must
    // still be complete; any residual holes are caught by the dense-matrix
    // build below.
    let mut kept_assays: Vec<String> = Vec::new();
    let mut n_assays_dropped_missingness = 0usize;
    for assay in assay_meta.keys() {
        let present = sample_order
            .iter()
            .filter(|s| abundance_by_key.contains_key(&(assay.clone(), (*s).clone())))
            .count();
        let missing_fraction = 1.0 - present as f64 / n_samples as f64;
        if missing_fraction <= args.max_missing_fraction + 1e-12 {
            kept_assays.push(assay.clone());
        } else {
            n_assays_dropped_missingness += 1;
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
        "decompose nmf: retained {} assays, dropped {} for exceeding --max-missing-fraction={}",
        kept_assays.len(),
        n_assays_dropped_missingness,
        args.max_missing_fraction
    );

    let assay_order: Vec<String> = kept_assays;
    let n_assays = assay_order.len();
    let n_assays_retained = n_assays;

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
                    sample,
                    assay
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
        other => bail!("--stability-metric {:?}: expected `jaccard-top20`", other),
    };

    // ── build NmfConfig base (k will be finalised after k-selection) ─────────
    let nmf_cfg_base = NmfConfig {
        k: if k_sel == KSelection::Fixed {
            k_fixed
        } else {
            k_min
        }, // temp; overridden below
        beta_loss,
        init,
        max_iter: args.max_iter,
        tol: args.tol,
        seed: if args.n_seeds == 1 {
            args.seed
        } else {
            args.seed_base
        },
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
        eprintln!("decompose nmf: selected_k={}", sel_result.selected_k,);
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
        h: Vec<Vec<f64>>, // k × p
        w: Vec<Vec<f64>>, // n × k
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
                if cnt > 0.0 {
                    total / cnt
                } else {
                    1.0
                }
            };
            let scale = (mean_x / k_surv.max(1) as f64).sqrt().max(1e-9);
            let mut rng = atman_core::ica::Xoshiro256pp::new(args.seed_base.wrapping_add(9999));
            let mut w: Vec<Vec<f64>> = (0..n_samples)
                .map(|_| {
                    (0..k_surv)
                        .map(|_| rng.next_normal().abs() * scale)
                        .collect()
                })
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
                        // `j` indexes two rows (`h_fixed[a][j]` and `h_fixed[b][j]`) in lockstep.
                        #[allow(clippy::needless_range_loop)]
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
                        if w[i][a] < 0.0 {
                            w[i][a] = 0.0;
                        }
                    }
                }
            }
            w
        };

        RunOutput::MultiSeed {
            all_rows,
            surviving,
            w_refit,
        }
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
            RunOutput::MultiSeed {
                all_rows,
                surviving,
                ..
            } => {
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
    if let RunOutput::MultiSeed {
        ref all_rows,
        ref surviving,
        ..
    } = run_output
    {
        if let Some(ref stab_path) = args.output_stability {
            if let Some(parent) = stab_path.parent() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("creating stability dir {:?}", parent))?;
            }
            let mut out = String::from("program\tstable_seed_fraction\tn_seeds_present\n");
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
                    .map(format_float)
                    .unwrap_or_else(|| "NA".to_string());
                let selected = if row.k == sweep_result.selected_k {
                    "1"
                } else {
                    "0"
                };
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
            "max-missing-fraction": args.max_missing_fraction,
            "n_assays_retained": n_assays_retained,
            "n_assays_dropped_missingness": n_assays_dropped_missingness,
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
