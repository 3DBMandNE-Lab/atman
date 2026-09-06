use anyhow::{bail, Context, Result};
use atman_core::align::{build_archetypes, summarize_archetypes, AlignMetric, AlignedProgram};
use atman_core::align_bootstrap::{
    align_bootstrap, BootstrapParams, BootstrapRow, CohortMatrix, Decomposition,
};
use atman_core::align_project::{project, Atlas, ProjectionMethod};
use atman_core::compositional::{apply_transform, Transform};
use atman_core::nmf::{
    apply_transform as apply_nmf_transform, BetaLoss, Init, Transform as NmfTransform,
    TransformRecord as NmfTransformRecord, DEFAULT_EXP2_CLIP_CLAMP,
};
use clap::{Args as ClapArgs, Subcommand, ValueEnum};
use csv::ReaderBuilder;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{
    atomic_write, hash_labeled_inputs, read_measurements_long, sidecar_path_for, write_run_sidecar,
};

/// Sniff the `command` field of the neighboring `*.run.json` sidecar for a
/// loadings TSV and map it to a decomposition method label.
///
/// Returns `"nmf"`, `"ica"`, or `"unknown"` for a single loadings file.
/// When multiple files are provided, returns `"mixed"` if they disagree,
/// otherwise the common method. Emits a warning to stderr when no sidecar
/// is readable.
fn sniff_decomposition_method(loadings_paths: &[PathBuf]) -> String {
    let mut methods: Vec<String> = Vec::new();
    for path in loadings_paths {
        let sidecar = sidecar_path_for(path);
        let method = match std::fs::read_to_string(&sidecar) {
            Ok(text) => match serde_json::from_str::<serde_json::Value>(&text) {
                Ok(json) => {
                    let cmd = json["command"].as_str().unwrap_or("");
                    if cmd.contains("decompose nmf") {
                        "nmf".to_string()
                    } else if cmd.contains("decompose ica") {
                        "ica".to_string()
                    } else if cmd.is_empty() {
                        eprintln!(
                            "align: warning: sidecar {:?} has no 'command' field; \
                             recording decomposition_method=unknown",
                            sidecar
                        );
                        "unknown".to_string()
                    } else {
                        eprintln!(
                            "align: warning: sidecar {:?} has unrecognised command {:?}; \
                             recording decomposition_method=unknown",
                            sidecar, cmd
                        );
                        "unknown".to_string()
                    }
                }
                Err(_) => {
                    eprintln!(
                        "align: warning: could not parse sidecar {:?}; \
                         recording decomposition_method=unknown",
                        sidecar
                    );
                    "unknown".to_string()
                }
            },
            Err(_) => {
                eprintln!(
                    "align: warning: no sidecar found at {:?}; \
                     recording decomposition_method=unknown",
                    sidecar
                );
                "unknown".to_string()
            }
        };
        methods.push(method);
    }
    if methods.is_empty() {
        return "unknown".to_string();
    }
    let first = methods[0].clone();
    if methods.iter().all(|m| *m == first) {
        first
    } else {
        "mixed".to_string()
    }
}

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Cross-cohort program alignment and sweep mode.
    Programs(ProgramsArgs),
    /// Subject-level bootstrap of cross-cohort archetype alignment
    /// producing per-archetype universality probabilities.
    Bootstrap(BootstrapArgs),
    /// Project a new cohort's subject × protein matrix onto a trained
    /// atlas of cross-cohort archetype loading vectors, emitting
    /// per-subject activations without re-running the alignment.
    Project(ProjectArgs),
}

#[derive(ClapArgs, Debug)]
pub struct ProgramsArgs {
    /// Comma-separated per-cohort loadings TSV paths.
    #[arg(long)]
    loadings: String,

    /// Comma-separated cohort labels (same order as --loadings). Defaults to file stems.
    #[arg(long)]
    cohorts: Option<String>,

    /// Comma-separated per-cohort annotations TSV paths (same order as --loadings).
    #[arg(long)]
    annotations: Option<String>,

    /// Label column in the loadings TSV: `gene_symbol` (default), `assay_id`, or `protein`.
    #[arg(long, default_value = "gene_symbol")]
    label_col: String,

    /// Similarity metric (single-mode). Ignored when --sweep is set.
    #[arg(long, value_enum, default_value_t = SingleMetric::Jaccard)]
    metric: SingleMetric,

    /// Top-N loadings for Jaccard (single-mode).
    #[arg(long, default_value_t = 40)]
    top_n: usize,

    /// Similarity threshold (single-mode).
    #[arg(long, default_value_t = 0.15)]
    tau: f64,

    /// Require reciprocal-best matches (single-mode and sweep).
    #[arg(long, default_value_t = false)]
    reciprocal_best: bool,

    /// Require same annotation category on both sides (single-mode).
    #[arg(long, default_value_t = false)]
    category_constraint: bool,

    /// Archetype output TSV (single-mode).
    #[arg(long)]
    output: Option<PathBuf>,

    /// Enable sweep mode; expects --metrics and --output-matrix.
    #[arg(long, default_value_t = false)]
    sweep: bool,

    /// Sweep spec: e.g. `jaccard:top_n=[20,40]:tau=[0.10,0.15],cosine:tau=[0.20,0.30]`.
    #[arg(long)]
    metrics: Option<String>,

    /// Emit both category-constrained and unconstrained rows for each sweep cell.
    #[arg(long, default_value_t = false)]
    compare_constrained_vs_unconstrained: bool,

    /// Sweep output matrix TSV.
    #[arg(long)]
    output_matrix: Option<PathBuf>,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum SingleMetric {
    Jaccard,
    Cosine,
    /// Cosine after mean-centring each loading vector. Use this when
    /// comparing non-negative loadings (NMF) against signed ones (ICA):
    /// plain cosine has a positivity floor on non-negative vectors, so
    /// one shared threshold is selective for signed loadings and inert
    /// for non-negative ones.
    #[value(name = "cosine-centered", alias = "cosine-centred")]
    CosineCentered,
    Spearman,
}

impl From<SingleMetric> for AlignMetric {
    fn from(m: SingleMetric) -> Self {
        match m {
            SingleMetric::Jaccard => AlignMetric::Jaccard,
            SingleMetric::Cosine => AlignMetric::Cosine,
            SingleMetric::CosineCentered => AlignMetric::CosineCentered,
            SingleMetric::Spearman => AlignMetric::Spearman,
        }
    }
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Programs(args) => run_programs(args),
        Command::Bootstrap(args) => run_bootstrap(args),
        Command::Project(args) => run_project(args),
    }
}

#[derive(ClapArgs, Debug)]
pub struct BootstrapArgs {
    /// Comma-separated list of canonical Atman directories, one per
    /// cohort. Each must contain `measurements.tsv`, `samples.tsv`,
    /// `proteins.tsv`. All cohorts share the same protein universe
    /// (intersection across cohorts is enforced).
    #[arg(long)]
    cohorts: String,

    /// Comma-separated cohort labels (same order as --cohorts).
    /// Defaults to directory basenames.
    #[arg(long)]
    labels: Option<String>,

    /// Number of ICA components per cohort.
    #[arg(long)]
    k: usize,

    /// Number of bootstrap iterations.
    #[arg(long, default_value_t = 100)]
    n_boot: usize,

    /// Top-level seed. Per-iteration sub-seeds derive deterministically
    /// from SplitMix64`(seed, iter)`. The default `20260418` is an
    /// arbitrary fixed integer chosen for reproducibility — change it
    /// only when you intend to reseed an analysis (any other value is
    /// equally valid; keeping it pinned means re-runs without `--seed`
    /// reproduce prior outputs bit-for-bit).
    #[arg(long, default_value_t = 20260418)]
    seed: u64,

    /// Cosine similarity threshold for the alignment step that groups
    /// bootstrap programs into archetypes.
    #[arg(long, default_value_t = 0.30)]
    cosine_tau: f64,

    /// Floor cosine similarity between a point-estimate archetype's
    /// representative loading and a bootstrap program for that
    /// bootstrap archetype to be counted as a match.
    #[arg(long, default_value_t = 0.50)]
    match_tau: f64,

    /// Top-N loadings used by the Jaccard branch of the alignment
    /// (unused for cosine metric; kept for API symmetry).
    #[arg(long, default_value_t = 40)]
    top_n: usize,

    /// FastICA max iterations per seed per cohort.
    #[arg(long, default_value_t = 300)]
    max_iter: usize,

    /// FastICA convergence tolerance.
    #[arg(long, default_value_t = 1e-4)]
    tol: f64,

    /// Refuse when any cohort has fewer than this many subjects.
    #[arg(long, default_value_t = 20)]
    min_subjects: usize,

    /// Similarity metric used both to group bootstrap programs into
    /// archetypes and to match bootstrap/jackknife archetypes back
    /// to the point estimate. `cosine` (default) is fastest and
    /// sign-invariant; `jaccard` uses top-`--top-n` loading sets;
    /// `spearman` uses absolute rank correlation.
    #[arg(long, default_value = "cosine")]
    metric: String,

    /// Drop assays with more than this fraction of missing samples
    /// per cohort (matches `decompose ica`).
    #[arg(long, default_value_t = 0.0)]
    max_missing_fraction: f64,

    /// Imputation for residual missingness: `none` (fail), `mean`.
    #[arg(long, default_value = "none")]
    impute: String,

    /// Two-sided alpha for the percentile and BCa CIs on the
    /// bootstrap `n_cohorts` distribution. Default `0.05` ⇒ 95% CI.
    /// Must lie in `(0, 1)`.
    #[arg(long, default_value_t = 0.05)]
    ci_alpha: f64,

    /// Per-resample decomposition method. `ica` (default) is FastICA,
    /// byte-identical to prior releases. `nmf` runs single-seed
    /// multiplicative-updates NMF per resample instead — applied
    /// identically to the point estimate, every bootstrap resample,
    /// and every jackknife replicate.
    #[arg(long, default_value = "ica")]
    decomposition: String,

    /// Beta-divergence loss for `--decomposition nmf`: `frobenius`
    /// (default) or `kullback-leibler` (alias `kl`). Ignored for `ica`.
    #[arg(long, default_value = "frobenius")]
    beta_loss: String,

    /// NMF initialization strategy for `--decomposition nmf`: `random`
    /// (default, seeded from the same SplitMix64 derivation the ICA
    /// branch uses) or `nndsvda`. Ignored for `ica`.
    #[arg(long, default_value = "random")]
    init: String,

    /// NMF max multiplicative-update iterations per seed per cohort.
    /// Only used with `--decomposition nmf`.
    #[arg(long, default_value_t = 500)]
    nmf_max_iter: usize,

    /// NMF convergence tolerance. Only used with `--decomposition nmf`.
    #[arg(long, default_value_t = 1e-5)]
    nmf_tol: f64,

    /// Worker threads for the bootstrap and jackknife loops. `0`
    /// (default) uses one thread per available core. Output bytes do
    /// not depend on this value — every iteration has its own sub-seed
    /// and results are folded in iteration order — so it is safe to
    /// vary between a pinned run and its replay; the sidecar records it.
    #[arg(long, default_value_t = 0)]
    threads: usize,

    /// Pre-decomposition transform for `--decomposition nmf`: `none`
    /// (default; requires non-negative input — rejected loudly
    /// otherwise), `exp2-clip`, or `shift-min`. Re-applied fresh to
    /// every per-resample matrix (point estimate, each bootstrap
    /// resample, each jackknife replicate) rather than cached once.
    /// Ignored for `ica`.
    #[arg(long, default_value = "none")]
    transform: String,

    /// Clamp radius `c` for `--transform exp2-clip`
    /// (`2^clamp(x, -c, +c)`). Only valid together with `--transform
    /// exp2-clip`; passing it with any other `--transform` is a hard
    /// error. Default 6.0 when omitted.
    #[arg(long)]
    transform_clamp: Option<f64>,

    /// Output TSV path. Summary with one row per point-estimate
    /// archetype.
    #[arg(long)]
    output: PathBuf,
}

fn run_bootstrap(args: BootstrapArgs) -> Result<()> {
    let started_at = SystemTime::now();
    if args.k == 0 {
        bail!("--k must be >= 1");
    }
    if args.n_boot == 0 {
        bail!("--n-boot must be >= 1");
    }
    if args.min_subjects < 2 {
        bail!("--min-subjects must be >= 2");
    }
    if !(0.0..=1.0).contains(&args.cosine_tau) {
        bail!("--cosine-tau must be in [0, 1]");
    }
    if !(0.0..=1.0).contains(&args.match_tau) {
        bail!("--match-tau must be in [0, 1]");
    }
    if !(args.ci_alpha.is_finite() && args.ci_alpha > 0.0 && args.ci_alpha < 1.0) {
        bail!("--ci-alpha must lie in (0, 1)");
    }
    let cohort_dirs: Vec<PathBuf> = args
        .cohorts
        .split(',')
        .map(|s| PathBuf::from(s.trim()))
        .filter(|p| !p.as_os_str().is_empty())
        .collect();
    if cohort_dirs.len() < 2 {
        bail!("--cohorts must list at least 2 cohort directories");
    }
    let labels: Vec<String> = match &args.labels {
        Some(s) => s.split(',').map(|s| s.trim().to_string()).collect(),
        None => cohort_dirs
            .iter()
            .map(|p| file_stem(p).to_string())
            .collect(),
    };
    if labels.len() != cohort_dirs.len() {
        bail!(
            "--labels has {} entries but --cohorts has {}",
            labels.len(),
            cohort_dirs.len()
        );
    }
    let impute_mean = match args.impute.as_str() {
        "none" => false,
        "mean" => true,
        other => bail!("--impute {other:?}; expected none or mean"),
    };

    let matrices: Vec<CohortMatrix> = cohort_dirs
        .iter()
        .zip(labels.iter())
        .map(|(dir, label)| {
            // bootstrap aligns programs across cohorts and does not
            // label individual samples, so the per-cohort sample_order
            // is not needed here.
            load_cohort_matrix(dir, label, args.max_missing_fraction, impute_mean)
                .map(|(m, _sample_order)| m)
        })
        .collect::<Result<Vec<_>>>()?;
    // Enforce common protein universe across cohorts by intersecting
    // their label sets and restricting each matrix to the
    // intersection in a canonical order.
    let matrices = intersect_cohorts(matrices)?;

    eprintln!(
        "align bootstrap: cohorts={} k={} n_boot={} cosine_tau={} match_tau={} seed={} threads={}",
        matrices.len(),
        args.k,
        args.n_boot,
        args.cosine_tau,
        args.match_tau,
        args.seed,
        if args.threads == 0 {
            "all".to_string()
        } else {
            args.threads.to_string()
        }
    );

    let metric = match args.metric.as_str() {
        "cosine" => AlignMetric::Cosine,
        "cosine-centered" | "cosine-centred" => AlignMetric::CosineCentered,
        "jaccard" => AlignMetric::Jaccard,
        "spearman" => AlignMetric::Spearman,
        other => bail!("--metric {other:?}; expected cosine, jaccard, or spearman"),
    };
    // Resolved effective `--transform-clamp`, for the sidecar's provenance
    // record (the RESOLVED parameter set — `Some(DEFAULT_EXP2_CLIP_CLAMP)`
    // when `--transform exp2-clip` is used without an explicit clamp, not
    // the raw possibly-`None` flag value). Stays `None` for `--transform
    // none`/`shift-min` (no clamp is in effect) and for `--decomposition
    // ica` (transform is not applied at all).
    let mut resolved_transform_clamp: Option<f64> = None;
    let decomposition = match args.decomposition.as_str() {
        "ica" => Decomposition::Ica {
            max_iter: args.max_iter,
            tol: args.tol,
        },
        "nmf" => {
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
            // Same `--transform-clamp` validation rule as `decompose nmf`:
            // explicit value rejected unless `--transform exp2-clip`.
            if args.transform_clamp.is_some() && args.transform != "exp2-clip" {
                bail!(
                    "--transform-clamp is only valid together with --transform exp2-clip \
                     (got --transform {:?})",
                    args.transform
                );
            }
            if let Some(c) = args.transform_clamp {
                if !(c.is_finite() && c > 0.0) {
                    bail!("--transform-clamp must be finite and > 0, got {}", c);
                }
            }
            let transform = match args.transform.as_str() {
                "none" => NmfTransform::None,
                "exp2-clip" => {
                    let clamp = args.transform_clamp.unwrap_or(DEFAULT_EXP2_CLIP_CLAMP);
                    resolved_transform_clamp = Some(clamp);
                    NmfTransform::Exp2Clip { clamp }
                }
                "shift-min" => NmfTransform::ShiftMin,
                other => bail!(
                    "--transform {:?}: expected `none`, `exp2-clip`, or `shift-min`",
                    other
                ),
            };
            Decomposition::Nmf {
                beta_loss,
                init,
                max_iter: args.nmf_max_iter,
                tol: args.nmf_tol,
                transform,
            }
        }
        other => bail!("--decomposition {other:?}; expected ica or nmf"),
    };
    let params = BootstrapParams {
        k: args.k,
        n_boot: args.n_boot,
        seed: args.seed,
        top_n: args.top_n,
        cosine_tau: args.cosine_tau,
        match_tau: args.match_tau,
        decomposition,
        min_subjects: args.min_subjects,
        metric,
        ci_alpha: args.ci_alpha,
        threads: args.threads,
    };
    let rows: Vec<BootstrapRow> =
        align_bootstrap(&matrices, params).map_err(|e| anyhow::anyhow!(e))?;

    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output dir {:?}", parent))?;
    }
    write_bootstrap_summary(&args.output, &rows)?;
    eprintln!(
        "align bootstrap: wrote {} archetypes to {}",
        rows.len(),
        args.output.display()
    );

    let finished_at = SystemTime::now();
    // Hash each cohort's canonical input directory separately.
    let mut labeled: Vec<(String, PathBuf)> = Vec::new();
    for (dir, label) in cohort_dirs.iter().zip(labels.iter()) {
        let file = "measurements.tsv";
        labeled.push((format!("{}_{}", label, file), dir.join(file)));
        labeled.push((format!("{}_samples", label), dir.join("samples.tsv")));
        labeled.push((format!("{}_proteins", label), dir.join("proteins.tsv")));
    }
    let refs: Vec<(&str, &Path)> = labeled
        .iter()
        .map(|(l, p)| (l.as_str(), p.as_path()))
        .collect();
    let inputs_sha256 = hash_labeled_inputs(&refs)?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "align bootstrap",
        json!({
            "cohorts": args.cohorts,
            "labels": args.labels,
            "k": args.k,
            "n-boot": args.n_boot,
            "seed": args.seed,
            "metric": args.metric,
            "cosine-tau": args.cosine_tau,
            "match-tau": args.match_tau,
            "top-n": args.top_n,
            "max-iter": args.max_iter,
            "tol": args.tol,
            "min-subjects": args.min_subjects,
            "max-missing-fraction": args.max_missing_fraction,
            "impute": args.impute,
            "ci-alpha": args.ci_alpha,
            "beta-loss": args.beta_loss,
            "init": args.init,
            "nmf-max-iter": args.nmf_max_iter,
            "nmf-tol": args.nmf_tol,
            "threads": args.threads,
            "transform": args.transform,
            // Resolved effective clamp (6.0 default when --transform
            // exp2-clip is used without an explicit --transform-clamp),
            // not the raw flag — the sidecar records the resolved
            // parameter set. null for none/shift-min/ica.
            "transform-clamp": resolved_transform_clamp,
            "output": args.output.display().to_string(),
            // bootstrap always runs its chosen decomposition (ica|nmf)
            // internally per resample — no external loadings TSV to sniff.
            "decomposition_method": args.decomposition,
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("align bootstrap: sidecar={}", sidecar.display());
    Ok(())
}

/// Returns the cohort abundance matrix together with the `sample_order`
/// that its rows are indexed by. `sample_order` is the single source of
/// truth for which `sample_id` each matrix row (and any downstream
/// projected activation) belongs to: it is built from the MEASUREMENTS
/// file in first-seen order over non-QC, finite rows — NOT from
/// samples.tsv. Callers that need to label rows MUST use this vector
/// rather than re-deriving an order from samples.tsv.
fn load_cohort_matrix(
    dir: &Path,
    label: &str,
    max_missing_fraction: f64,
    impute_mean: bool,
) -> Result<(CohortMatrix, Vec<String>)> {
    let file = "measurements.tsv";
    let records = read_measurements_long(&dir.join(file))?;
    if records.is_empty() {
        bail!("no measurements in {:?}", dir.join(file));
    }
    // Abundance by (assay, sample). Only rows with finite values not
    // dropped by QC contribute.
    let mut abundance: BTreeMap<(String, String), f64> = BTreeMap::new();
    let mut sample_order: Vec<String> = Vec::new();
    let mut seen_samples: BTreeSet<String> = BTreeSet::new();
    let mut assays: BTreeSet<String> = BTreeSet::new();
    for r in &records {
        if r.dropped_by_qc {
            continue;
        }
        let v = r.abundance.as_f64();
        if !v.is_finite() {
            continue;
        }
        let sid = r.sample_id.clone();
        let aid = r.assay_id.0.clone();
        if seen_samples.insert(sid.clone()) {
            sample_order.push(sid.clone());
        }
        assays.insert(aid.clone());
        abundance.insert((aid, sid), v);
    }
    if sample_order.len() < 2 {
        bail!("cohort {label:?}: need >= 2 samples");
    }
    let n = sample_order.len();
    let mut kept: Vec<String> = Vec::new();
    for a in &assays {
        let present = sample_order
            .iter()
            .filter(|s| abundance.contains_key(&(a.clone(), (*s).clone())))
            .count();
        let missing_frac = 1.0 - present as f64 / n as f64;
        if missing_frac <= max_missing_fraction + 1e-12 {
            kept.push(a.clone());
        }
    }
    if kept.is_empty() {
        bail!(
            "cohort {label:?}: no assays retained at --max-missing-fraction={max_missing_fraction}"
        );
    }

    let mut data = vec![vec![0.0_f64; kept.len()]; n];
    for (j, a) in kept.iter().enumerate() {
        let mut vals: Vec<Option<f64>> = Vec::with_capacity(n);
        let mut sum = 0.0_f64;
        let mut count = 0usize;
        for s in &sample_order {
            let v = abundance.get(&(a.clone(), s.clone())).copied();
            if let Some(x) = v {
                sum += x;
                count += 1;
            }
            vals.push(v);
        }
        let mean = if count > 0 { sum / count as f64 } else { 0.0 };
        for (i, v) in vals.iter().enumerate() {
            match v {
                Some(x) => data[i][j] = *x,
                None => {
                    if impute_mean {
                        data[i][j] = mean;
                    } else {
                        bail!(
                            "cohort {label:?}: assay {} missing sample {}; rerun with --impute mean",
                            a,
                            sample_order[i]
                        );
                    }
                }
            }
        }
    }
    Ok((
        CohortMatrix {
            label: label.to_string(),
            data,
            protein_labels: kept,
        },
        sample_order,
    ))
}

/// Restrict each cohort's matrix to the intersection of protein label
/// sets, in a canonical (sorted) order.
fn intersect_cohorts(mut matrices: Vec<CohortMatrix>) -> Result<Vec<CohortMatrix>> {
    if matrices.is_empty() {
        return Ok(matrices);
    }
    let mut universe: BTreeSet<String> = matrices[0].protein_labels.iter().cloned().collect();
    for m in matrices.iter().skip(1) {
        let s: BTreeSet<String> = m.protein_labels.iter().cloned().collect();
        universe = universe.intersection(&s).cloned().collect();
    }
    if universe.is_empty() {
        bail!("cohorts share zero proteins after intersection");
    }
    let ordered: Vec<String> = {
        let mut v: Vec<String> = universe.into_iter().collect();
        v.sort();
        v
    };
    for m in &mut matrices {
        let index_by_label: BTreeMap<String, usize> = m
            .protein_labels
            .iter()
            .enumerate()
            .map(|(i, l)| (l.clone(), i))
            .collect();
        let mut new_data = vec![vec![0.0_f64; ordered.len()]; m.data.len()];
        for (j_new, label) in ordered.iter().enumerate() {
            let j_old = index_by_label[label];
            for (i, row) in m.data.iter().enumerate() {
                new_data[i][j_new] = row[j_old];
            }
        }
        m.data = new_data;
        m.protein_labels = ordered.clone();
    }
    Ok(matrices)
}

fn write_bootstrap_summary(path: &Path, rows: &[BootstrapRow]) -> Result<()> {
    let mut out = String::from(
        "archetype_id\tobserved_n_cohorts\tobserved_cohorts\t\
         bootstrap_mean_n_cohorts\tbootstrap_prob_universal\tbootstrap_prob_multi\t\
         ci_lower_n_cohorts\tci_upper_n_cohorts\tbootstrap_match_rate\t\
         alignment_entropy\tbca_lower_n_cohorts\tbca_upper_n_cohorts\t\
         bca_fallback_to_percentile\n",
    );
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{:.6}\t\
             {:.6}\t{:.6}\t{:.6}\t{}\n",
            r.archetype_id,
            r.observed_n_cohorts,
            r.observed_cohorts.join(","),
            r.bootstrap_mean_n_cohorts,
            r.bootstrap_prob_universal,
            r.bootstrap_prob_multi,
            r.ci_lower_n_cohorts,
            r.ci_upper_n_cohorts,
            r.bootstrap_match_rate,
            r.alignment_entropy,
            r.bca_lower_n_cohorts,
            r.bca_upper_n_cohorts,
            if r.bca_fallback_to_percentile { 1 } else { 0 },
        ));
    }
    atomic_write(path, out.as_bytes())
}

fn run_programs(args: ProgramsArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let paths: Vec<PathBuf> = args
        .loadings
        .split(',')
        .map(|s| PathBuf::from(s.trim()))
        .collect();
    if paths.is_empty() {
        bail!("--loadings must list at least one TSV");
    }
    let cohort_labels: Vec<String> = match &args.cohorts {
        Some(c) => c.split(',').map(|s| s.trim().to_string()).collect(),
        None => paths.iter().map(|p| file_stem(p).to_string()).collect(),
    };
    if cohort_labels.len() != paths.len() {
        bail!(
            "--cohorts has {} labels but --loadings has {} paths",
            cohort_labels.len(),
            paths.len()
        );
    }
    let annotation_paths: Vec<Option<PathBuf>> = match &args.annotations {
        Some(a) => {
            let parts: Vec<PathBuf> = a.split(',').map(|s| PathBuf::from(s.trim())).collect();
            if parts.len() != paths.len() {
                bail!(
                    "--annotations has {} paths but --loadings has {} paths",
                    parts.len(),
                    paths.len()
                );
            }
            parts.into_iter().map(Some).collect()
        }
        None => vec![None; paths.len()],
    };

    // Pass 1: build union of labels across all cohorts.
    let mut loadings_raw: Vec<Vec<(String, String, f64)>> = Vec::with_capacity(paths.len());
    for path in &paths {
        loadings_raw.push(read_loadings(path, &args.label_col)?);
    }
    let mut label_index: BTreeMap<String, usize> = BTreeMap::new();
    for cohort in &loadings_raw {
        for (_prog, label, _val) in cohort {
            if !label_index.contains_key(label) {
                let idx = label_index.len();
                label_index.insert(label.clone(), idx);
            }
        }
    }
    let universe_size = label_index.len();
    if universe_size == 0 {
        bail!("no labels found across cohorts");
    }

    // Assemble AlignedProgram structs across cohorts.
    let mut programs: Vec<AlignedProgram> = Vec::new();
    for (ci, cohort_rows) in loadings_raw.iter().enumerate() {
        let cohort_label = &cohort_labels[ci];
        let annotation_map: BTreeMap<String, String> = match &annotation_paths[ci] {
            Some(ap) => read_annotations(ap)?,
            None => BTreeMap::new(),
        };
        let mut by_program: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        for (program, label, value) in cohort_rows {
            let vec = by_program
                .entry(program.clone())
                .or_insert_with(|| vec![0.0_f64; universe_size]);
            if let Some(idx) = label_index.get(label) {
                vec[*idx] = *value;
            }
        }
        for (program, values) in by_program {
            programs.push(AlignedProgram {
                cohort: cohort_label.clone(),
                program: program.clone(),
                category: annotation_map.get(&program).cloned(),
                values,
            });
        }
    }

    let primary_output = if args.sweep {
        let matrix_path = args
            .output_matrix
            .as_ref()
            .context("sweep mode requires --output-matrix")?;
        let metrics_spec = args
            .metrics
            .as_deref()
            .context("sweep mode requires --metrics")?;
        let grid = parse_sweep_spec(metrics_spec)?;
        run_sweep(
            &programs,
            &grid,
            args.reciprocal_best,
            args.compare_constrained_vs_unconstrained,
            matrix_path,
        )?;
        matrix_path.clone()
    } else {
        let output = args
            .output
            .as_ref()
            .context("single-mode requires --output")?;
        // A tau that admits nearly every candidate pair is not
        // thresholding anything: the structure is then decided entirely
        // by reciprocal-best matching, and a sweep over tau looks robust
        // when it is only saturated. Non-negative loadings (NMF) hit
        // this because cosine has a positivity floor on them.
        let (n_pairs, n_above, median_sim) =
            atman_core::align::tau_selectivity(&programs, args.metric.into(), args.top_n, args.tau);
        if n_pairs > 0 {
            let frac = n_above as f64 / n_pairs as f64;
            eprintln!(
                "align programs: tau selectivity: {n_above}/{n_pairs} cross-cohort pairs \
                 ({:.1}%) are at or above --tau {}; median similarity {:.3}",
                100.0 * frac,
                args.tau,
                median_sim,
            );
            if frac > 0.90 {
                eprintln!(
                    "align programs: warning: --tau {} admits {:.1}% of cross-cohort pairs, so it \
                     is not thresholding — the archetype structure is determined by \
                     reciprocal-best matching alone and a sweep over tau will look robust \
                     because it is saturated. Non-negative loadings hit this because cosine \
                     cannot go below zero on them; --metric cosine-centered removes that floor.",
                    args.tau,
                    100.0 * frac,
                );
            }
        }
        let labels = build_archetypes(
            &programs,
            args.metric.into(),
            args.top_n,
            args.tau,
            args.reciprocal_best,
            args.category_constraint,
        );
        write_archetypes(output, &programs, &labels)?;
        let summary = summarize_archetypes(&programs, &labels);
        eprintln!(
            "align programs: archetypes={} multi_cohort={} universal={} category_recovery={:.3}",
            summary.n_archetypes,
            summary.n_multi_cohort,
            summary.n_universal,
            summary.category_recovery
        );
        output.clone()
    };

    let finished_at = SystemTime::now();
    let mut labeled: Vec<(String, PathBuf)> = Vec::new();
    for (ci, path) in paths.iter().enumerate() {
        labeled.push((format!("loadings_{}", cohort_labels[ci]), path.clone()));
    }
    for (ci, ap) in annotation_paths.iter().enumerate() {
        if let Some(p) = ap {
            labeled.push((format!("annotations_{}", cohort_labels[ci]), p.clone()));
        }
    }
    let refs: Vec<(&str, &Path)> = labeled
        .iter()
        .map(|(l, p)| (l.as_str(), p.as_path()))
        .collect();
    let inputs_sha256 = hash_labeled_inputs(&refs)?;
    let sidecar = sidecar_path_for(&primary_output);
    let metric = match args.metric {
        SingleMetric::Jaccard => "jaccard",
        SingleMetric::Cosine => "cosine",
        SingleMetric::CosineCentered => "cosine-centered",
        SingleMetric::Spearman => "spearman",
    };
    let decomposition_method = sniff_decomposition_method(&paths);
    write_run_sidecar(
        &sidecar,
        "align programs",
        json!({
            "loadings": args.loadings,
            "cohorts": args.cohorts,
            "annotations": args.annotations,
            "label-col": args.label_col,
            "metric": metric,
            "top-n": args.top_n,
            "tau": args.tau,
            "reciprocal-best": args.reciprocal_best,
            "category-constraint": args.category_constraint,
            "output": args.output.as_ref().map(|p| p.display().to_string()),
            "sweep": args.sweep,
            "metrics": args.metrics,
            "compare-constrained-vs-unconstrained": args.compare_constrained_vs_unconstrained,
            "output-matrix": args.output_matrix.as_ref().map(|p| p.display().to_string()),
            "decomposition_method": decomposition_method,
        }),
        &inputs_sha256,
        std::slice::from_ref(&primary_output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("align programs: sidecar={}", sidecar.display());
    Ok(())
}

fn read_loadings(path: &Path, label_col: &str) -> Result<Vec<(String, String, f64)>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let program_col = headers
        .iter()
        .position(|h| h == "program")
        .with_context(|| format!("{:?} missing 'program' column", path))?;
    let loading_col = headers
        .iter()
        .position(|h| h == "loading")
        .with_context(|| format!("{:?} missing 'loading' column", path))?;
    let label_idx = headers
        .iter()
        .position(|h| h == label_col)
        .or_else(|| headers.iter().position(|h| h == "gene_symbol"))
        .or_else(|| headers.iter().position(|h| h == "assay_id"))
        .or_else(|| headers.iter().position(|h| h == "protein"))
        .with_context(|| format!("{:?} missing label column", path))?;
    let mut out = Vec::new();
    for row in reader.records() {
        let row = row?;
        let value: f64 = row[loading_col]
            .parse()
            .with_context(|| format!("parsing loading in {:?}", path))?;
        if !value.is_finite() {
            continue;
        }
        let label = row[label_idx].trim();
        if label.is_empty() {
            continue;
        }
        out.push((row[program_col].to_string(), label.to_string(), value));
    }
    if out.is_empty() {
        bail!("no usable loadings in {:?}", path);
    }
    Ok(out)
}

fn read_annotations(path: &Path) -> Result<BTreeMap<String, String>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let program_col = headers
        .iter()
        .position(|h| h == "program")
        .with_context(|| format!("{:?} missing 'program' column", path))?;
    let category_col = headers
        .iter()
        .position(|h| h == "category")
        .with_context(|| format!("{:?} missing 'category' column", path))?;
    let mut out = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let program = row[program_col].to_string();
        let category = row[category_col].trim().to_string();
        if !category.is_empty() {
            out.insert(program, category);
        }
    }
    Ok(out)
}

fn file_stem(path: &Path) -> &str {
    path.file_stem()
        .and_then(|os| os.to_str())
        .unwrap_or("cohort")
}

fn write_archetypes(path: &Path, programs: &[AlignedProgram], labels: &[usize]) -> Result<()> {
    // Assign compact archetype ids: sort members by (cohort, program), group by label.
    let mut by_label: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (i, lab) in labels.iter().enumerate() {
        by_label.entry(*lab).or_default().push(i);
    }
    let mut archetype_id_map: BTreeMap<usize, usize> = BTreeMap::new();
    let mut next_id = 1usize;
    // Assign multi-member archetypes first in deterministic order.
    for (lab, members) in by_label.iter() {
        if members.len() >= 2 {
            archetype_id_map.insert(*lab, next_id);
            next_id += 1;
        }
    }
    for (lab, members) in by_label.iter() {
        if members.len() < 2 && !archetype_id_map.contains_key(lab) {
            archetype_id_map.insert(*lab, next_id);
            next_id += 1;
        }
    }
    let mut out = String::from(
        "archetype_id\tcohort\tprogram\tcategory\tn_members\tn_cohorts\tis_singleton\n",
    );
    let mut archetype_sizes: BTreeMap<usize, (usize, BTreeSet<String>)> = BTreeMap::new();
    for (lab, members) in by_label.iter() {
        let cohorts: BTreeSet<String> = members
            .iter()
            .map(|i| programs[*i].cohort.clone())
            .collect();
        archetype_sizes.insert(*lab, (members.len(), cohorts));
    }
    let mut rows: Vec<(usize, &AlignedProgram, usize, usize, bool)> = Vec::new();
    for (lab, members) in by_label.iter() {
        let id = archetype_id_map[lab];
        let (n_members, cohorts) = &archetype_sizes[lab];
        for &i in members {
            rows.push((id, &programs[i], *n_members, cohorts.len(), *n_members < 2));
        }
    }
    rows.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then(a.1.cohort.cmp(&b.1.cohort))
            .then(a.1.program.cmp(&b.1.program))
    });
    for (id, prog, n_members, n_cohorts, singleton) in rows {
        out.push_str(&format!(
            "A{:04}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            id,
            prog.cohort,
            prog.program,
            prog.category.clone().unwrap_or_default(),
            n_members,
            n_cohorts,
            if singleton { 1 } else { 0 }
        ));
    }
    atomic_write(path, out.as_bytes())
}

#[derive(Debug, Clone)]
struct SweepCell {
    metric: AlignMetric,
    top_n: usize,
    tau: f64,
}

fn split_top_level(input: &str, separator: char) -> Vec<String> {
    let mut out = Vec::new();
    let mut depth = 0i32;
    let mut current = String::new();
    for c in input.chars() {
        match c {
            '[' | '(' | '{' => {
                depth += 1;
                current.push(c);
            }
            ']' | ')' | '}' => {
                depth -= 1;
                current.push(c);
            }
            other if other == separator && depth == 0 => {
                out.push(std::mem::take(&mut current));
            }
            other => current.push(other),
        }
    }
    if !current.is_empty() {
        out.push(current);
    }
    out
}

fn parse_sweep_spec(spec: &str) -> Result<Vec<SweepCell>> {
    let mut out = Vec::new();
    for block in split_top_level(spec, ',') {
        let block = block.trim();
        if block.is_empty() {
            continue;
        }
        let parts = split_top_level(block, ':');
        let metric_name = parts.first().cloned().unwrap_or_default();
        let metric = AlignMetric::parse(&metric_name)
            .with_context(|| format!("unknown metric {:?} in --metrics", metric_name))?;
        let mut top_ns = vec![40usize];
        let mut taus: Vec<f64> = Vec::new();
        for kv in parts.iter().skip(1) {
            let kv = kv.trim();
            if kv.is_empty() {
                continue;
            }
            let (key, raw_values) = kv
                .split_once('=')
                .with_context(|| format!("expected key=value in {:?}", kv))?;
            let inner = raw_values
                .trim()
                .trim_start_matches('[')
                .trim_end_matches(']');
            match key.trim() {
                "top_n" => {
                    top_ns = inner
                        .split(',')
                        .map(|s| s.trim().parse())
                        .collect::<std::result::Result<Vec<usize>, _>>()
                        .with_context(|| format!("parsing top_n list in {:?}", kv))?;
                }
                "tau" => {
                    taus = inner
                        .split(',')
                        .map(|s| s.trim().parse())
                        .collect::<std::result::Result<Vec<f64>, _>>()
                        .with_context(|| format!("parsing tau list in {:?}", kv))?;
                }
                other => bail!("unknown sweep key {:?}", other),
            }
        }
        if taus.is_empty() {
            bail!("metric {:?} has no tau values", metric_name);
        }
        for &tn in &top_ns {
            for &tau in &taus {
                out.push(SweepCell {
                    metric,
                    top_n: tn,
                    tau,
                });
            }
        }
    }
    if out.is_empty() {
        bail!("--metrics produced no sweep cells");
    }
    Ok(out)
}

fn run_sweep(
    programs: &[AlignedProgram],
    grid: &[SweepCell],
    reciprocal_best: bool,
    compare_constrained: bool,
    output: &Path,
) -> Result<()> {
    let mut out = String::from(
        "metric\ttop_n\ttau\treciprocal_best\tcategory_constraint\tn_archetypes\tn_multi_cohort\tn_universal\tcategory_recovery\n",
    );
    let modes: Vec<(bool, &str)> = if compare_constrained {
        vec![(true, "1"), (false, "0")]
    } else {
        vec![(false, "0")]
    };
    for cell in grid {
        for (cat_flag, cat_flag_str) in &modes {
            let labels = build_archetypes(
                programs,
                cell.metric,
                cell.top_n,
                cell.tau,
                reciprocal_best,
                *cat_flag,
            );
            let summary = summarize_archetypes(programs, &labels);
            let cat_recovery = if summary.category_recovery.is_finite() {
                format!("{:.6}", summary.category_recovery)
            } else {
                "NaN".to_string()
            };
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                cell.metric.as_str(),
                cell.top_n,
                cell.tau,
                if reciprocal_best { 1 } else { 0 },
                cat_flag_str,
                summary.n_archetypes,
                summary.n_multi_cohort,
                summary.n_universal,
                cat_recovery,
            ));
        }
    }
    atomic_write(output, out.as_bytes())?;
    eprintln!(
        "align programs sweep: cells={} output={}",
        grid.len() * modes.len(),
        output.display()
    );
    Ok(())
}

#[derive(ClapArgs, Debug)]
pub struct ProjectArgs {
    /// Path to the archetype TSV emitted by a prior `align programs`.
    #[arg(long)]
    atlas_archetypes: PathBuf,

    /// Comma-separated `cohort_label=loadings_tsv` list — the same
    /// per-cohort loading TSVs that were fed to the `align programs`
    /// run that produced `--atlas-archetypes`.
    #[arg(long)]
    atlas_loadings: String,

    /// Label column in each atlas loading TSV: `gene_symbol` (default),
    /// `assay_id`, or `protein`.
    #[arg(long, default_value = "gene_symbol")]
    label_col: String,

    /// Canonical Atman directory for the cohort to project.
    #[arg(long)]
    cohort_dir: PathBuf,

    /// Transform applied to the cohort matrix before projection. Must match
    /// the transform used at atlas training time (caller's responsibility to
    /// verify — the atlas's own `transform_applied.json` records the
    /// training-time choice). `none`/`clr`/`alr`/`ratio-anchor` are
    /// compositional transforms; `exp2-clip`/`shift-min` are the NMF input
    /// transforms (`2^clamp(x, -c, +c)` and `x - min(X)` respectively) shared
    /// with `decompose nmf`. `shift-min` always recomputes its shift on THIS
    /// cohort's matrix — it does not reuse the training-time shift.
    #[arg(long, default_value = "none")]
    transform: String,

    /// Reference gene symbol for `alr` / `ratio-anchor` transforms.
    #[arg(long)]
    alr_reference: Option<String>,

    /// Clamp radius `c` for `--transform exp2-clip` (`2^clamp(x, -c, +c)`).
    /// Only valid together with `--transform exp2-clip`; passing it with any
    /// other `--transform` is a hard error. Default 6.0 when omitted.
    #[arg(long)]
    transform_clamp: Option<f64>,

    /// Projection method: `ls` (plain least-squares — fails when the
    /// atlas is rank-deficient) or `ridge` (default, with
    /// `--ridge-lambda`).
    #[arg(long, default_value = "ridge")]
    projection: String,

    /// Ridge penalty on the normal equations. Ignored when
    /// `--projection ls`.
    #[arg(long, default_value_t = 0.01)]
    ridge_lambda: f64,

    /// Drop assays with more than this fraction of missing samples
    /// in the new cohort (before projection).
    #[arg(long, default_value_t = 0.0)]
    max_missing_fraction: f64,

    /// Imputation for residual missingness: `none` (fail) or `mean`.
    #[arg(long, default_value = "none")]
    impute: String,

    /// Output directory for `projected_activations.tsv`,
    /// `projection_qc.tsv`, and the run sidecar.
    #[arg(long)]
    output_dir: PathBuf,
}

fn parse_atlas_loadings(spec: &str) -> Result<Vec<(String, PathBuf)>> {
    let mut out = Vec::new();
    for chunk in spec.split(',') {
        let chunk = chunk.trim();
        if chunk.is_empty() {
            continue;
        }
        let (cohort, path) = chunk
            .split_once('=')
            .with_context(|| format!("--atlas-loadings chunk {chunk:?} is not cohort=path"))?;
        let cohort = cohort.trim().to_string();
        let path = PathBuf::from(path.trim());
        if cohort.is_empty() {
            bail!("--atlas-loadings chunk {chunk:?} has empty cohort label");
        }
        out.push((cohort, path));
    }
    if out.is_empty() {
        bail!("--atlas-loadings resolved to zero entries");
    }
    Ok(out)
}

struct AtlasArchetypeRow {
    archetype_id: String,
    cohort: String,
    program: String,
    n_members: usize,
    n_cohorts: usize,
}

fn read_atlas_archetypes(path: &Path) -> Result<Vec<AtlasArchetypeRow>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let col = |name: &str| {
        headers
            .iter()
            .position(|h| h == name)
            .with_context(|| format!("{:?} missing {name} column", path))
    };
    let id = col("archetype_id")?;
    let cohort = col("cohort")?;
    let program = col("program")?;
    let n_members = col("n_members")?;
    let n_cohorts = col("n_cohorts")?;
    let mut out = Vec::new();
    for row in reader.records() {
        let row = row?;
        out.push(AtlasArchetypeRow {
            archetype_id: row[id].to_string(),
            cohort: row[cohort].to_string(),
            program: row[program].to_string(),
            n_members: row[n_members].parse().unwrap_or(0),
            n_cohorts: row[n_cohorts].parse().unwrap_or(0),
        });
    }
    Ok(out)
}

fn parse_transform(
    name: &str,
    alr_reference: Option<&str>,
    labels: &[String],
) -> Result<Transform> {
    match name {
        "none" => Ok(Transform::None),
        "clr" => Ok(Transform::Clr),
        "alr" => {
            let reference =
                alr_reference.context("--transform alr requires --alr-reference <gene>")?;
            let idx = labels
                .iter()
                .position(|l| l == reference)
                .with_context(|| format!("--alr-reference {reference:?} not in cohort labels"))?;
            Ok(Transform::Alr {
                reference_index: idx,
            })
        }
        "ratio-anchor" => {
            let reference = alr_reference
                .context("--transform ratio-anchor requires --alr-reference <gene>")?;
            let idx = labels
                .iter()
                .position(|l| l == reference)
                .with_context(|| format!("--alr-reference {reference:?} not in cohort labels"))?;
            Ok(Transform::RatioAnchor {
                reference_index: idx,
            })
        }
        "ilr" => Ok(Transform::Ilr),
        other => bail!(
            "unknown --transform {other:?}; expected \
             none|clr|alr|ilr|ratio-anchor|exp2-clip|shift-min"
        ),
    }
}

fn run_project(args: ProjectArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let impute_mean = match args.impute.as_str() {
        "none" => false,
        "mean" => true,
        other => bail!("--impute {other:?}; expected none or mean"),
    };
    let method = match args.projection.as_str() {
        "ls" | "least-squares" => ProjectionMethod::LeastSquares,
        "ridge" => {
            if !(args.ridge_lambda.is_finite() && args.ridge_lambda >= 0.0) {
                bail!("--ridge-lambda must be finite and >= 0");
            }
            ProjectionMethod::Ridge(args.ridge_lambda)
        }
        other => bail!("--projection {other:?}; expected ls or ridge"),
    };

    // Parse atlas-loadings, read each per-cohort loadings TSV into
    // (program → Vec<(label, value)>). Reuse the existing
    // read_loadings helper used by `align programs`.
    let loadings_spec = parse_atlas_loadings(&args.atlas_loadings)?;
    let mut per_cohort_programs: BTreeMap<String, BTreeMap<String, BTreeMap<String, f64>>> =
        BTreeMap::new();
    for (cohort_label, path) in &loadings_spec {
        let rows = read_loadings(path, &args.label_col)?;
        let by_program = per_cohort_programs.entry(cohort_label.clone()).or_default();
        for (program, label, value) in rows {
            by_program.entry(program).or_default().insert(label, value);
        }
    }

    // Read atlas archetypes; for each archetype, keep only rows whose
    // cohort was passed on `--atlas-loadings`.
    let archetype_rows = read_atlas_archetypes(&args.atlas_archetypes)?;
    let known_cohorts: BTreeSet<String> = per_cohort_programs.keys().cloned().collect();
    let mut by_archetype: BTreeMap<String, Vec<(String, String, usize, usize)>> = BTreeMap::new();
    for row in archetype_rows {
        if !known_cohorts.contains(&row.cohort) {
            // Silently drop — the archetype may have member programs
            // from cohorts we didn't pass, in which case we'd be
            // averaging an incomplete set.
            continue;
        }
        by_archetype.entry(row.archetype_id).or_default().push((
            row.cohort,
            row.program,
            row.n_members,
            row.n_cohorts,
        ));
    }
    if by_archetype.is_empty() {
        bail!(
            "no archetypes survived the cohort filter — does --atlas-loadings \
             cover the cohorts in --atlas-archetypes?"
        );
    }

    // Intersection protein universe across the participating atlas
    // cohorts. Any archetype whose members span cohorts outside the
    // intersection contributes only on the intersection labels.
    let mut universe: Option<BTreeSet<String>> = None;
    for programs in per_cohort_programs.values() {
        let mut cohort_labels: BTreeSet<String> = BTreeSet::new();
        for label_map in programs.values() {
            for l in label_map.keys() {
                cohort_labels.insert(l.clone());
            }
        }
        universe = Some(match universe {
            Some(u) => u.intersection(&cohort_labels).cloned().collect(),
            None => cohort_labels,
        });
    }
    let universe: Vec<String> = universe
        .map(|u| {
            let mut v: Vec<String> = u.into_iter().collect();
            v.sort();
            v
        })
        .unwrap_or_default();
    if universe.is_empty() {
        bail!("atlas cohorts share zero protein labels after intersection");
    }

    // Build atlas: one loading row per multi-member archetype, mean
    // across member programs' loading vectors on the intersection
    // universe.
    let mut atlas = Atlas {
        protein_labels: universe.clone(),
        archetype_ids: Vec::new(),
        loadings: Vec::new(),
        n_members: Vec::new(),
        n_cohorts: Vec::new(),
    };
    for (id, members) in &by_archetype {
        // Skip singletons — projection only makes sense for
        // archetypes that are reproducible across cohorts.
        if members.iter().all(|(_, _, n, _)| *n < 2) {
            continue;
        }
        let mut sum = vec![0.0_f64; universe.len()];
        let mut n_present = 0usize;
        for (cohort, program, _, _) in members {
            let labels = match per_cohort_programs.get(cohort).and_then(|p| p.get(program)) {
                Some(m) => m,
                None => continue,
            };
            for (i, label) in universe.iter().enumerate() {
                if let Some(&v) = labels.get(label) {
                    sum[i] += v;
                }
            }
            n_present += 1;
        }
        if n_present == 0 {
            continue;
        }
        let mean: Vec<f64> = sum.iter().map(|s| s / n_present as f64).collect();
        let n_cohort_set: BTreeSet<&str> = members.iter().map(|(c, _, _, _)| c.as_str()).collect();
        atlas.archetype_ids.push(id.clone());
        atlas.loadings.push(mean);
        atlas.n_members.push(n_present);
        atlas.n_cohorts.push(n_cohort_set.len());
    }
    if atlas.archetype_ids.is_empty() {
        bail!(
            "no multi-member archetypes assembled from the atlas inputs — \
             every archetype looked like a singleton on the intersection universe"
        );
    }
    eprintln!(
        "align project: atlas k={} p={} cohorts_used={}",
        atlas.archetype_ids.len(),
        atlas.protein_labels.len(),
        per_cohort_programs.len()
    );

    // Load the new cohort's canonical abundance matrix (subjects × proteins).
    // `load_cohort_matrix` returns proteins keyed by assay_id; the
    // atlas is keyed by `--label-col` (default gene_symbol), so remap
    // the cohort's protein labels via the cohort's proteins.tsv.
    let (mut cohort_matrix, sample_order) = load_cohort_matrix(
        &args.cohort_dir,
        "cohort",
        args.max_missing_fraction,
        impute_mean,
    )?;
    {
        use std::collections::HashMap;
        let proteins_path = args.cohort_dir.join("proteins.tsv");
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(&proteins_path)
            .with_context(|| format!("opening {:?}", proteins_path))?;
        let headers = reader.headers()?.clone();
        let assay_idx = headers
            .iter()
            .position(|h| h == "assay_id")
            .context("cohort proteins.tsv missing assay_id")?;
        let label_idx = headers
            .iter()
            .position(|h| h == args.label_col.as_str())
            .or_else(|| headers.iter().position(|h| h == "gene_symbol"))
            .or_else(|| headers.iter().position(|h| h == "protein"))
            .with_context(|| {
                format!(
                    "cohort proteins.tsv missing label column {:?}",
                    args.label_col
                )
            })?;
        let mut assay_to_label: HashMap<String, String> = HashMap::new();
        for row in reader.records() {
            let row = row?;
            assay_to_label.insert(
                row[assay_idx].to_string(),
                row[label_idx].trim().to_string(),
            );
        }
        let mut remapped: Vec<String> = Vec::with_capacity(cohort_matrix.protein_labels.len());
        let mut keep_mask: Vec<bool> = Vec::with_capacity(cohort_matrix.protein_labels.len());
        for assay_id in &cohort_matrix.protein_labels {
            match assay_to_label.get(assay_id) {
                Some(label) if !label.is_empty() => {
                    remapped.push(label.clone());
                    keep_mask.push(true);
                }
                _ => {
                    // No label resolvable — drop this column from the
                    // cohort universe; the atlas intersection will
                    // naturally exclude it.
                    keep_mask.push(false);
                }
            }
        }
        if keep_mask.iter().all(|k| !k) {
            bail!(
                "cohort proteins.tsv has no rows with a non-empty {:?} column",
                args.label_col
            );
        }
        // Apply the mask to both labels and every row of the matrix.
        cohort_matrix.protein_labels = remapped;
        for row in &mut cohort_matrix.data {
            let filtered: Vec<f64> = row
                .iter()
                .zip(keep_mask.iter())
                .filter_map(|(v, k)| if *k { Some(*v) } else { None })
                .collect();
            *row = filtered;
        }
    }
    // `--transform-clamp` only applies to `exp2-clip`; reject it otherwise
    // (mirrors `decompose nmf --transform-clamp`).
    if args.transform_clamp.is_some() && args.transform != "exp2-clip" {
        bail!(
            "--transform-clamp is only valid together with --transform exp2-clip \
             (got --transform {:?})",
            args.transform
        );
    }
    if let Some(c) = args.transform_clamp {
        if !(c.is_finite() && c > 0.0) {
            bail!("--transform-clamp must be finite and > 0, got {}", c);
        }
    }

    // Apply the transform in the cohort's own protein ordering.
    //
    // `exp2-clip` / `shift-min` are the NMF input transforms (shared with
    // `decompose nmf` via `atman_core::nmf::apply_transform`); every other
    // name is a compositional transform (`atman_core::compositional`). The
    // two families use different core functions (in-place mutation +
    // `TransformRecord` vs. an out-of-place `Result`), so they are branched
    // here rather than unified into one call.
    let (transformed, transform_record): (Vec<Vec<f64>>, NmfTransformRecord) =
        match args.transform.as_str() {
            "exp2-clip" => {
                let clamp = args.transform_clamp.unwrap_or(DEFAULT_EXP2_CLIP_CLAMP);
                let mut data = cohort_matrix.data.clone();
                let record = apply_nmf_transform(&mut data, &NmfTransform::Exp2Clip { clamp });
                (data, record)
            }
            "shift-min" => {
                // Recomputes the shift on THIS cohort's matrix; the plan's
                // intended semantics — no reuse of a training-time shift.
                let mut data = cohort_matrix.data.clone();
                let record = apply_nmf_transform(&mut data, &NmfTransform::ShiftMin);
                (data, record)
            }
            other => {
                let transform = parse_transform(
                    other,
                    args.alr_reference.as_deref(),
                    &cohort_matrix.protein_labels,
                )?;
                let transformed = apply_transform(&cohort_matrix.data, transform)
                    .map_err(|e| anyhow::anyhow!("transform failed: {e}"))?;
                (
                    transformed,
                    NmfTransformRecord {
                        name: other.to_string(),
                        clamp: None,
                        shift: None,
                    },
                )
            }
        };
    // For ILR the ordering changes; we emit a clear error because
    // post-transform labels won't match the atlas universe.
    if args.transform == "ilr" {
        bail!(
            "--transform ilr is not supported in `align project` v1: ILR \
             changes the coordinate ordering so atlas labels would not \
             line up with the transformed cohort columns"
        );
    }
    // Subjects are the original sample_ids for the cohort. Use the
    // `sample_order` returned by `load_cohort_matrix` as the single
    // source of truth: it is exactly the order the matrix rows (and
    // hence the projected activations) are built from. Re-deriving an
    // order from samples.tsv would risk attaching activations to the
    // wrong sample whenever the file order differs from the
    // measurement first-seen order — even when the counts coincide.
    let subject_ids = sample_order;
    if subject_ids.len() != transformed.len() {
        bail!(
            "sample_id count {} != transformed matrix rows {}; \
             this should not happen — report as a bug",
            subject_ids.len(),
            transformed.len()
        );
    }
    // Consistency guard: every sample the matrix was built from must
    // be declared in samples.tsv. This catches the inverse error (a
    // measured sample missing from the sample sheet) without letting
    // samples.tsv dictate the activation labels.
    {
        use std::io::Read;
        let path = args.cohort_dir.join("samples.tsv");
        let mut text = String::new();
        std::fs::File::open(&path)
            .with_context(|| format!("opening {:?}", path))?
            .read_to_string(&mut text)?;
        let mut lines = text.lines();
        let header = lines.next().context("samples.tsv empty")?;
        let cols: Vec<&str> = header.split('\t').collect();
        let sid_idx = cols
            .iter()
            .position(|c| *c == "sample_id")
            .context("samples.tsv missing sample_id")?;
        let declared: BTreeSet<String> = lines
            .filter_map(|line| {
                let row: Vec<&str> = line.split('\t').collect();
                row.get(sid_idx).map(|s| s.to_string())
            })
            .collect();
        if let Some(missing) = subject_ids.iter().find(|s| !declared.contains(*s)) {
            bail!(
                "sample {:?} has measurements but is absent from samples.tsv; \
                 cannot label projected activations",
                missing
            );
        }
    }

    let result = project(
        &atlas,
        &cohort_matrix.protein_labels,
        &subject_ids,
        &transformed,
        method,
    )
    .map_err(|e| anyhow::anyhow!(e))?;

    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating {:?}", args.output_dir))?;
    let activations_path = args.output_dir.join("projected_activations.tsv");
    write_projected_activations(&activations_path, &result)?;
    let qc_path = args.output_dir.join("projection_qc.tsv");
    write_projection_qc(&qc_path, &result)?;
    // Loading-weighted coverage, written beside the per-sample QC.
    let arch_cov_path = qc_path
        .parent()
        .unwrap_or(Path::new("."))
        .join("projection_archetype_coverage.tsv");
    write_archetype_coverage(&arch_cov_path, &result)?;
    eprintln!(
        "align project: wrote {} subjects × {} archetypes; avg coverage = {:.3}",
        result.subject_ids.len(),
        result.archetype_ids.len(),
        if result.qc.is_empty() {
            0.0
        } else {
            result.qc.iter().map(|q| q.coverage_fraction).sum::<f64>() / result.qc.len() as f64
        }
    );
    // Unweighted coverage stays high when every low-loading protein is
    // present and every defining one is absent. Flag the archetypes
    // whose own mass did not survive the join.
    for c in &result.archetype_coverage {
        if c.weighted_coverage < 0.5 || c.top20_present < 10 {
            eprintln!(
                "align project: warning: archetype {} has {:.1}% of its loading mass present in \
                 this cohort and {} of its top-20 defining proteins; the largest absent loading \
                 is {:.0}% of its maximum. The projected score is computed from what remains and \
                 may not be the same quantity the archetype names.",
                c.archetype_id,
                100.0 * c.weighted_coverage,
                c.top20_present,
                100.0 * c.largest_absent_loading_share,
            );
        }
    }

    let finished_at = SystemTime::now();
    let mut labeled: Vec<(String, PathBuf)> = Vec::new();
    labeled.push(("atlas_archetypes".into(), args.atlas_archetypes.clone()));
    for (cohort, path) in &loadings_spec {
        labeled.push((format!("atlas_loadings_{}", cohort), path.clone()));
    }
    let file = "measurements.tsv";
    labeled.push((format!("cohort_{file}"), args.cohort_dir.join(file)));
    labeled.push(("cohort_samples".into(), args.cohort_dir.join("samples.tsv")));
    labeled.push((
        "cohort_proteins".into(),
        args.cohort_dir.join("proteins.tsv"),
    ));
    let refs: Vec<(&str, &Path)> = labeled
        .iter()
        .map(|(l, p)| (l.as_str(), p.as_path()))
        .collect();
    let inputs_sha256 = hash_labeled_inputs(&refs)?;
    let sidecar = sidecar_path_for(&activations_path);
    let atlas_loadings_paths: Vec<PathBuf> = loadings_spec.iter().map(|(_, p)| p.clone()).collect();
    let decomposition_method = sniff_decomposition_method(&atlas_loadings_paths);
    write_run_sidecar(
        &sidecar,
        "align project",
        json!({
            "atlas-archetypes": args.atlas_archetypes.display().to_string(),
            "atlas-loadings": args.atlas_loadings,
            "label-col": args.label_col,
            "cohort-dir": args.cohort_dir.display().to_string(),
            "transform": args.transform,
            "transform_clamp": transform_record.clamp,
            "transform_shift": transform_record.shift,
            "alr-reference": args.alr_reference,
            "projection": args.projection,
            "ridge-lambda": args.ridge_lambda,
            "max-missing-fraction": args.max_missing_fraction,
            "impute": args.impute,
            "atlas-proteins-missing-in-cohort": result.atlas_proteins_missing_in_cohort,
            "atlas-k": atlas.archetype_ids.len(),
            "atlas-p": atlas.protein_labels.len(),
            "decomposition_method": decomposition_method,
        }),
        &inputs_sha256,
        &[activations_path.clone(), qc_path.clone()],
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("align project: sidecar={}", sidecar.display());
    Ok(())
}

fn write_projected_activations(
    path: &Path,
    result: &atman_core::align_project::ProjectionResult,
) -> Result<()> {
    let mut buf = String::from("sample_id");
    for a in &result.archetype_ids {
        buf.push('\t');
        buf.push_str(a);
    }
    buf.push('\n');
    for (si, sid) in result.subject_ids.iter().enumerate() {
        buf.push_str(sid);
        for &v in &result.activations[si] {
            buf.push('\t');
            buf.push_str(&format!("{v:.6}"));
        }
        buf.push('\n');
    }
    atomic_write(path, buf.as_bytes())
}

/// Per-archetype loading-weighted coverage.
///
/// The per-sample `coverage_fraction` counts atlas proteins and weights
/// them equally, so it cannot distinguish a cohort missing a hundred
/// irrelevant proteins from one missing the handful that define an
/// archetype. This table answers the second question.
fn write_archetype_coverage(
    path: &Path,
    result: &atman_core::align_project::ProjectionResult,
) -> Result<()> {
    let mut buf = String::from(
        "archetype_id\tweighted_coverage\ttop20_present\tlargest_absent_loading_share\n",
    );
    for c in &result.archetype_coverage {
        buf.push_str(&format!(
            "{}\t{:.6}\t{}\t{:.6}\n",
            c.archetype_id, c.weighted_coverage, c.top20_present, c.largest_absent_loading_share,
        ));
    }
    atomic_write(path, buf.as_bytes())
}

fn write_projection_qc(
    path: &Path,
    result: &atman_core::align_project::ProjectionResult,
) -> Result<()> {
    let mut buf =
        String::from("sample_id\tresidual_norm\tcoverage_fraction\tn_present\tn_missing\n");
    for q in &result.qc {
        buf.push_str(&format!(
            "{}\t{:.6}\t{:.6}\t{}\t{}\n",
            q.subject_id, q.residual_norm, q.coverage_fraction, q.n_present, q.n_missing,
        ));
    }
    atomic_write(path, buf.as_bytes())
}
