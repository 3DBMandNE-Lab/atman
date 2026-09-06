use super::{k_selection_bound_hit, warn_if_k_selection_hit_bound};
use anyhow::{bail, Context, Result};
use atman_core::compositional::{apply_transform, Transform};
use atman_core::decompose_unmix::{
    bootstrap_ci, ora_enrichment, select_k_auto, unmix, AbundanceCi, AbundanceMethod,
    AnnotationRow, EndmemberMethod, KSweepRow, LoadingCi, UnmixConfig, UnmixResult,
};
use clap::Args as ClapArgs;
use serde_json::json;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{
    atomic_write, hash_canonical_inputs, read_measurements_long, sidecar_path_for,
    write_run_sidecar,
};

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

    /// Endmember extraction method.
    ///
    /// `spa` (default): Successive Projection Algorithm — greedy
    /// max-residual vertex search with no RNG, so the endmember set is
    /// a function of the data alone and `--seed` does not affect it.
    /// `vca`: Vertex Component Analysis, which probes the reduced space
    /// with random directions; correct, but its vertex set can still
    /// move with `--seed` on noisy data, so sweep seeds before
    /// reporting. `nfindr`: iterative simplex-volume maximization,
    /// initialized from VCA and therefore also seed-dependent.
    #[arg(long, default_value = "spa")]
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

pub(super) fn run_unmix(args: UnmixArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let endmember_method = match args.method.as_str() {
        "spa" => EndmemberMethod::Spa,
        "vca" => EndmemberMethod::Vca,
        "nfindr" => EndmemberMethod::Nfindr {
            max_passes: args.nfindr_max_passes,
        },
        other => bail!("--method {other:?}; expected `spa`, `vca` or `nfindr`"),
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
    if args.k == "auto" {
        warn_if_k_selection_hit_bound(
            "decompose unmix",
            "auto (elbow)",
            effective_k,
            args.k_min,
            args.k_max,
        );
    }
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
    if result.n_abundance_not_converged > 0 {
        eprintln!(
            "decompose unmix: warning: {} of {} per-subject abundance solves hit \
             --fcls-max-iter={} without reaching --fcls-tol={:.1e}. Those abundance vectors are \
             where the projected-gradient loop stopped, not the constrained least-squares \
             solution. Raise --fcls-max-iter or relax --fcls-tol.",
            result.n_abundance_not_converged,
            subject_ids.len(),
            args.fcls_max_iter,
            args.fcls_tol,
        );
    }

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
            // `spa` ignores the seed entirely; recording that here means
            // a reader does not have to know which methods are stochastic.
            "seed-affects-endmembers": args.method != "spa",
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
            // Per-subject abundance solves that ran out of iterations.
            // Zero for `--abundance ucls`, which is closed-form.
            extras.insert(
                "abundance_n_not_converged".into(),
                serde_json::json!(result.n_abundance_not_converged),
            );
            extras.insert(
                "k_selection_bound_hit".into(),
                serde_json::json!(k_selection_bound_hit(
                    args.k != "auto",
                    effective_k,
                    args.k_min,
                    args.k_max,
                )),
            );
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
    let mut reader = csv::ReaderBuilder::new()
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
