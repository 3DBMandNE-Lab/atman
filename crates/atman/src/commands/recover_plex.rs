//! Recover TMT plex assignment from the absence pattern in canonical
//! `measurements.tsv`. Useful when ingest dropped `plate_id` but plex-level
//! missingness is plex-categorical (a protein is absent for every sample in
//! its plex). Validates honestly: if absence is not plex-coherent, the
//! quality summary's `diagnostic` column says so and downstream uses must
//! treat the recovered ids as unreliable.

use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;
use serde_json::json;
use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;
use std::time::SystemTime;

use crate::commands::jaccard_cluster::{
    find, format_f, jaccard_distance, median_f, quality_path_for, union,
};
use crate::io::{
    atomic_write, hash_canonical_inputs, read_measurements_long, read_proteins, read_samples,
    sidecar_path_for, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output TSV path for per-sample plex assignments.
    #[arg(long)]
    output: PathBuf,

    /// Jaccard distance below which two samples are linked into the same component.
    /// Default 0.05 — strict, expects within-plex sharing of nearly all absent proteins.
    #[arg(long, default_value_t = 0.05)]
    threshold: f64,

    /// Expected plex multiplicity (TMT-11 typically 9–11 samples per plex). Used only
    /// for diagnostic warnings on cluster size.
    #[arg(long, default_value_t = 11)]
    expected_cluster_size: usize,

    /// If set, write samples.tsv augmented with a `recovered_plex_id` column to this
    /// path. Existing DE/detectability commands can then reference it via `--design`.
    /// Only emitted when diagnostic == plex_coherent, to avoid downstream poisoning.
    #[arg(long)]
    augmented_samples_output: Option<PathBuf>,

    /// Infer patient-level pairing within each plex: matches each sample of condition A
    /// to a sample of condition B by minimal sample_id-stem distance, accepting only
    /// stem-distance ≤ `--pair-stem-distance`. Format: "A-B" (e.g. "tumor-paired_non_tumor").
    /// When inference is unique and complete, a `patient_id` column is added to the
    /// augmented samples output (if requested) and a companion pairs TSV
    /// (`<output-stem>_pairs.tsv`) is written. Requires numeric sample-id stems.
    #[arg(long)]
    infer_pairs: Option<String>,

    /// Maximum allowed |numeric stem difference| between paired samples in --infer-pairs.
    #[arg(long, default_value_t = 1)]
    pair_stem_distance: u64,
}

#[derive(Debug, Clone)]
struct PerSample {
    sample_id: String,
    bitset: Vec<u64>,
}

#[derive(Debug)]
struct ClusterRow {
    sample_id: String,
    plex_id: usize,
    cluster_size: usize,
    mean_within_jaccard_distance: f64,
}

#[derive(Debug)]
struct QualityRow {
    n_samples: usize,
    n_proteins_universe: usize,
    n_clusters: usize,
    n_singletons: usize,
    median_cluster_size: f64,
    mean_within_jaccard_distance: f64,
    mean_between_jaccard_distance: f64,
    separation_ratio: f64,
    expected_cluster_size: usize,
    threshold: f64,
    diagnostic: String,
}

pub fn run(args: Args) -> Result<()> {
    if !(0.0..=1.0).contains(&args.threshold) {
        bail!("--threshold must be in [0, 1]");
    }
    if args.expected_cluster_size == 0 {
        bail!("--expected-cluster-size must be > 0");
    }
    let started_at = SystemTime::now();

    let samples = read_samples(&args.input_dir.join("samples.tsv"))?;
    let proteins = read_proteins(&args.input_dir.join("proteins.tsv"))?;
    let measurements = read_measurements_long(&args.input_dir.join("measurements.tsv"))?;

    let mut sample_index: HashMap<String, usize> = HashMap::new();
    for (i, s) in samples.iter().enumerate() {
        sample_index.insert(s.sample_id.clone(), i);
    }

    let mut assay_index: HashMap<(String, String), usize> = HashMap::new();
    for (i, p) in proteins.iter().enumerate() {
        assay_index.insert(protein_key(p.platform.as_str(), &p.assay_id.0), i);
    }

    let n_samples = samples.len();
    let n_proteins = proteins.len();
    if n_samples == 0 {
        bail!("samples.tsv has no samples");
    }
    if n_proteins == 0 {
        bail!("proteins.tsv has no assays");
    }

    let n_words = n_proteins.div_ceil(64);
    // observed_bits[s] has bit p set if (sample s, protein p) was observed in measurements.tsv.
    let mut observed: Vec<Vec<u64>> = vec![vec![0u64; n_words]; n_samples];
    for m in &measurements {
        let Some(&si) = sample_index.get(&m.sample_id) else {
            continue;
        };
        let Some(&pi) = assay_index.get(&protein_key(m.platform.as_str(), &m.assay_id.0)) else {
            continue;
        };
        let w = pi / 64;
        let b = pi % 64;
        observed[si][w] |= 1u64 << b;
    }

    // absence[s] = !observed[s], masked to the n_proteins-sized universe.
    let mut per_sample: Vec<PerSample> = Vec::with_capacity(n_samples);
    let trailing = n_proteins % 64;
    let last_mask: u64 = if trailing == 0 {
        !0u64
    } else {
        (1u64 << trailing) - 1
    };
    for (si, s) in samples.iter().enumerate() {
        let mut bitset = vec![0u64; n_words];
        for w in 0..n_words {
            let mut bits = !observed[si][w];
            if w == n_words - 1 {
                bits &= last_mask;
            }
            bitset[w] = bits;
        }
        per_sample.push(PerSample {
            sample_id: s.sample_id.clone(),
            bitset,
        });
    }

    // Pairwise Jaccard distance, connected-component clustering.
    let mut parent: Vec<usize> = (0..n_samples).collect();

    let mut all_distances: Vec<f64> = Vec::with_capacity(n_samples * (n_samples - 1) / 2);
    let mut edges: Vec<(usize, usize, f64)> = Vec::new();
    for i in 0..n_samples {
        for j in (i + 1)..n_samples {
            let d = jaccard_distance(&per_sample[i].bitset, &per_sample[j].bitset);
            all_distances.push(d);
            if d <= args.threshold {
                edges.push((i, j, d));
                union(&mut parent, i, j);
            }
        }
    }

    // Group samples by component root.
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n_samples {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }

    // Assign deterministic plex ids by min(sample_id) within group.
    let mut group_keys: Vec<(String, Vec<usize>)> = groups
        .into_values()
        .map(|members| {
            let key = members
                .iter()
                .map(|&m| per_sample[m].sample_id.clone())
                .min()
                .unwrap();
            (key, members)
        })
        .collect();
    group_keys.sort_by(|a, b| a.0.cmp(&b.0));

    // Per-cluster mean within-Jaccard distance.
    let mut sample_to_plex: HashMap<usize, (usize, usize, f64)> = HashMap::new();
    let mut all_within: Vec<f64> = Vec::new();
    for (plex_idx, (_key, members)) in group_keys.iter().enumerate() {
        let plex_id = plex_idx + 1;
        let size = members.len();
        let mean_within = if size <= 1 {
            f64::NAN
        } else {
            let mut sum = 0.0f64;
            let mut k = 0usize;
            for a in 0..members.len() {
                for b in (a + 1)..members.len() {
                    let d = jaccard_distance(
                        &per_sample[members[a]].bitset,
                        &per_sample[members[b]].bitset,
                    );
                    sum += d;
                    all_within.push(d);
                    k += 1;
                }
            }
            sum / k as f64
        };
        for &m in members {
            sample_to_plex.insert(m, (plex_id, size, mean_within));
        }
    }

    // Build per-sample output rows in samples.tsv canonical order.
    let mut rows: Vec<ClusterRow> = Vec::with_capacity(n_samples);
    for (si, s) in samples.iter().enumerate() {
        let (plex_id, size, mean_within) = sample_to_plex[&si];
        rows.push(ClusterRow {
            sample_id: s.sample_id.clone(),
            plex_id,
            cluster_size: size,
            mean_within_jaccard_distance: mean_within,
        });
    }
    write_rows(&args.output, &rows)?;

    // Mean between-cluster distance: pairs whose samples are in different components.
    let mut between_sum = 0.0f64;
    let mut between_n = 0usize;
    for i in 0..n_samples {
        for j in (i + 1)..n_samples {
            if find(&mut parent, i) != find(&mut parent, j) {
                between_sum += jaccard_distance(&per_sample[i].bitset, &per_sample[j].bitset);
                between_n += 1;
            }
        }
    }
    let mean_between = if between_n == 0 {
        f64::NAN
    } else {
        between_sum / between_n as f64
    };
    let mean_within_global = if all_within.is_empty() {
        f64::NAN
    } else {
        all_within.iter().sum::<f64>() / all_within.len() as f64
    };

    let cluster_sizes: Vec<usize> = group_keys.iter().map(|(_, m)| m.len()).collect();
    let n_clusters = cluster_sizes.len();
    let n_singletons = cluster_sizes.iter().filter(|&&s| s == 1).count();
    let median_size = median_f(&cluster_sizes.iter().map(|&s| s as f64).collect::<Vec<_>>());
    let separation =
        if mean_within_global.is_nan() || mean_between.is_nan() || mean_within_global == 0.0 {
            if mean_between > 0.0 {
                f64::INFINITY
            } else {
                f64::NAN
            }
        } else {
            mean_between / mean_within_global
        };

    let diagnostic = classify(
        n_singletons,
        n_samples,
        median_size,
        mean_within_global,
        mean_between,
        separation,
        args.expected_cluster_size,
    );

    let quality = QualityRow {
        n_samples,
        n_proteins_universe: n_proteins,
        n_clusters,
        n_singletons,
        median_cluster_size: median_size,
        mean_within_jaccard_distance: mean_within_global,
        mean_between_jaccard_distance: mean_between,
        separation_ratio: separation,
        expected_cluster_size: args.expected_cluster_size,
        threshold: args.threshold,
        diagnostic: diagnostic.clone(),
    };
    let quality_path = quality_path_for(&args.output, "recovered_plex");
    write_quality(&quality_path, &quality)?;

    if diagnostic != "plex_coherent" {
        eprintln!(
            "recover_plex: WARNING diagnostic={} separation={:.2} median_within={:.4} median_between={:.4}",
            diagnostic, separation, mean_within_global, mean_between
        );
    }

    let mut optional_outputs: Vec<PathBuf> = Vec::new();
    let pair_assignments: Option<Vec<PairRow>> = if let Some(spec) = args.infer_pairs.as_ref() {
        if diagnostic == "plex_coherent" {
            // The plex is coherent and the user explicitly asked for pairs, so an
            // inference failure is a hard error under strict-failure: returning
            // exit 0 with no pairs file would silently hand back a "success" that
            // omits exactly what was requested. (The plex-incoherent SKIP below is
            // a legitimate, documented non-error path.)
            let pairs = infer_pairs(spec, args.pair_stem_distance, &samples, &rows)
                .with_context(|| format!("--infer-pairs {:?} failed", spec))?;
            let pairs_path = pairs_path_for(&args.output);
            write_pairs(&pairs_path, &pairs)?;
            optional_outputs.push(pairs_path.clone());
            eprintln!(
                "recover_plex: inferred {} patient pairs ({}); wrote {}",
                pairs.len(),
                spec,
                pairs_path.display()
            );
            Some(pairs)
        } else {
            eprintln!(
                "recover_plex: skipping --infer-pairs (diagnostic={}); plex assignment is unreliable",
                diagnostic
            );
            None
        }
    } else {
        None
    };

    if let Some(aug_path) = args.augmented_samples_output.as_ref() {
        if diagnostic == "plex_coherent" {
            write_augmented_samples(
                &args.input_dir.join("samples.tsv"),
                aug_path,
                &rows,
                pair_assignments.as_deref(),
            )?;
            optional_outputs.push(aug_path.clone());
            eprintln!(
                "recover_plex: wrote augmented samples to {}",
                aug_path.display()
            );
        } else {
            eprintln!(
                "recover_plex: skipping --augmented-samples-output (diagnostic={}); refusing to emit unreliable plex ids",
                diagnostic
            );
        }
    }
    eprintln!(
        "recover_plex: samples={} proteins={} clusters={} singletons={} median_size={} diagnostic={}",
        n_samples, n_proteins, n_clusters, n_singletons, median_size as usize, diagnostic
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &["measurements.tsv", "samples.tsv", "proteins.tsv"],
    )?;
    let sidecar = sidecar_path_for(&args.output);
    let mut outputs = vec![args.output.clone(), quality_path.clone()];
    outputs.extend(optional_outputs);
    write_run_sidecar(
        &sidecar,
        "recover-plex",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output": args.output.display().to_string(),
            "quality-output": quality_path.display().to_string(),
            "augmented-samples-output": args.augmented_samples_output.as_ref().map(|p| p.display().to_string()),
            "infer-pairs": args.infer_pairs,
            "pair-stem-distance": args.pair_stem_distance,
            "threshold": args.threshold,
            "expected-cluster-size": args.expected_cluster_size,
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        None,
    )?;

    let _ = edges; // edges retained only for potential future debug emit
    let _ = all_distances;

    Ok(())
}

fn protein_key(platform: &str, assay_id: &str) -> (String, String) {
    (platform.to_string(), assay_id.to_string())
}

fn classify(
    n_singletons: usize,
    n_samples: usize,
    median_size: f64,
    mean_within: f64,
    mean_between: f64,
    separation: f64,
    expected: usize,
) -> String {
    if mean_between.is_nan() {
        return "single_cluster".to_string();
    }
    let frac_singleton = n_singletons as f64 / n_samples as f64;
    if frac_singleton > 0.5 {
        return "absence_not_plex_coherent".to_string();
    }
    if !mean_within.is_finite() || separation < 2.0 {
        return "weak_plex_signal".to_string();
    }
    let lo = (expected as f64 * 0.5).max(2.0);
    let hi = expected as f64 * 2.0;
    if median_size < lo || median_size > hi {
        return "atypical_cluster_size".to_string();
    }
    "plex_coherent".to_string()
}

fn write_rows(path: &std::path::Path, rows: &[ClusterRow]) -> Result<()> {
    let mut buf =
        String::from("sample_id\trecovered_plex_id\tcluster_size\tmean_within_jaccard_distance\n");
    for r in rows {
        let mw = if r.mean_within_jaccard_distance.is_nan() {
            String::new()
        } else {
            format!("{:.6}", r.mean_within_jaccard_distance)
        };
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            r.sample_id, r.plex_id, r.cluster_size, mw
        ));
    }
    atomic_write(path, buf.as_bytes()).with_context(|| format!("writing {:?}", path))
}

fn write_quality(path: &std::path::Path, q: &QualityRow) -> Result<()> {
    let mut buf = String::new();
    buf.push_str(
        "n_samples\tn_proteins_universe\tn_clusters\tn_singletons\tmedian_cluster_size\tmean_within_jaccard_distance\tmean_between_jaccard_distance\tseparation_ratio\texpected_cluster_size\tthreshold\tdiagnostic\n",
    );
    buf.push_str(&format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        q.n_samples,
        q.n_proteins_universe,
        q.n_clusters,
        q.n_singletons,
        format_f(q.median_cluster_size),
        format_f(q.mean_within_jaccard_distance),
        format_f(q.mean_between_jaccard_distance),
        format_f(q.separation_ratio),
        q.expected_cluster_size,
        format_f(q.threshold),
        q.diagnostic,
    ));
    atomic_write(path, buf.as_bytes()).with_context(|| format!("writing {:?}", path))
}

fn write_augmented_samples(
    samples_path: &std::path::Path,
    out_path: &std::path::Path,
    rows: &[ClusterRow],
    pairs: Option<&[PairRow]>,
) -> Result<()> {
    let plex_for: HashMap<&str, usize> = rows
        .iter()
        .map(|r| (r.sample_id.as_str(), r.plex_id))
        .collect();
    let patient_for: HashMap<&str, &str> = pairs
        .map(|ps| {
            let mut m: HashMap<&str, &str> = HashMap::new();
            for p in ps {
                m.insert(p.sample_a.as_str(), p.patient_id.as_str());
                m.insert(p.sample_b.as_str(), p.patient_id.as_str());
            }
            m
        })
        .unwrap_or_default();
    let text = std::fs::read_to_string(samples_path)
        .with_context(|| format!("reading {:?}", samples_path))?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("samples.tsv has no header"))?;
    let mut buf = String::new();
    buf.push_str(header);
    buf.push_str("\trecovered_plex_id");
    if pairs.is_some() {
        buf.push_str("\tpatient_id");
    }
    buf.push('\n');
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let sample_id = line.split('\t').next().unwrap_or_default();
        let plex = plex_for
            .get(sample_id)
            .ok_or_else(|| anyhow::anyhow!("sample {} not in recovered plex table", sample_id))?;
        buf.push_str(line);
        // Pad to 4 digits so categorical encoding is lex-sorted = numeric-sorted.
        buf.push_str(&format!("\tplex_{:04}", plex));
        if pairs.is_some() {
            let pid = patient_for.get(sample_id).copied().unwrap_or("");
            buf.push('\t');
            buf.push_str(pid);
        }
        buf.push('\n');
    }
    atomic_write(out_path, buf.as_bytes()).with_context(|| format!("writing {:?}", out_path))
}

#[derive(Debug, Clone)]
struct PairRow {
    patient_id: String,
    sample_a: String,
    sample_b: String,
    plex_id: usize,
    stem_diff: u64,
}

fn parse_pair_spec(spec: &str) -> Result<(String, String)> {
    let mut parts = spec.split('-').map(str::trim);
    let a = parts.next().unwrap_or_default();
    let b = parts.next().unwrap_or_default();
    if a.is_empty() || b.is_empty() || parts.next().is_some() {
        bail!(
            "--infer-pairs must be `A-B` (e.g. `tumor-paired_non_tumor`); got {:?}",
            spec
        );
    }
    Ok((a.to_string(), b.to_string()))
}

fn numeric_stem(sample_id: &str) -> Option<u64> {
    let digits: String = sample_id
        .chars()
        .skip_while(|c| !c.is_ascii_digit())
        .take_while(|c| c.is_ascii_digit())
        .collect();
    digits.parse::<u64>().ok()
}

fn infer_pairs(
    spec: &str,
    max_stem_diff: u64,
    samples: &[atman_core::Sample],
    rows: &[ClusterRow],
) -> Result<Vec<PairRow>> {
    let (cond_a, cond_b) = parse_pair_spec(spec)?;
    let plex_for: HashMap<&str, usize> = rows
        .iter()
        .map(|r| (r.sample_id.as_str(), r.plex_id))
        .collect();
    let cond_for: HashMap<&str, &str> = samples
        .iter()
        .map(|s| {
            (
                s.sample_id.as_str(),
                s.condition.as_deref().unwrap_or_default(),
            )
        })
        .collect();

    // Group sample_ids in conditions A and B by plex.
    let mut by_plex_a: BTreeMap<usize, Vec<(u64, String)>> = BTreeMap::new();
    let mut by_plex_b: BTreeMap<usize, Vec<(u64, String)>> = BTreeMap::new();
    let mut a_total = 0usize;
    for s in samples {
        let Some(&plex) = plex_for.get(s.sample_id.as_str()) else {
            continue;
        };
        let Some(stem) = numeric_stem(&s.sample_id) else {
            bail!(
                "sample_id {:?} has no numeric stem; --infer-pairs requires numeric ids",
                s.sample_id
            );
        };
        match cond_for.get(s.sample_id.as_str()).copied().unwrap_or("") {
            c if c == cond_a => {
                by_plex_a
                    .entry(plex)
                    .or_default()
                    .push((stem, s.sample_id.clone()));
                a_total += 1;
            }
            c if c == cond_b => {
                by_plex_b
                    .entry(plex)
                    .or_default()
                    .push((stem, s.sample_id.clone()));
            }
            _ => {}
        }
    }

    let mut pairs: Vec<PairRow> = Vec::new();
    for (plex, a_list) in &by_plex_a {
        let mut a_sorted = a_list.clone();
        a_sorted.sort();
        let mut b_sorted = by_plex_b.get(plex).cloned().unwrap_or_default();
        b_sorted.sort();
        let mut used = vec![false; b_sorted.len()];
        for (a_stem, a_id) in &a_sorted {
            // Pick the B with smallest |stem-diff| that is still unused and within tolerance.
            let mut best: Option<(u64, usize)> = None;
            for (i, (b_stem, _)) in b_sorted.iter().enumerate() {
                if used[i] {
                    continue;
                }
                let diff = b_stem.abs_diff(*a_stem);
                if diff > max_stem_diff {
                    continue;
                }
                if best.is_none_or(|(d, _)| diff < d) {
                    best = Some((diff, i));
                }
            }
            let Some((diff, idx)) = best else {
                bail!(
                    "no eligible {} partner for {} within stem-distance {} in plex {}",
                    cond_b,
                    a_id,
                    max_stem_diff,
                    plex
                );
            };
            used[idx] = true;
            let (_b_stem, b_id) = &b_sorted[idx];
            pairs.push(PairRow {
                patient_id: format!("pair_{:04}", pairs.len() + 1),
                sample_a: a_id.clone(),
                sample_b: b_id.clone(),
                plex_id: *plex,
                stem_diff: diff,
            });
        }
    }
    if pairs.len() != a_total {
        bail!(
            "pair inference incomplete: paired {} of {} {} samples",
            pairs.len(),
            a_total,
            cond_a
        );
    }
    Ok(pairs)
}

fn pairs_path_for(out: &std::path::Path) -> PathBuf {
    let stem = out
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "recovered_plex".to_string());
    let ext = out
        .extension()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "tsv".to_string());
    out.with_file_name(format!("{}_pairs.{}", stem, ext))
}

fn write_pairs(path: &std::path::Path, pairs: &[PairRow]) -> Result<()> {
    let mut buf = String::from("patient_id\tsample_a\tsample_b\tplex_id\tstem_diff\n");
    for p in pairs {
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            p.patient_id, p.sample_a, p.sample_b, p.plex_id, p.stem_diff
        ));
    }
    atomic_write(path, buf.as_bytes()).with_context(|| format!("writing {:?}", path))
}
