//! Cluster proteins by absence pattern across samples — the protein-side
//! dual of `recover-plex`. Two proteins that fail to quantify in the same
//! set of samples (e.g. share a peptide-level fragility, or always drop out
//! of the same TMT plexes) end up in the same cluster.
//!
//! Proteins with negligible absence (default ≤ `--min-absent`) are not
//! clustered; they are reported as a single "complete" group, since they
//! carry no signal for protein-side topology. Like `recover-plex`, the
//! quality summary fails honestly when the absence pattern is not
//! protein-coherent (separation_ratio < 2× → `weak_signal`, frac_singleton
//! > 0.5 → `absence_not_protein_coherent`).

use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;
use serde_json::json;
use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;
use std::time::SystemTime;

use crate::io::{
    atomic_write, hash_canonical_inputs, read_measurements_long, read_proteins, read_samples,
    sidecar_path_for, write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing canonical Atman TSV files.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output TSV path for per-protein cluster assignments.
    #[arg(long)]
    output: PathBuf,

    /// Jaccard distance below which two proteins are linked. Default 0.05 —
    /// strict, expects two proteins in the same cluster to share nearly all
    /// samples in which they're absent.
    #[arg(long, default_value_t = 0.05)]
    threshold: f64,

    /// Minimum number of absent samples a protein must have to be considered
    /// for clustering. Proteins with fewer absences are reported as the single
    /// "complete" cluster.
    #[arg(long, default_value_t = 5)]
    min_absent: usize,
}

#[derive(Debug, Clone)]
struct PerProtein {
    assay_id: String,
    gene_symbol: String,
    bitset: Vec<u64>,
    n_absent: usize,
}

#[derive(Debug)]
struct ClusterRow {
    assay_id: String,
    gene_symbol: String,
    cluster_id: String,
    cluster_size: usize,
    n_absent: usize,
    mean_within_jaccard_distance: f64,
}

#[derive(Debug)]
struct QualityRow {
    n_proteins: usize,
    n_samples_universe: usize,
    n_complete: usize,
    n_clusterable: usize,
    n_clusters: usize,
    n_singletons: usize,
    median_cluster_size: f64,
    mean_within_jaccard_distance: f64,
    mean_between_jaccard_distance: f64,
    separation_ratio: f64,
    threshold: f64,
    min_absent: usize,
    diagnostic: String,
}

pub fn run(args: Args) -> Result<()> {
    if !(0.0..=1.0).contains(&args.threshold) {
        bail!("--threshold must be in [0, 1]");
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

    let n_words = n_samples.div_ceil(64);
    // observed[p] has bit s set if (sample s, protein p) was observed.
    let mut observed: Vec<Vec<u64>> = vec![vec![0u64; n_words]; n_proteins];
    for m in &measurements {
        let Some(&si) = sample_index.get(&m.sample_id) else {
            continue;
        };
        let Some(&pi) = assay_index.get(&protein_key(m.platform.as_str(), &m.assay_id.0)) else {
            continue;
        };
        let w = si / 64;
        let b = si % 64;
        observed[pi][w] |= 1u64 << b;
    }

    // absence[p] = !observed[p], masked to n_samples-sized universe.
    let trailing = n_samples % 64;
    let last_mask: u64 = if trailing == 0 {
        !0u64
    } else {
        (1u64 << trailing) - 1
    };
    let mut per_protein: Vec<PerProtein> = Vec::with_capacity(n_proteins);
    for (pi, p) in proteins.iter().enumerate() {
        let mut bitset = vec![0u64; n_words];
        let mut n_absent = 0usize;
        for w in 0..n_words {
            let mut bits = !observed[pi][w];
            if w == n_words - 1 {
                bits &= last_mask;
            }
            bitset[w] = bits;
            n_absent += bits.count_ones() as usize;
        }
        per_protein.push(PerProtein {
            assay_id: p.assay_id.0.clone(),
            gene_symbol: p.gene_symbol.clone().unwrap_or_default(),
            bitset,
            n_absent,
        });
    }

    // Partition into "complete" (n_absent < min_absent) and "clusterable".
    let mut clusterable_idx: Vec<usize> = Vec::new();
    let mut complete_idx: Vec<usize> = Vec::new();
    for (i, p) in per_protein.iter().enumerate() {
        if p.n_absent < args.min_absent {
            complete_idx.push(i);
        } else {
            clusterable_idx.push(i);
        }
    }

    // Connected-components clustering on clusterable proteins.
    let m = clusterable_idx.len();
    let mut parent: Vec<usize> = (0..m).collect();
    fn find(parent: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while parent[r] != r {
            r = parent[r];
        }
        let mut cur = x;
        while parent[cur] != r {
            let nxt = parent[cur];
            parent[cur] = r;
            cur = nxt;
        }
        r
    }
    fn union(parent: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            if ra < rb {
                parent[rb] = ra;
            } else {
                parent[ra] = rb;
            }
        }
    }

    for i in 0..m {
        for j in (i + 1)..m {
            let pi = clusterable_idx[i];
            let pj = clusterable_idx[j];
            let d = jaccard_distance(&per_protein[pi].bitset, &per_protein[pj].bitset);
            if d <= args.threshold {
                union(&mut parent, i, j);
            }
        }
    }

    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..m {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(clusterable_idx[i]);
    }

    // Deterministic cluster ids: sort groups by min(assay_id) within group.
    let mut group_keys: Vec<(String, Vec<usize>)> = groups
        .into_values()
        .map(|members| {
            let key = members
                .iter()
                .map(|&pi| per_protein[pi].assay_id.clone())
                .min()
                .unwrap();
            (key, members)
        })
        .collect();
    group_keys.sort_by(|a, b| a.0.cmp(&b.0));

    // Per-cluster mean within-Jaccard distance.
    let mut all_within: Vec<f64> = Vec::new();
    let mut protein_to_cluster: HashMap<usize, (String, usize, f64)> = HashMap::new();
    for (idx, (_key, members)) in group_keys.iter().enumerate() {
        let cid = format!("absclust_{:05}", idx + 1);
        let size = members.len();
        let mean_within = if size <= 1 {
            f64::NAN
        } else {
            let mut sum = 0.0f64;
            let mut k = 0usize;
            for a in 0..members.len() {
                for b in (a + 1)..members.len() {
                    let d = jaccard_distance(
                        &per_protein[members[a]].bitset,
                        &per_protein[members[b]].bitset,
                    );
                    sum += d;
                    all_within.push(d);
                    k += 1;
                }
            }
            sum / k as f64
        };
        for &pi in members {
            protein_to_cluster.insert(pi, (cid.clone(), size, mean_within));
        }
    }

    // Mean between-cluster Jaccard distance.
    let mut between_sum = 0.0f64;
    let mut between_n = 0usize;
    for i in 0..m {
        for j in (i + 1)..m {
            let pi = clusterable_idx[i];
            let pj = clusterable_idx[j];
            if find(&mut parent, i) != find(&mut parent, j) {
                between_sum += jaccard_distance(&per_protein[pi].bitset, &per_protein[pj].bitset);
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

    // Build per-protein output rows.
    let mut rows: Vec<ClusterRow> = Vec::with_capacity(n_proteins);
    for (i, p) in per_protein.iter().enumerate() {
        if p.n_absent < args.min_absent {
            rows.push(ClusterRow {
                assay_id: p.assay_id.clone(),
                gene_symbol: p.gene_symbol.clone(),
                cluster_id: "complete".to_string(),
                cluster_size: complete_idx.len(),
                n_absent: p.n_absent,
                mean_within_jaccard_distance: f64::NAN,
            });
            continue;
        }
        let (cid, size, mean_within) = protein_to_cluster
            .get(&i)
            .cloned()
            .unwrap_or_else(|| ("singleton".to_string(), 1, f64::NAN));
        rows.push(ClusterRow {
            assay_id: p.assay_id.clone(),
            gene_symbol: p.gene_symbol.clone(),
            cluster_id: cid,
            cluster_size: size,
            n_absent: p.n_absent,
            mean_within_jaccard_distance: mean_within,
        });
    }
    write_rows(&args.output, &rows)?;

    let cluster_sizes: Vec<usize> = group_keys.iter().map(|(_, mm)| mm.len()).collect();
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
        clusterable_idx.len(),
        mean_within_global,
        mean_between,
        separation,
    );

    let quality = QualityRow {
        n_proteins,
        n_samples_universe: n_samples,
        n_complete: complete_idx.len(),
        n_clusterable: clusterable_idx.len(),
        n_clusters,
        n_singletons,
        median_cluster_size: median_size,
        mean_within_jaccard_distance: mean_within_global,
        mean_between_jaccard_distance: mean_between,
        separation_ratio: separation,
        threshold: args.threshold,
        min_absent: args.min_absent,
        diagnostic: diagnostic.clone(),
    };
    let quality_path = quality_path_for(&args.output);
    write_quality(&quality_path, &quality)?;

    if diagnostic != "protein_coherent" {
        eprintln!(
            "absence_topology: WARNING diagnostic={} separation={:.2} mean_within={:.4} mean_between={:.4}",
            diagnostic, separation, mean_within_global, mean_between
        );
    }
    eprintln!(
        "absence_topology: proteins={} complete={} clusterable={} clusters={} singletons={} median_size={} diagnostic={}",
        n_proteins,
        complete_idx.len(),
        clusterable_idx.len(),
        n_clusters,
        n_singletons,
        median_size as usize,
        diagnostic
    );

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_canonical_inputs(
        &args.input_dir,
        &["measurements.tsv", "samples.tsv", "proteins.tsv"],
    )?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "absence-topology",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output": args.output.display().to_string(),
            "quality-output": quality_path.display().to_string(),
            "threshold": args.threshold,
            "min-absent": args.min_absent,
        }),
        &inputs_sha256,
        &[args.output.clone(), quality_path.clone()],
        started_at,
        finished_at,
        None,
    )?;

    Ok(())
}

fn protein_key(platform: &str, assay_id: &str) -> (String, String) {
    (platform.to_string(), assay_id.to_string())
}

fn jaccard_distance(a: &[u64], b: &[u64]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    let mut inter = 0u64;
    let mut uni = 0u64;
    for (x, y) in a.iter().zip(b.iter()) {
        inter += (x & y).count_ones() as u64;
        uni += (x | y).count_ones() as u64;
    }
    if uni == 0 {
        0.0
    } else {
        1.0 - (inter as f64) / (uni as f64)
    }
}

fn median_f(v: &[f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        0.5 * (s[n / 2 - 1] + s[n / 2])
    }
}

fn classify(
    n_singletons: usize,
    n_clusterable: usize,
    mean_within: f64,
    mean_between: f64,
    separation: f64,
) -> String {
    if n_clusterable == 0 {
        return "no_clusterable_proteins".to_string();
    }
    if mean_between.is_nan() {
        return "single_cluster".to_string();
    }
    let frac_singleton = n_singletons as f64 / n_clusterable as f64;
    if frac_singleton > 0.5 {
        return "absence_not_protein_coherent".to_string();
    }
    if !mean_within.is_finite() || separation < 2.0 {
        return "weak_signal".to_string();
    }
    "protein_coherent".to_string()
}

fn write_rows(path: &std::path::Path, rows: &[ClusterRow]) -> Result<()> {
    let mut buf = String::from(
        "assay_id\tgene_symbol\tcluster_id\tcluster_size\tn_absent\tmean_within_jaccard_distance\n",
    );
    for r in rows {
        let mw = if r.mean_within_jaccard_distance.is_nan() {
            String::new()
        } else {
            format!("{:.6}", r.mean_within_jaccard_distance)
        };
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            r.assay_id, r.gene_symbol, r.cluster_id, r.cluster_size, r.n_absent, mw
        ));
    }
    atomic_write(path, buf.as_bytes()).with_context(|| format!("writing {:?}", path))
}

fn quality_path_for(out: &std::path::Path) -> PathBuf {
    let stem = out
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "absence_topology".to_string());
    let ext = out
        .extension()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "tsv".to_string());
    out.with_file_name(format!("{}_quality.{}", stem, ext))
}

fn write_quality(path: &std::path::Path, q: &QualityRow) -> Result<()> {
    let mut buf = String::from(
        "n_proteins\tn_samples_universe\tn_complete\tn_clusterable\tn_clusters\tn_singletons\tmedian_cluster_size\tmean_within_jaccard_distance\tmean_between_jaccard_distance\tseparation_ratio\tthreshold\tmin_absent\tdiagnostic\n",
    );
    buf.push_str(&format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        q.n_proteins,
        q.n_samples_universe,
        q.n_complete,
        q.n_clusterable,
        q.n_clusters,
        q.n_singletons,
        format_f(q.median_cluster_size),
        format_f(q.mean_within_jaccard_distance),
        format_f(q.mean_between_jaccard_distance),
        format_f(q.separation_ratio),
        format_f(q.threshold),
        q.min_absent,
        q.diagnostic,
    ));
    atomic_write(path, buf.as_bytes()).with_context(|| format!("writing {:?}", path))
}

fn format_f(x: f64) -> String {
    if x.is_nan() {
        String::new()
    } else if x.is_infinite() {
        if x > 0.0 {
            "inf".to_string()
        } else {
            "-inf".to_string()
        }
    } else {
        format!("{:.6}", x)
    }
}
