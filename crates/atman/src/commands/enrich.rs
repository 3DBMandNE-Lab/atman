use anyhow::{bail, Context, Result};
use atman_core::de::bh_fdr;
use atman_core::gsea::{gsea, GseaConfig};
use clap::{Args as ClapArgs, Subcommand, ValueEnum};
use csv::ReaderBuilder;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use super::enrich_gprofiler::{run_gprofiler, GprofilerArgs};
use crate::io::{
    atomic_write, format_f64, hash_labeled_inputs, need_col, sidecar_path_for, write_run_sidecar,
};
use serde_json::json;
use std::path::Path as StdPath;
use std::time::SystemTime;

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Over-representation analysis for significant DE genes.
    Ora(OraArgs),
    /// Pre-ranked gene-set enrichment analysis (Subramanian / fgseaSimple).
    Gsea(GseaCliArgs),
    /// Live g:Profiler REST wrapper with on-disk cache.
    Gprofiler(GprofilerArgs),
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum RankBy {
    /// `sign(mean_diff) * -log10(p_value)` — combines effect direction and significance.
    SignedLog10P,
    /// Raw t-statistic (signed by direction).
    T,
    /// Raw `mean_diff` (log2 fold-change).
    Log2Fc,
}

impl RankBy {
    fn key(&self) -> &'static str {
        match self {
            RankBy::SignedLog10P => "signed_log10_p",
            RankBy::T => "t",
            RankBy::Log2Fc => "log2_fc",
        }
    }
}

#[derive(ClapArgs, Debug)]
pub struct OraArgs {
    /// Path to de_results.tsv produced by `atman de`.
    #[arg(long)]
    de_results: PathBuf,

    /// Gene-set TSV with columns `set_name` and `gene_symbol`.
    #[arg(long)]
    gene_sets: PathBuf,

    /// Output TSV path.
    #[arg(long)]
    output: PathBuf,

    /// Optional comparison to analyze. If omitted, all comparisons are pooled.
    #[arg(long)]
    comparison: Option<String>,

    /// BH-q threshold for defining DE hits.
    #[arg(long, default_value_t = 0.05)]
    q: f64,

    /// Optional measured-universe file. Accepts a `gene_symbol` column or one gene per line.
    #[arg(long)]
    universe: Option<PathBuf>,
}

#[derive(Debug)]
struct OraRow {
    set_name: String,
    universe_size: usize,
    hit_count: usize,
    set_size: usize,
    overlap_size: usize,
    odds_ratio: f64,
    p_value: f64,
    bh_q: Option<f64>,
    overlap_genes: String,
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Ora(args) => run_ora(args),
        Command::Gsea(args) => run_gsea(args),
        Command::Gprofiler(args) => run_gprofiler(args),
    }
}

#[derive(ClapArgs, Debug)]
pub struct GseaCliArgs {
    /// Path to de_results.tsv produced by `atman de`.
    #[arg(long)]
    de_results: PathBuf,

    /// Gene-set TSV with columns `set_name` and `gene_symbol`.
    #[arg(long)]
    gene_sets: PathBuf,

    /// Output TSV path.
    #[arg(long)]
    output: PathBuf,

    /// Optional comparison to analyze. Required when de_results contains more than one.
    #[arg(long)]
    comparison: Option<String>,

    /// Ranking statistic.
    #[arg(long, value_enum, default_value_t = RankBy::SignedLog10P)]
    rank_by: RankBy,

    /// Number of gene-membership permutations for the null.
    #[arg(long, default_value_t = 1000)]
    n_permutations: usize,

    /// PRNG seed for the permutation null.
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Weighting exponent applied to |stat| during the running-sum walk.
    #[arg(long, default_value_t = 1.0)]
    weight_p: f64,

    /// Minimum set size after intersection with the ranked universe.
    #[arg(long, default_value_t = 5)]
    min_set_size: usize,

    /// Maximum set size after intersection with the ranked universe.
    #[arg(long, default_value_t = 500)]
    max_set_size: usize,
}

#[derive(Debug)]
struct GseaWriteRow {
    set_name: String,
    set_size: usize,
    es: f64,
    nes: f64,
    p_value: f64,
    bh_q: Option<f64>,
    leading_edge: String,
}

fn run_gsea(args: GseaCliArgs) -> Result<()> {
    let started_at = SystemTime::now();
    if args.n_permutations == 0 {
        bail!("--n-permutations must be > 0");
    }
    if !args.weight_p.is_finite() || args.weight_p < 0.0 {
        bail!("--weight-p must be finite and >= 0");
    }
    if args.min_set_size == 0 {
        bail!("--min-set-size must be >= 1");
    }
    if args.max_set_size < args.min_set_size {
        bail!(
            "--max-set-size ({}) must be >= --min-set-size ({})",
            args.max_set_size,
            args.min_set_size
        );
    }
    let ranked = read_ranked_genes(&args.de_results, args.comparison.as_deref(), args.rank_by)?;
    if ranked.is_empty() {
        bail!(
            "no ranked genes parsed from {:?} (empty DE results, all NaN, or filter excluded all rows)",
            args.de_results
        );
    }
    let gene_sets = read_gene_sets(&args.gene_sets)?;
    if gene_sets.is_empty() {
        bail!("no gene sets found in {:?}", args.gene_sets);
    }

    let cfg = GseaConfig {
        n_permutations: args.n_permutations,
        seed: args.seed,
        weight_p: args.weight_p,
        min_set_size: args.min_set_size,
        max_set_size: args.max_set_size,
    };
    let raw = gsea(&ranked, &gene_sets, &cfg);
    let p_values: Vec<Option<f64>> = raw.iter().map(|r| Some(r.p_value)).collect();
    let qs = bh_fdr(&p_values);
    let mut rows: Vec<GseaWriteRow> = raw
        .into_iter()
        .zip(qs)
        .map(|(r, q)| GseaWriteRow {
            set_name: r.set_name,
            set_size: r.set_size,
            es: r.es,
            nes: r.nes,
            p_value: r.p_value,
            bh_q: q,
            leading_edge: r.leading_edge.join(","),
        })
        .collect();
    rows.sort_by(|a, b| {
        match (a.bh_q, b.bh_q) {
            (Some(qa), Some(qb)) => qa.partial_cmp(&qb).unwrap_or(std::cmp::Ordering::Equal),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => std::cmp::Ordering::Equal,
        }
        .then_with(|| {
            b.nes
                .abs()
                .partial_cmp(&a.nes.abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .then_with(|| a.set_name.cmp(&b.set_name))
    });

    write_gsea_rows(&args.output, &rows)?;
    eprintln!(
        "enrich gsea: sets={} ranked_genes={} perms={} rank_by={}",
        rows.len(),
        ranked.len(),
        cfg.n_permutations,
        args.rank_by.key()
    );

    let finished_at = SystemTime::now();
    let labeled: Vec<(&str, &StdPath)> = vec![
        ("de_results", args.de_results.as_path()),
        ("gene_sets", args.gene_sets.as_path()),
    ];
    let inputs_sha256 = hash_labeled_inputs(&labeled)?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "enrich gsea",
        json!({
            "de-results": args.de_results.display().to_string(),
            "gene-sets": args.gene_sets.display().to_string(),
            "output": args.output.display().to_string(),
            "comparison": args.comparison,
            "rank-by": args.rank_by.key(),
            "n-permutations": args.n_permutations,
            "seed": args.seed,
            "weight-p": args.weight_p,
            "min-set-size": args.min_set_size,
            "max-set-size": args.max_set_size,
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("enrich gsea: sidecar={}", sidecar.display());
    Ok(())
}

/// Read DE results into a ranked gene list. Filters out rows with empty
/// gene symbols, non-finite ranking statistics, or non-empty `skip_reason`.
/// When multiple rows resolve to the same gene symbol (multiple panels),
/// the row whose absolute statistic is largest wins. Returns the resulting
/// list sorted descending by statistic, with stable secondary order by
/// gene symbol so ties in the statistic are deterministic.
fn read_ranked_genes(
    path: &Path,
    comparison: Option<&str>,
    rank_by: RankBy,
) -> Result<Vec<(String, f64)>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let gene_col = need_col(&headers, "gene_symbol", path)?;
    let comparison_col = need_col(&headers, "comparison", path)?;
    let skip_col = headers.iter().position(|h| h == "skip_reason");
    let stat_col = match rank_by {
        RankBy::T => need_col(&headers, "t", path)?,
        RankBy::Log2Fc => need_col(&headers, "mean_diff", path)?,
        RankBy::SignedLog10P => need_col(&headers, "p_value", path)?,
    };
    let direction_col = if matches!(rank_by, RankBy::SignedLog10P) {
        Some(need_col(&headers, "mean_diff", path)?)
    } else {
        None
    };

    let mut observed_comparisons: BTreeSet<String> = BTreeSet::new();
    let mut by_gene: BTreeMap<String, f64> = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let comp = row[comparison_col].trim().to_string();
        observed_comparisons.insert(comp.clone());
        if comparison.map(|wanted| comp != wanted).unwrap_or(false) {
            continue;
        }
        if let Some(idx) = skip_col {
            if !row[idx].trim().is_empty() {
                continue;
            }
        }
        let gene = row[gene_col].trim().to_string();
        if gene.is_empty() {
            continue;
        }
        let stat: f64 = match rank_by {
            RankBy::T | RankBy::Log2Fc => match row[stat_col].trim().parse() {
                Ok(v) => v,
                Err(_) => continue,
            },
            RankBy::SignedLog10P => {
                let p: f64 = match row[stat_col].trim().parse() {
                    Ok(v) => v,
                    Err(_) => continue,
                };
                let dir: f64 = match row[direction_col.unwrap()].trim().parse() {
                    Ok(v) => v,
                    Err(_) => continue,
                };
                if !(p.is_finite() && (0.0..=1.0).contains(&p) && dir.is_finite()) {
                    continue;
                }
                let p_floor = p.max(1e-300);
                dir.signum() * -p_floor.log10()
            }
        };
        if !stat.is_finite() {
            continue;
        }
        let entry = by_gene.entry(gene).or_insert(stat);
        if stat.abs() > entry.abs() {
            *entry = stat;
        }
    }
    if comparison.is_none() && observed_comparisons.len() > 1 {
        bail!(
            "{:?} contains {} distinct comparisons ({:?}); pass --comparison to choose one",
            path,
            observed_comparisons.len(),
            observed_comparisons
        );
    }
    let mut out: Vec<(String, f64)> = by_gene.into_iter().collect();
    out.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.0.cmp(&b.0))
    });
    Ok(out)
}

fn write_gsea_rows(path: &Path, rows: &[GseaWriteRow]) -> Result<()> {
    let mut out = String::from("set_name\tset_size\tes\tnes\tp_value\tbh_q\tleading_edge\n");
    for row in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            row.set_name,
            row.set_size,
            format_f64(row.es),
            format_f64(row.nes),
            format_f64(row.p_value),
            row.bh_q.map(format_f64).unwrap_or_default(),
            row.leading_edge,
        ));
    }
    atomic_write(path, out.as_bytes())
}

fn run_ora(args: OraArgs) -> Result<()> {
    let started_at = SystemTime::now();
    if !(args.q.is_finite() && (0.0..=1.0).contains(&args.q)) {
        bail!("--q must satisfy 0 <= q <= 1");
    }
    let de = read_de_genes(&args.de_results, args.comparison.as_deref(), args.q)?;
    if de.universe.is_empty() {
        bail!("no measured genes found in {:?}", args.de_results);
    }
    let universe = if let Some(path) = &args.universe {
        read_universe(path)?
    } else {
        de.universe
    };
    let hits: BTreeSet<String> = de.hits.intersection(&universe).cloned().collect();
    let gene_sets = read_gene_sets(&args.gene_sets)?;
    if gene_sets.is_empty() {
        bail!("no gene sets found in {:?}", args.gene_sets);
    }

    let mut rows = Vec::new();
    let mut p_values = Vec::new();
    for (set_name, genes) in gene_sets {
        let set_genes: BTreeSet<String> = genes.intersection(&universe).cloned().collect();
        if set_genes.is_empty() {
            continue;
        }
        let overlap: BTreeSet<String> = set_genes.intersection(&hits).cloned().collect();
        let a = overlap.len();
        let b = hits.len().saturating_sub(a);
        let c = set_genes.len().saturating_sub(a);
        let d = universe.len().saturating_sub(a + b + c);
        let p = fisher_upper_tail(a, hits.len(), set_genes.len(), universe.len());
        p_values.push(Some(p));
        rows.push(OraRow {
            set_name,
            universe_size: universe.len(),
            hit_count: hits.len(),
            set_size: set_genes.len(),
            overlap_size: a,
            odds_ratio: odds_ratio(a, b, c, d),
            p_value: p,
            bh_q: None,
            overlap_genes: overlap.into_iter().collect::<Vec<_>>().join(","),
        });
    }
    let qs = bh_fdr(&p_values);
    for (row, q) in rows.iter_mut().zip(qs) {
        row.bh_q = q;
    }
    rows.sort_by(|a, b| {
        match (a.bh_q, b.bh_q) {
            (Some(qa), Some(qb)) => qa.partial_cmp(&qb).unwrap_or(std::cmp::Ordering::Equal),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => std::cmp::Ordering::Equal,
        }
        .then_with(|| b.overlap_size.cmp(&a.overlap_size))
        .then_with(|| a.set_name.cmp(&b.set_name))
    });
    write_rows(&args.output, &rows)?;
    eprintln!(
        "enrich ora: sets={} hits={} universe={} q={}",
        rows.len(),
        hits.len(),
        universe.len(),
        args.q
    );

    let finished_at = SystemTime::now();
    let mut labeled: Vec<(&str, &StdPath)> = vec![
        ("de_results", args.de_results.as_path()),
        ("gene_sets", args.gene_sets.as_path()),
    ];
    if let Some(p) = &args.universe {
        labeled.push(("universe", p.as_path()));
    }
    let inputs_sha256 = hash_labeled_inputs(&labeled)?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "enrich ora",
        json!({
            "de-results": args.de_results.display().to_string(),
            "gene-sets": args.gene_sets.display().to_string(),
            "output": args.output.display().to_string(),
            "comparison": args.comparison,
            "q": args.q,
            "universe": args.universe.as_ref().map(|p| p.display().to_string()),
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("enrich ora: sidecar={}", sidecar.display());
    Ok(())
}

struct DeGenes {
    universe: BTreeSet<String>,
    hits: BTreeSet<String>,
}

fn read_de_genes(path: &Path, comparison: Option<&str>, q_threshold: f64) -> Result<DeGenes> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let gene_col = need_col(&headers, "gene_symbol", path)?;
    let comparison_col = need_col(&headers, "comparison", path)?;
    let q_col = need_col(&headers, "bh_q", path)?;
    let mut universe = BTreeSet::new();
    let mut hits = BTreeSet::new();
    for row in reader.records() {
        let row = row?;
        if comparison
            .map(|wanted| row[comparison_col].trim() != wanted)
            .unwrap_or(false)
        {
            continue;
        }
        let gene = row[gene_col].trim();
        if gene.is_empty() {
            continue;
        }
        universe.insert(gene.to_string());
        if let Ok(q) = row[q_col].trim().parse::<f64>() {
            if q.is_finite() && q <= q_threshold {
                hits.insert(gene.to_string());
            }
        }
    }
    Ok(DeGenes { universe, hits })
}

fn read_gene_sets(path: &Path) -> Result<BTreeMap<String, BTreeSet<String>>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let set_col = need_col(&headers, "set_name", path)?;
    let gene_col = need_col(&headers, "gene_symbol", path)?;
    let mut out: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let set = row[set_col].trim();
        let gene = row[gene_col].trim();
        if set.is_empty() || gene.is_empty() {
            continue;
        }
        out.entry(set.to_string())
            .or_default()
            .insert(gene.to_string());
    }
    Ok(out)
}

fn read_universe(path: &Path) -> Result<BTreeSet<String>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {:?}", path))?;
    let first = text.lines().next().unwrap_or("");
    if first.split('\t').any(|h| h == "gene_symbol") {
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_reader(text.as_bytes());
        let headers = reader.headers()?.clone();
        let gene_col = need_col(&headers, "gene_symbol", path)?;
        let mut out = BTreeSet::new();
        for row in reader.records() {
            let row = row?;
            let gene = row[gene_col].trim();
            if !gene.is_empty() {
                out.insert(gene.to_string());
            }
        }
        Ok(out)
    } else {
        Ok(text
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .map(str::to_string)
            .collect())
    }
}

fn fisher_upper_tail(overlap: usize, hits: usize, set_size: usize, universe: usize) -> f64 {
    let max_k = hits.min(set_size);
    (overlap..=max_k)
        .map(|k| hypergeom_pmf(k, hits, set_size, universe))
        .sum::<f64>()
        .min(1.0)
}

fn hypergeom_pmf(k: usize, hits: usize, set_size: usize, universe: usize) -> f64 {
    if k > hits || k > set_size || set_size - k > universe - hits {
        return 0.0;
    }
    (log_choose(hits, k) + log_choose(universe - hits, set_size - k)
        - log_choose(universe, set_size))
    .exp()
}

fn log_choose(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    let k = k.min(n - k);
    (1..=k)
        .map(|i| ((n + 1 - i) as f64).ln() - (i as f64).ln())
        .sum()
}

fn odds_ratio(a: usize, b: usize, c: usize, d: usize) -> f64 {
    if b == 0 || c == 0 {
        if a == 0 {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        (a as f64 * d as f64) / (b as f64 * c as f64)
    }
}

fn write_rows(path: &Path, rows: &[OraRow]) -> Result<()> {
    let mut out = String::from(
        "set_name\tuniverse_size\thit_count\tset_size\toverlap_size\todds_ratio\tp_value\tbh_q\toverlap_genes\n",
    );
    for row in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            row.set_name,
            row.universe_size,
            row.hit_count,
            row.set_size,
            row.overlap_size,
            format_f64(row.odds_ratio),
            format_f64(row.p_value),
            row.bh_q.map(format_f64).unwrap_or_default(),
            row.overlap_genes,
        ));
    }
    atomic_write(path, out.as_bytes())
}
