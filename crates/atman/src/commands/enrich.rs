use anyhow::{bail, Context, Result};
use atman_core::de::bh_fdr;
use clap::{Args as ClapArgs, Subcommand};
use csv::ReaderBuilder;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use super::enrich_gprofiler::{run_gprofiler, GprofilerArgs};
use crate::io::{atomic_write, need_col};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Over-representation analysis for significant DE genes.
    Ora(OraArgs),
    /// Live g:Profiler REST wrapper with on-disk cache.
    Gprofiler(GprofilerArgs),
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
        Command::Gprofiler(args) => run_gprofiler(args),
    }
}

fn run_ora(args: OraArgs) -> Result<()> {
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
            row.odds_ratio,
            row.p_value,
            row.bh_q.map(|q| q.to_string()).unwrap_or_default(),
            row.overlap_genes,
        ));
    }
    atomic_write(path, out.as_bytes())
}
