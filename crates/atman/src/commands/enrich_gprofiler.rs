//! Live g:Profiler REST wrapper with on-disk query cache.
//!
//! Query payloads are canonicalized to JSON and hashed (SHA-256). The hash is
//! the cache key. Hits are read back from disk; misses POST to g:Profiler and
//! persist the raw JSON response before parsing. Pinning `--ontology-version`
//! into the cache key makes re-runs deterministic under ontology updates:
//! if you set `--ontology-version 2026-04-01`, any upstream ontology change
//! produces a different cache key (forcing a re-fetch) while the same
//! `--ontology-version` always resolves to the same cached bytes.

use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;
use csv::ReaderBuilder;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use crate::io::{atomic_write, escape_tsv, sha256_hex};

const DEFAULT_ENDPOINT: &str = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/";

#[derive(ClapArgs, Debug)]
pub struct GprofilerArgs {
    /// Optional DE-results TSV (produced by `atman de`). Yields a single query of
    /// genes whose `bh_q <= --user-threshold` for the selected comparison.
    #[arg(long, conflicts_with = "query_tsv")]
    pub de_results: Option<PathBuf>,

    /// Comparison to filter when reading `--de-results`.
    #[arg(long)]
    pub comparison: Option<String>,

    /// Optional per-query TSV with columns `query\tgene_symbol`. Each distinct
    /// query string becomes its own g:Profiler submission. Typical usage is
    /// `query = program_01 ...` for ICA program annotation.
    #[arg(long, conflicts_with = "de_results")]
    pub query_tsv: Option<PathBuf>,

    /// Background universe TSV (one gene per line or a `gene_symbol` column).
    /// Sent to g:Profiler as `background`.
    #[arg(long)]
    pub background: PathBuf,

    /// Organism tag (e.g. `hsapiens`).
    #[arg(long, default_value = "hsapiens")]
    pub organism: String,

    /// Comma-separated source list (e.g. `GO:BP,GO:MF,KEGG,REAC,WP`).
    #[arg(long, default_value = "GO:BP,GO:MF,GO:CC,KEGG,REAC,WP")]
    pub sources: String,

    /// Significance threshold method: `fdr`, `bonferroni`, or `g_SCS`.
    #[arg(long, default_value = "fdr")]
    pub threshold_method: String,

    /// User threshold on DE q-values when building the DE query, and on
    /// g:Profiler `user_threshold` for significance filtering in the response.
    #[arg(long, default_value_t = 0.05)]
    pub user_threshold: f64,

    /// Cache directory. Created on first run; persistent across runs.
    #[arg(long)]
    pub cache_dir: PathBuf,

    /// Ontology version string, baked into the cache key so an ontology bump
    /// forces a re-fetch while identical runs stay deterministic.
    #[arg(long, default_value = "pinned")]
    pub ontology_version: String,

    /// Disable live POST; fail on cache miss. Useful for reproducible CI runs
    /// against a pre-populated cache directory.
    #[arg(long, default_value_t = false)]
    pub offline: bool,

    /// Override endpoint (primarily for tests).
    #[arg(long, default_value = DEFAULT_ENDPOINT)]
    pub endpoint: String,

    /// Output TSV path.
    #[arg(long)]
    pub output: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CanonicalRequest {
    organism: String,
    query_name: String,
    genes: Vec<String>,
    background: Vec<String>,
    sources: Vec<String>,
    threshold_method: String,
    user_threshold: f64,
    ontology_version: String,
}

pub fn run_gprofiler(args: GprofilerArgs) -> Result<()> {
    if !(args.user_threshold.is_finite() && (0.0..=1.0).contains(&args.user_threshold)) {
        bail!("--user-threshold must be in [0, 1]");
    }
    let background = read_gene_list(&args.background)?
        .into_iter()
        .collect::<BTreeSet<_>>();
    if background.is_empty() {
        bail!("background {:?} yielded zero genes", args.background);
    }
    let sources: Vec<String> = args
        .sources
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    if sources.is_empty() {
        bail!("--sources must list at least one source");
    }

    let queries = if let Some(path) = &args.de_results {
        let genes = read_de_hits(path, args.comparison.as_deref(), args.user_threshold)?;
        if genes.is_empty() {
            bail!(
                "no genes with bh_q <= {} in {:?}",
                args.user_threshold,
                path
            );
        }
        vec![("de_hits".to_string(), genes)]
    } else if let Some(path) = &args.query_tsv {
        read_query_tsv(path)?
    } else {
        bail!("provide either --de-results or --query-tsv");
    };

    std::fs::create_dir_all(&args.cache_dir)
        .with_context(|| format!("creating cache dir {:?}", args.cache_dir))?;

    let mut all_rows: Vec<GprofilerRow> = Vec::new();
    for (query_name, genes) in queries {
        let genes_sorted: Vec<String> = genes.iter().cloned().collect::<BTreeSet<_>>().into_iter().collect();
        let background_sorted: Vec<String> = background.iter().cloned().collect();
        let canonical = CanonicalRequest {
            organism: args.organism.clone(),
            query_name: query_name.clone(),
            genes: genes_sorted.clone(),
            background: background_sorted.clone(),
            sources: sources.clone(),
            threshold_method: args.threshold_method.clone(),
            user_threshold: args.user_threshold,
            ontology_version: args.ontology_version.clone(),
        };
        let key = cache_key(&canonical);
        let cache_path = args.cache_dir.join(format!("{key}.json"));
        let response_json = if cache_path.exists() {
            std::fs::read_to_string(&cache_path)
                .with_context(|| format!("reading cache {:?}", cache_path))?
        } else if args.offline {
            bail!(
                "cache miss for query {:?} (key {}) in --offline mode; pre-populate {:?}",
                query_name,
                key,
                cache_path
            );
        } else {
            let body = build_request_body(&canonical);
            let response = post_gprofiler(&args.endpoint, &body)?;
            atomic_write(&cache_path, response.as_bytes())?;
            response
        };
        let parsed: Value = serde_json::from_str(&response_json)
            .with_context(|| format!("parsing g:Profiler JSON from {:?}", cache_path))?;
        let rows = parse_response(&query_name, &parsed)?;
        all_rows.extend(rows);
    }

    write_rows(&args.output, &all_rows)?;
    eprintln!(
        "enrich gprofiler: queries={} rows={} cache_dir={} output={}",
        all_rows
            .iter()
            .map(|r| r.query.clone())
            .collect::<BTreeSet<_>>()
            .len(),
        all_rows.len(),
        args.cache_dir.display(),
        args.output.display()
    );
    Ok(())
}

fn build_request_body(req: &CanonicalRequest) -> String {
    let body = json!({
        "organism": req.organism,
        "query": { req.query_name.clone(): req.genes.clone() },
        "sources": req.sources,
        "background": req.background,
        "user_threshold": req.user_threshold,
        "significance_threshold_method": req.threshold_method,
        "all_results": false,
        "ordered": false,
        "no_iea": false,
        "combined": false,
        "measure_underrepresentation": false,
        "domain_scope": "custom",
        "numeric_ns": "",
        "output": "json",
    });
    body.to_string()
}

fn post_gprofiler(endpoint: &str, body: &str) -> Result<String> {
    let response = ureq::post(endpoint)
        .set("Content-Type", "application/json")
        .set("User-Agent", "atman/1.0")
        .send_string(body)
        .with_context(|| format!("POST to g:Profiler endpoint {}", endpoint))?;
    response
        .into_string()
        .with_context(|| "reading g:Profiler response body".to_string())
}

fn cache_key(req: &CanonicalRequest) -> String {
    // Serialize canonically (sorted maps, deterministic order via BTreeMap).
    let mut map: BTreeMap<&str, Value> = BTreeMap::new();
    map.insert("organism", Value::String(req.organism.clone()));
    map.insert("query_name", Value::String(req.query_name.clone()));
    map.insert(
        "genes",
        Value::Array(req.genes.iter().map(|g| Value::String(g.clone())).collect()),
    );
    map.insert(
        "background",
        Value::Array(
            req.background
                .iter()
                .map(|g| Value::String(g.clone()))
                .collect(),
        ),
    );
    map.insert(
        "sources",
        Value::Array(req.sources.iter().map(|g| Value::String(g.clone())).collect()),
    );
    map.insert(
        "threshold_method",
        Value::String(req.threshold_method.clone()),
    );
    // Encode f64 with to_string for stable representation.
    map.insert("user_threshold", Value::String(req.user_threshold.to_string()));
    map.insert(
        "ontology_version",
        Value::String(req.ontology_version.clone()),
    );
    let canonical = serde_json::to_string(&map).expect("serialize cache key");
    sha256_hex(canonical.as_bytes())
}

fn read_gene_list(path: &Path) -> Result<Vec<String>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {:?}", path))?;
    let first = text.lines().next().unwrap_or("");
    if first.split('\t').any(|h| h == "gene_symbol") {
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_reader(text.as_bytes());
        let headers = reader.headers()?.clone();
        let gene_col = headers
            .iter()
            .position(|h| h == "gene_symbol")
            .with_context(|| format!("missing gene_symbol column in {:?}", path))?;
        let mut out = Vec::new();
        for row in reader.records() {
            let row = row?;
            let gene = row[gene_col].trim();
            if !gene.is_empty() {
                out.push(gene.to_string());
            }
        }
        Ok(out)
    } else {
        Ok(text
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty())
            .map(str::to_string)
            .collect())
    }
}

fn read_de_hits(path: &Path, comparison: Option<&str>, q_threshold: f64) -> Result<Vec<String>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let gene_col = headers
        .iter()
        .position(|h| h == "gene_symbol")
        .with_context(|| format!("missing gene_symbol in {:?}", path))?;
    let comparison_col = headers
        .iter()
        .position(|h| h == "comparison")
        .with_context(|| format!("missing comparison in {:?}", path))?;
    let q_col = headers
        .iter()
        .position(|h| h == "bh_q")
        .with_context(|| format!("missing bh_q in {:?}", path))?;
    let mut out = BTreeSet::new();
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
        if let Ok(q) = row[q_col].trim().parse::<f64>() {
            if q.is_finite() && q <= q_threshold {
                out.insert(gene.to_string());
            }
        }
    }
    Ok(out.into_iter().collect())
}

fn read_query_tsv(path: &Path) -> Result<Vec<(String, Vec<String>)>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let query_col = headers
        .iter()
        .position(|h| h == "query" || h == "program")
        .with_context(|| format!("missing query/program column in {:?}", path))?;
    let gene_col = headers
        .iter()
        .position(|h| h == "gene_symbol")
        .with_context(|| format!("missing gene_symbol column in {:?}", path))?;
    let mut by_query: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let query = row[query_col].trim();
        let gene = row[gene_col].trim();
        if query.is_empty() || gene.is_empty() {
            continue;
        }
        by_query
            .entry(query.to_string())
            .or_default()
            .insert(gene.to_string());
    }
    if by_query.is_empty() {
        bail!("no rows in query TSV {:?}", path);
    }
    Ok(by_query
        .into_iter()
        .map(|(k, v)| (k, v.into_iter().collect()))
        .collect())
}

#[derive(Debug, Clone)]
struct GprofilerRow {
    query: String,
    source: String,
    native: String,
    name: String,
    p_value: f64,
    intersection_size: i64,
    query_size: i64,
    term_size: i64,
    effective_domain_size: i64,
    intersections: String,
}

fn parse_response(query_name: &str, root: &Value) -> Result<Vec<GprofilerRow>> {
    let results = root
        .get("result")
        .and_then(|v| v.as_array())
        .with_context(|| format!("g:Profiler response missing `result` array for {query_name}"))?;
    let mut out = Vec::with_capacity(results.len());
    for entry in results {
        let query = entry
            .get("query")
            .and_then(|v| v.as_str())
            .unwrap_or(query_name);
        let source = entry
            .get("source")
            .and_then(|v| v.as_str())
            .unwrap_or("");
        let native = entry.get("native").and_then(|v| v.as_str()).unwrap_or("");
        let name = entry.get("name").and_then(|v| v.as_str()).unwrap_or("");
        let p_value = entry.get("p_value").and_then(|v| v.as_f64()).unwrap_or(f64::NAN);
        let intersection_size = entry
            .get("intersection_size")
            .and_then(|v| v.as_i64())
            .unwrap_or(0);
        let query_size = entry.get("query_size").and_then(|v| v.as_i64()).unwrap_or(0);
        let term_size = entry.get("term_size").and_then(|v| v.as_i64()).unwrap_or(0);
        let effective_domain_size = entry
            .get("effective_domain_size")
            .and_then(|v| v.as_i64())
            .unwrap_or(0);
        let intersections = entry
            .get("intersections")
            .and_then(|v| v.as_array())
            .map(|arr| {
                arr.iter()
                    .filter_map(|v| v.as_array())
                    .flat_map(|inner| inner.iter().filter_map(|s| s.as_str()))
                    .collect::<BTreeSet<&str>>()
                    .into_iter()
                    .collect::<Vec<_>>()
                    .join(";")
            })
            .unwrap_or_default();
        out.push(GprofilerRow {
            query: query.to_string(),
            source: source.to_string(),
            native: native.to_string(),
            name: name.to_string(),
            p_value,
            intersection_size,
            query_size,
            term_size,
            effective_domain_size,
            intersections,
        });
    }
    Ok(out)
}

fn write_rows(path: &Path, rows: &[GprofilerRow]) -> Result<()> {
    let mut out = String::from(
        "query\tsource\tnative\tname\tp_value\tintersection_size\tquery_size\tterm_size\teffective_domain_size\tintersections\n",
    );
    for row in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            row.query,
            row.source,
            row.native,
            escape_tsv(&row.name),
            row.p_value,
            row.intersection_size,
            row.query_size,
            row.term_size,
            row.effective_domain_size,
            row.intersections,
        ));
    }
    atomic_write(path, out.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cache_key_is_stable_for_same_request() {
        let a = CanonicalRequest {
            organism: "hsapiens".into(),
            query_name: "program_01".into(),
            genes: vec!["AQP4".into(), "GFAP".into()],
            background: vec!["AQP4".into(), "GFAP".into(), "SYN1".into()],
            sources: vec!["GO:BP".into()],
            threshold_method: "fdr".into(),
            user_threshold: 0.05,
            ontology_version: "pinned".into(),
        };
        let b = a.clone();
        assert_eq!(cache_key(&a), cache_key(&b));
    }

    #[test]
    fn cache_key_changes_with_ontology_version() {
        let a = CanonicalRequest {
            organism: "hsapiens".into(),
            query_name: "p".into(),
            genes: vec!["A".into()],
            background: vec!["A".into(), "B".into()],
            sources: vec!["GO:BP".into()],
            threshold_method: "fdr".into(),
            user_threshold: 0.05,
            ontology_version: "v1".into(),
        };
        let mut b = a.clone();
        b.ontology_version = "v2".into();
        assert_ne!(cache_key(&a), cache_key(&b));
    }

    #[test]
    fn parse_response_reads_canonical_columns() {
        let raw = json!({
            "result": [
                {
                    "query": "program_01",
                    "source": "GO:BP",
                    "native": "GO:0001234",
                    "name": "some biological process",
                    "p_value": 1.2e-5,
                    "intersection_size": 3,
                    "query_size": 10,
                    "term_size": 50,
                    "effective_domain_size": 20000,
                    "intersections": [["AQP4", "GFAP"], ["S100B"]]
                }
            ]
        });
        let rows = parse_response("program_01", &raw).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].native, "GO:0001234");
        assert!((rows[0].p_value - 1.2e-5).abs() < 1e-18);
        assert_eq!(rows[0].intersection_size, 3);
        assert_eq!(rows[0].query_size, 10);
        // Intersections flattened and deduplicated.
        let parts: BTreeSet<&str> = rows[0].intersections.split(';').collect();
        assert_eq!(parts, BTreeSet::from(["AQP4", "GFAP", "S100B"]));
    }
}
