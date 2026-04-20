//! `atman bench` — cross-tool benchmark harness.
//!
//! v1 scope: single `decompose` subcommand. Scores atman's FastICA
//! decomposition against a planted-archetype fixture and (optionally)
//! the outputs of named external tools invoked via thin shell
//! adapters under `bench/adapters/<tool>.sh`.
//!
//! Adapter contract: each adapter receives two arguments — a fixture
//! directory (with `planted_loadings.tsv` + `abundance.tsv` + optional
//! `cohorts.tsv`) and an output directory where it must write
//! `recovered_loadings.tsv` (columns: `archetype_id`, `protein`,
//! `loading`). The atman harness only orchestrates and scores; it
//! never executes Python/R directly. Missing or failing adapters
//! emit a `tool_not_available` row rather than aborting the run.

use anyhow::{bail, Context, Result};
use atman_core::bench_decompose::score_against_planted;
use atman_core::ica::{canonicalize_ica, fast_ica};
use clap::{Args as ClapArgs, Subcommand};
use csv::ReaderBuilder;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Instant, SystemTime};

use crate::io::{atomic_write, hash_labeled_inputs, sidecar_path_for, write_run_sidecar};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: BenchCommand,
}

#[derive(Subcommand, Debug)]
enum BenchCommand {
    /// Score `atman decompose ica` (+ optional external tools) against
    /// a planted-archetype fixture.
    Decompose(DecomposeArgs),
}

#[derive(ClapArgs, Debug)]
pub struct DecomposeArgs {
    /// Fixture directory containing `planted_loadings.tsv`
    /// (`archetype_id`, `protein`, `loading`) and `abundance.tsv`
    /// (first column `sample_id`, remaining columns one protein each).
    #[arg(long)]
    fixture: PathBuf,

    /// Comma-separated list of tool labels to score. `atman` always
    /// works (native); any other label triggers
    /// `bench/adapters/<label>.sh <fixture> <tmp_out>`. Missing or
    /// failing adapters emit `tool_not_available = 1` rather than
    /// aborting the run.
    #[arg(long, default_value = "atman")]
    tools: String,

    /// Directory containing adapter shell scripts. Convention:
    /// `<adapters_dir>/<tool>.sh`.
    #[arg(long)]
    adapters_dir: Option<PathBuf>,

    /// Seed forwarded to atman's FastICA (and echoed to adapter
    /// scripts via the `ATMAN_BENCH_SEED` env var).
    #[arg(long, default_value_t = 20260420)]
    seed: u64,

    /// Number of components `k` to ask each tool for. Defaults to the
    /// number of planted archetypes.
    #[arg(long)]
    k: Option<usize>,

    /// Top-N for the Jaccard overlap metric.
    #[arg(long, default_value_t = 50)]
    top_n: usize,

    /// FastICA max iterations.
    #[arg(long, default_value_t = 300)]
    max_iter: usize,

    /// FastICA convergence tolerance.
    #[arg(long, default_value_t = 1e-4)]
    tol: f64,

    /// Whether to run each tool twice and score
    /// `determinism_score = 1.0` if byte-equal, else 0.0.
    #[arg(long, default_value_t = true)]
    check_determinism: bool,

    /// Output TSV with one row per (tool, planted_archetype).
    #[arg(long)]
    output: PathBuf,
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        BenchCommand::Decompose(args) => run_decompose(args),
    }
}

fn run_decompose(args: DecomposeArgs) -> Result<()> {
    let started_at = SystemTime::now();
    let fixture = &args.fixture;
    let (protein_labels, planted) = read_planted_loadings(fixture)?;
    let (sample_ids, abundance) = read_abundance(fixture, &protein_labels)?;
    let k = args.k.unwrap_or(planted.len());
    if k == 0 {
        bail!("fixture has zero planted archetypes");
    }

    let tools: Vec<String> = args
        .tools
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    if tools.is_empty() {
        bail!("--tools resolved to zero entries");
    }

    let mut rows: Vec<ToolRunResult> = Vec::new();
    for tool in &tools {
        let result = if tool == "atman" {
            run_atman_native(&abundance, k, args.seed, args.max_iter, args.tol, &protein_labels)
        } else {
            let adapters_dir = match &args.adapters_dir {
                Some(d) => d.clone(),
                None => {
                    rows.push(ToolRunResult::unavailable(
                        tool.clone(),
                        "missing --adapters-dir".into(),
                    ));
                    continue;
                }
            };
            run_external_adapter(tool, fixture, &adapters_dir, args.seed, k)
        };
        let mut tr = match result {
            Ok(t) => t,
            Err(e) => {
                rows.push(ToolRunResult::unavailable(tool.clone(), e.to_string()));
                continue;
            }
        };

        if args.check_determinism {
            // Re-run and byte-compare the recovered loadings.
            let second = if tool == "atman" {
                run_atman_native(
                    &abundance,
                    k,
                    args.seed,
                    args.max_iter,
                    args.tol,
                    &protein_labels,
                )
                .ok()
            } else if let Some(d) = &args.adapters_dir {
                run_external_adapter(tool, fixture, d, args.seed, k).ok()
            } else {
                None
            };
            tr.determinism_score = match &second {
                Some(s) if s.recovered_flat == tr.recovered_flat => 1.0,
                Some(_) => 0.0,
                None => f64::NAN,
            };
        }
        rows.push(tr);
    }

    // Score recovered vs planted per tool.
    let mut result_rows = Vec::new();
    for tool in &rows {
        if tool.tool_not_available {
            // Emit a single "unavailable" row so the benchmark TSV
            // records every requested tool.
            result_rows.push(BenchRow {
                tool: tool.tool.clone(),
                planted_archetype: 0,
                matched_recovered: None,
                archetype_correlation: f64::NAN,
                recovery_jaccard: f64::NAN,
                matched_cosine: f64::NAN,
                runtime_seconds: f64::NAN,
                determinism_score: f64::NAN,
                tool_not_available: 1,
                unavailable_reason: tool.unavailable_reason.clone(),
            });
            continue;
        }
        let scores = score_against_planted(&planted, &tool.recovered, args.top_n);
        for s in scores {
            result_rows.push(BenchRow {
                tool: tool.tool.clone(),
                planted_archetype: s.planted_archetype,
                matched_recovered: s.matched_recovered,
                archetype_correlation: s.archetype_correlation,
                recovery_jaccard: s.recovery_jaccard,
                matched_cosine: s.matched_cosine,
                runtime_seconds: tool.runtime_seconds,
                determinism_score: tool.determinism_score,
                tool_not_available: 0,
                unavailable_reason: String::new(),
            });
        }
    }

    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating {:?}", parent))?;
    }
    write_bench_results(&args.output, &result_rows)?;
    eprintln!(
        "bench decompose: tools={} planted_k={} samples={} proteins={}",
        tools.len(),
        planted.len(),
        sample_ids.len(),
        protein_labels.len(),
    );

    let finished_at = SystemTime::now();
    let planted_path = fixture.join("planted_loadings.tsv");
    let abundance_path = fixture.join("abundance.tsv");
    let input_dir_sha256 = hash_labeled_inputs(&[
        ("planted_loadings.tsv", planted_path.as_path()),
        ("abundance.tsv", abundance_path.as_path()),
    ])?;
    let sidecar = sidecar_path_for(&args.output);
    write_run_sidecar(
        &sidecar,
        "bench decompose",
        json!({
            "fixture": args.fixture.display().to_string(),
            "tools": args.tools,
            "adapters-dir": args.adapters_dir.as_ref().map(|p| p.display().to_string()),
            "seed": args.seed,
            "k": k,
            "top-n": args.top_n,
            "max-iter": args.max_iter,
            "tol": args.tol,
            "check-determinism": args.check_determinism,
        }),
        &input_dir_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
    )?;
    eprintln!("bench decompose: sidecar={}", sidecar.display());
    Ok(())
}

struct ToolRunResult {
    tool: String,
    recovered: Vec<Vec<f64>>,
    recovered_flat: Vec<u8>,
    runtime_seconds: f64,
    determinism_score: f64,
    tool_not_available: bool,
    unavailable_reason: String,
}

impl ToolRunResult {
    fn unavailable(tool: String, reason: String) -> Self {
        Self {
            tool,
            recovered: Vec::new(),
            recovered_flat: Vec::new(),
            runtime_seconds: f64::NAN,
            determinism_score: f64::NAN,
            tool_not_available: true,
            unavailable_reason: reason,
        }
    }
}

fn run_atman_native(
    abundance: &[Vec<f64>],
    k: usize,
    seed: u64,
    max_iter: usize,
    tol: f64,
    protein_labels: &[String],
) -> Result<ToolRunResult> {
    let start = Instant::now();
    let ica = fast_ica(abundance, k, seed, max_iter, tol);
    let canon = canonicalize_ica(&ica);
    let elapsed = start.elapsed().as_secs_f64();
    let recovered = canon.loadings.clone();
    // Stable byte representation for determinism-check: concatenate
    // every (archetype, protein, rounded-loading) triple in fixed order.
    let mut flat = Vec::with_capacity(recovered.len() * protein_labels.len() * 16);
    for (ai, row) in recovered.iter().enumerate() {
        for (pi, &v) in row.iter().enumerate() {
            flat.extend_from_slice(
                format!("{}\t{}\t{:.10}\n", ai, protein_labels[pi], v).as_bytes(),
            );
        }
    }
    Ok(ToolRunResult {
        tool: "atman".to_string(),
        recovered,
        recovered_flat: flat,
        runtime_seconds: elapsed,
        determinism_score: f64::NAN,
        tool_not_available: false,
        unavailable_reason: String::new(),
    })
}

fn run_external_adapter(
    tool: &str,
    fixture: &Path,
    adapters_dir: &Path,
    seed: u64,
    k: usize,
) -> Result<ToolRunResult> {
    let script = adapters_dir.join(format!("{tool}.sh"));
    if !script.exists() {
        bail!("adapter {script:?} does not exist");
    }
    let tmp = tempfile::tempdir().context("creating bench tmp dir")?;
    let out_dir = tmp.path().to_path_buf();
    let start = Instant::now();
    let status = Command::new(&script)
        .arg(fixture)
        .arg(&out_dir)
        .env("ATMAN_BENCH_SEED", seed.to_string())
        .env("ATMAN_BENCH_K", k.to_string())
        .status()
        .with_context(|| format!("running adapter {script:?}"))?;
    if !status.success() {
        bail!("adapter {tool} exited non-zero");
    }
    let elapsed = start.elapsed().as_secs_f64();
    let recovered_path = out_dir.join("recovered_loadings.tsv");
    if !recovered_path.exists() {
        bail!("adapter {tool} did not produce recovered_loadings.tsv");
    }
    let (_, recovered) = read_recovered_loadings(&recovered_path)?;
    let flat = std::fs::read(&recovered_path)?;
    Ok(ToolRunResult {
        tool: tool.to_string(),
        recovered,
        recovered_flat: flat,
        runtime_seconds: elapsed,
        determinism_score: f64::NAN,
        tool_not_available: false,
        unavailable_reason: String::new(),
    })
}

fn read_planted_loadings(fixture: &Path) -> Result<(Vec<String>, Vec<Vec<f64>>)> {
    let path = fixture.join("planted_loadings.tsv");
    let (proteins, loadings) = read_recovered_loadings(&path)
        .with_context(|| format!("reading {:?}", path))?;
    if loadings.is_empty() {
        bail!("planted_loadings.tsv has no archetypes");
    }
    Ok((proteins, loadings))
}

/// Read a loadings TSV in the universal bench format
/// (`archetype_id`, `protein`, `loading`) and return
/// `(ordered_protein_labels, [archetype][protein])`.
fn read_recovered_loadings(path: &Path) -> Result<(Vec<String>, Vec<Vec<f64>>)> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let arch_idx = headers
        .iter()
        .position(|h| h == "archetype_id")
        .context("loadings TSV missing archetype_id column")?;
    let prot_idx = headers
        .iter()
        .position(|h| h == "protein")
        .context("loadings TSV missing protein column")?;
    let load_idx = headers
        .iter()
        .position(|h| h == "loading")
        .context("loadings TSV missing loading column")?;

    // Preserve first-seen protein order.
    let mut protein_order: Vec<String> = Vec::new();
    let mut protein_idx: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();
    let mut by_arch: BTreeMap<String, BTreeMap<String, f64>> = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let arch = row[arch_idx].to_string();
        let protein = row[prot_idx].to_string();
        let value: f64 = row[load_idx].parse().unwrap_or(0.0);
        if !protein_idx.contains_key(&protein) {
            protein_idx.insert(protein.clone(), protein_order.len());
            protein_order.push(protein.clone());
        }
        by_arch.entry(arch).or_default().insert(protein, value);
    }
    let arch_keys: Vec<String> = by_arch.keys().cloned().collect();
    let mut out = Vec::with_capacity(arch_keys.len());
    for a in &arch_keys {
        let map = &by_arch[a];
        let mut row = vec![0.0_f64; protein_order.len()];
        for (p, v) in map {
            if let Some(&i) = protein_idx.get(p) {
                row[i] = *v;
            }
        }
        out.push(row);
    }
    Ok((protein_order, out))
}

fn read_abundance(
    fixture: &Path,
    protein_order: &[String],
) -> Result<(Vec<String>, Vec<Vec<f64>>)> {
    let path = fixture.join("abundance.tsv");
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(&path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let sid_idx = headers
        .iter()
        .position(|h| h == "sample_id")
        .context("abundance.tsv missing sample_id column")?;
    // Column index per planted-protein label.
    let mut col_map: Vec<usize> = Vec::with_capacity(protein_order.len());
    let proteins_in_file: BTreeSet<&str> = headers.iter().collect();
    for p in protein_order {
        let idx = headers.iter().position(|h| h == p).with_context(|| {
            format!(
                "abundance.tsv missing planted protein {p:?}; columns were {:?}",
                proteins_in_file
            )
        })?;
        col_map.push(idx);
    }
    let mut sample_ids = Vec::new();
    let mut data = Vec::new();
    for row in reader.records() {
        let row = row?;
        sample_ids.push(row[sid_idx].to_string());
        let values: Vec<f64> = col_map
            .iter()
            .map(|&c| row[c].parse::<f64>().unwrap_or(0.0))
            .collect();
        data.push(values);
    }
    Ok((sample_ids, data))
}

struct BenchRow {
    tool: String,
    planted_archetype: usize,
    matched_recovered: Option<usize>,
    archetype_correlation: f64,
    recovery_jaccard: f64,
    matched_cosine: f64,
    runtime_seconds: f64,
    determinism_score: f64,
    tool_not_available: u8,
    unavailable_reason: String,
}

fn write_bench_results(path: &Path, rows: &[BenchRow]) -> Result<()> {
    let mut buf = String::from(
        "tool\tplanted_archetype\tmatched_recovered\tarchetype_correlation\t\
         recovery_jaccard\tmatched_cosine\truntime_seconds\tdeterminism_score\t\
         tool_not_available\tunavailable_reason\n",
    );
    for r in rows {
        let fmt_f = |v: f64| if v.is_finite() { format!("{v:.6}") } else { "NA".into() };
        buf.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.tool,
            r.planted_archetype,
            r.matched_recovered
                .map(|i| i.to_string())
                .unwrap_or_else(|| "NA".into()),
            fmt_f(r.archetype_correlation),
            fmt_f(r.recovery_jaccard),
            fmt_f(r.matched_cosine),
            fmt_f(r.runtime_seconds),
            fmt_f(r.determinism_score),
            r.tool_not_available,
            r.unavailable_reason,
        ));
    }
    atomic_write(path, buf.as_bytes())
}
