use anyhow::{bail, Context, Result};
use atman_core::align::{
    build_archetypes, summarize_archetypes, AlignMetric, AlignedProgram,
};
use clap::{Args as ClapArgs, Subcommand, ValueEnum};
use csv::ReaderBuilder;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use crate::io::atomic_write;

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Cross-cohort program alignment and sweep mode.
    Programs(ProgramsArgs),
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
    Spearman,
}

impl From<SingleMetric> for AlignMetric {
    fn from(m: SingleMetric) -> Self {
        match m {
            SingleMetric::Jaccard => AlignMetric::Jaccard,
            SingleMetric::Cosine => AlignMetric::Cosine,
            SingleMetric::Spearman => AlignMetric::Spearman,
        }
    }
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Programs(args) => run_programs(args),
    }
}

fn run_programs(args: ProgramsArgs) -> Result<()> {
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
        None => paths
            .iter()
            .map(|p| file_stem(p).to_string())
            .collect(),
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

    if args.sweep {
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
    } else {
        let output = args
            .output
            .as_ref()
            .context("single-mode requires --output")?;
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
    }
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
        out.push((
            row[program_col].to_string(),
            label.to_string(),
            value,
        ));
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
        let parts = split_top_level(&block, ':');
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
