use anyhow::{bail, Context, Result};
use clap::{Args as ClapArgs, Subcommand};
use csv::ReaderBuilder;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use crate::io::{atomic_write, escape_tsv, find_col, format_float, need_col, optional_cell};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Flag ICA programs that fail annotation, loading, or keratin filters.
    Filter(FilterArgs),
}

#[derive(ClapArgs, Debug)]
pub struct FilterArgs {
    /// Program loading TSV with program, gene_symbol/protein/assay_id, and loading columns.
    #[arg(long)]
    loadings: PathBuf,

    /// Program annotation TSV with program and top_annotation_p_value columns.
    #[arg(long)]
    annotations: PathBuf,

    /// Maximum annotation p-value required for an interpretable program.
    #[arg(long, default_value_t = 0.05)]
    min_annotation_pvalue: f64,

    /// Minimum absolute top loading required for an interpretable program.
    #[arg(long, default_value_t = 0.5)]
    min_top_loading: f64,

    /// Maximum fraction of top loadings that may be keratin-like.
    #[arg(long, default_value_t = 0.2)]
    max_keratin_fraction: f64,

    /// Number of strongest absolute-loading proteins used for keratin fraction.
    #[arg(long, default_value_t = 20)]
    top_n: usize,

    /// Output program interpretability TSV.
    #[arg(long)]
    output: PathBuf,
}

#[derive(Debug, Clone)]
struct Loading {
    label: String,
    loading: f64,
}

#[derive(Debug, Clone)]
struct Annotation {
    category: String,
    top_annotation: String,
    top_annotation_p_value: f64,
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Filter(args) => run_filter(args),
    }
}

fn run_filter(args: FilterArgs) -> Result<()> {
    if !args.min_annotation_pvalue.is_finite() || args.min_annotation_pvalue < 0.0 {
        bail!("--min-annotation-pvalue must be a non-negative finite value");
    }
    if !args.min_top_loading.is_finite() || args.min_top_loading < 0.0 {
        bail!("--min-top-loading must be a non-negative finite value");
    }
    if !args.max_keratin_fraction.is_finite()
        || args.max_keratin_fraction < 0.0
        || args.max_keratin_fraction > 1.0
    {
        bail!("--max-keratin-fraction must be in [0, 1]");
    }
    if args.top_n == 0 {
        bail!("--top-n must be at least 1");
    }

    let loadings = read_loadings(&args.loadings)?;
    let annotations = read_annotations(&args.annotations)?;
    let mut out = String::from(
        "program\tn_loadings\tmax_abs_loading\tkeratin_fraction_top_n\tannotation_p_value\tcategory\ttop_annotation\tannotation_pass\tloading_pass\tkeratin_pass\tinterpretable\tfail_reasons\n",
    );

    for (program, rows) in &loadings {
        let n_loadings = rows.len();
        let max_abs_loading = rows
            .iter()
            .map(|row| row.loading.abs())
            .fold(0.0_f64, f64::max);
        let keratin_fraction = keratin_fraction(rows, args.top_n);
        let annotation = annotations.get(program);
        let annotation_p_value = annotation.map(|a| a.top_annotation_p_value);
        let annotation_pass = annotation_p_value
            .map(|p| p <= args.min_annotation_pvalue)
            .unwrap_or(false);
        let loading_pass = max_abs_loading >= args.min_top_loading;
        let keratin_pass = keratin_fraction <= args.max_keratin_fraction;
        let interpretable = annotation_pass && loading_pass && keratin_pass;
        let mut fail_reasons = Vec::new();
        if annotation.is_none() {
            fail_reasons.push("missing_annotation");
        } else if !annotation_pass {
            fail_reasons.push("annotation_p_value");
        }
        if !loading_pass {
            fail_reasons.push("diffuse_loading");
        }
        if !keratin_pass {
            fail_reasons.push("keratin_contamination");
        }

        let category = annotation.map(|a| a.category.as_str()).unwrap_or("NA");
        let top_annotation = annotation
            .map(|a| a.top_annotation.as_str())
            .unwrap_or("NA");
        let annotation_p = annotation_p_value
            .map(format_float)
            .unwrap_or_else(|| "NA".to_string());
        out.push_str(&format!(
            "{program}\t{n_loadings}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            format_float(max_abs_loading),
            format_float(keratin_fraction),
            annotation_p,
            escape_tsv(category),
            escape_tsv(top_annotation),
            flag(annotation_pass),
            flag(loading_pass),
            flag(keratin_pass),
            flag(interpretable),
            if fail_reasons.is_empty() {
                "none".to_string()
            } else {
                fail_reasons.join(";")
            }
        ));
    }

    atomic_write(&args.output, out.as_bytes())?;
    eprintln!(
        "programs filter: programs={} output={}",
        loadings.len(),
        args.output.display()
    );
    Ok(())
}

fn read_loadings(path: &Path) -> Result<BTreeMap<String, Vec<Loading>>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let program_col = need_col(&headers, "program", path)?;
    let label_col = find_col(&headers, &["gene_symbol", "protein", "assay_id"])
        .with_context(|| format!("{:?} missing gene_symbol/protein/assay_id column", path))?;
    let loading_col = need_col(&headers, "loading", path)?;
    let mut out: BTreeMap<String, Vec<Loading>> = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let loading: f64 = row[loading_col]
            .parse()
            .with_context(|| format!("parsing loading in {:?}", path))?;
        if !loading.is_finite() {
            continue;
        }
        out.entry(row[program_col].to_string())
            .or_default()
            .push(Loading {
                label: row[label_col].to_string(),
                loading,
            });
    }
    if out.is_empty() {
        bail!("no finite program loadings found in {:?}", path);
    }
    Ok(out)
}

fn read_annotations(path: &Path) -> Result<BTreeMap<String, Annotation>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let program_col = need_col(&headers, "program", path)?;
    let category_col = find_col(&headers, &["category"]);
    let top_annotation_col = find_col(&headers, &["top_annotation", "name", "term"]);
    let p_col = need_col(&headers, "top_annotation_p_value", path)?;
    let mut out = BTreeMap::new();
    for row in reader.records() {
        let row = row?;
        let p_value: f64 = row[p_col]
            .parse()
            .with_context(|| format!("parsing top_annotation_p_value in {:?}", path))?;
        if !p_value.is_finite() {
            continue;
        }
        out.insert(
            row[program_col].to_string(),
            Annotation {
                category: optional_cell(&row, category_col)
                    .unwrap_or("NA")
                    .to_string(),
                top_annotation: optional_cell(&row, top_annotation_col)
                    .unwrap_or("NA")
                    .to_string(),
                top_annotation_p_value: p_value,
            },
        );
    }
    Ok(out)
}

fn keratin_fraction(rows: &[Loading], top_n: usize) -> f64 {
    let mut ranked = rows.to_vec();
    ranked.sort_by(|a, b| {
        b.loading
            .abs()
            .partial_cmp(&a.loading.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let keep = ranked.len().min(top_n);
    if keep == 0 {
        return 0.0;
    }
    let keratin = ranked
        .iter()
        .take(keep)
        .filter(|row| is_keratin_like(&row.label))
        .count();
    keratin as f64 / keep as f64
}

fn is_keratin_like(label: &str) -> bool {
    let upper = label.to_ascii_uppercase();
    upper.starts_with("KRT") || upper.contains("KERATIN")
}

fn flag(value: bool) -> &'static str {
    if value {
        "1"
    } else {
        "0"
    }
}
