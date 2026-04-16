use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;
use csv::ReaderBuilder;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::PathBuf;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Baseline de_results.tsv.
    #[arg(long)]
    baseline: PathBuf,

    /// Comma-separated list of LOO de_results.tsv files.
    #[arg(long)]
    loo: String,

    /// Output directory for robustness tables.
    #[arg(long)]
    output_dir: PathBuf,

    /// Top-k for ranking overlap.
    #[arg(long, default_value_t = 20)]
    top_k: usize,
}

#[derive(Debug, Clone)]
struct DeLite {
    mean_diff: Option<f64>,
    bh_q: Option<f64>,
}

pub fn run(args: Args) -> Result<()> {
    if args.top_k == 0 {
        bail!("top-k must be >= 1");
    }
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating {:?}", args.output_dir))?;

    let loo_files = parse_csv_list(&args.loo)?;
    let base = read_de(&args.baseline)?;
    let loos: Vec<_> = loo_files
        .iter()
        .map(|p| read_de(&PathBuf::from(p)))
        .collect::<Result<Vec<_>>>()?;

    let mut sign_out = String::from(
        "comparison\tpanel\tassay_id\tn_loo\t\
         n_sign_match\tsign_match_rate\t\
         baseline_q\tretained_q_lt_05\tn_retained_q_lt_05\tn_retained_q_lt_10\n",
    );
    for (cmp, base_rows) in &base {
        for (key, b) in base_rows {
            let mut n_loo = 0usize;
            let mut n_sign_match = 0usize;
            let mut n_ret05 = 0usize;
            let mut n_ret10 = 0usize;
            for loo in &loos {
                if let Some(r) = loo.get(cmp).and_then(|m| m.get(key)) {
                    if let (Some(bd), Some(ld)) = (b.mean_diff, r.mean_diff) {
                        n_loo += 1;
                        if bd != 0.0 && ld != 0.0 && bd.signum() == ld.signum() {
                            n_sign_match += 1;
                        }
                    }
                    if matches!(r.bh_q, Some(q) if q < 0.05) {
                        n_ret05 += 1;
                    }
                    if matches!(r.bh_q, Some(q) if q < 0.10) {
                        n_ret10 += 1;
                    }
                }
            }
            let sign_rate = if n_loo == 0 {
                None
            } else {
                Some(n_sign_match as f64 / n_loo as f64)
            };
            let base_q = b.bh_q;
            let retained_q_lt_05 = matches!(base_q, Some(q) if q < 0.05);
            let (panel, assay) = split_key(key);
            sign_out.push_str(cmp);
            sign_out.push('\t');
            sign_out.push_str(panel);
            sign_out.push('\t');
            sign_out.push_str(assay);
            sign_out.push('\t');
            sign_out.push_str(&n_loo.to_string());
            sign_out.push('\t');
            sign_out.push_str(&n_sign_match.to_string());
            sign_out.push('\t');
            push_opt_f64(&mut sign_out, sign_rate);
            sign_out.push('\t');
            push_opt_f64(&mut sign_out, base_q);
            sign_out.push('\t');
            sign_out.push_str(if retained_q_lt_05 { "true" } else { "false" });
            sign_out.push('\t');
            sign_out.push_str(&n_ret05.to_string());
            sign_out.push('\t');
            sign_out.push_str(&n_ret10.to_string());
            sign_out.push('\n');
        }
    }

    let mut rank_out =
        String::from("comparison\tloo_index\ttop_k\toverlap_count\tjaccard_index\n");
    for (cmp, base_rows) in &base {
        let base_top = top_k_set(base_rows, args.top_k);
        for (i, loo) in loos.iter().enumerate() {
            let loo_rows = match loo.get(cmp) {
                Some(v) => v,
                None => continue,
            };
            let loo_top = top_k_set(loo_rows, args.top_k);
            let inter = base_top.intersection(&loo_top).count();
            let union = base_top.union(&loo_top).count();
            let j = if union == 0 {
                None
            } else {
                Some(inter as f64 / union as f64)
            };
            rank_out.push_str(cmp);
            rank_out.push('\t');
            rank_out.push_str(&(i + 1).to_string());
            rank_out.push('\t');
            rank_out.push_str(&args.top_k.to_string());
            rank_out.push('\t');
            rank_out.push_str(&inter.to_string());
            rank_out.push('\t');
            push_opt_f64(&mut rank_out, j);
            rank_out.push('\n');
        }
    }

    std::fs::write(args.output_dir.join("loo_sign_stability.tsv"), sign_out.as_bytes())
        .with_context(|| "writing loo_sign_stability.tsv")?;
    std::fs::write(args.output_dir.join("rank_stability.tsv"), rank_out.as_bytes())
        .with_context(|| "writing rank_stability.tsv")?;
    eprintln!(
        "robustness: wrote {:?} and {:?}",
        args.output_dir.join("loo_sign_stability.tsv"),
        args.output_dir.join("rank_stability.tsv")
    );
    Ok(())
}

fn parse_csv_list(raw: &str) -> Result<Vec<String>> {
    let mut out = Vec::new();
    for t in raw.split(',') {
        let v = t.trim();
        if v.is_empty() {
            bail!("invalid empty entry in --loo {:?}", raw);
        }
        out.push(v.to_string());
    }
    if out.is_empty() {
        bail!("no LOO files provided");
    }
    Ok(out)
}

fn read_de(path: &PathBuf) -> Result<BTreeMap<String, HashMap<String, DeLite>>> {
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = rdr.headers()?.clone();
    let need = |name: &str| -> Result<usize> {
        headers
            .iter()
            .position(|h| h == name)
            .with_context(|| format!("missing column {:?} in {:?}", name, path))
    };
    let c_cmp = need("comparison")?;
    let c_panel = need("panel")?;
    let c_assay = need("assay_id")?;
    let c_diff = need("mean_diff")?;
    let c_q = need("bh_q")?;

    let mut out: BTreeMap<String, HashMap<String, DeLite>> = BTreeMap::new();
    for rec in rdr.records() {
        let row = rec?;
        let cmp = row[c_cmp].to_string();
        let key = format!("{}::{}", row[c_panel].trim(), row[c_assay].trim());
        out.entry(cmp).or_default().insert(
            key,
            DeLite {
                mean_diff: parse_opt_f64(&row[c_diff])?,
                bh_q: parse_opt_f64(&row[c_q])?,
            },
        );
    }
    Ok(out)
}

fn parse_opt_f64(s: &str) -> Result<Option<f64>> {
    if s.trim().is_empty() {
        return Ok(None);
    }
    Ok(Some(s.trim().parse().with_context(|| format!("parse f64 {:?}", s))?))
}

fn top_k_set(rows: &HashMap<String, DeLite>, k: usize) -> BTreeSet<String> {
    let mut v: Vec<(&String, f64)> = rows
        .iter()
        .filter_map(|(k, r)| r.mean_diff.map(|d| (k, d.abs())))
        .collect();
    v.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    v.into_iter().take(k).map(|(k, _)| k.clone()).collect()
}

fn split_key(key: &str) -> (&str, &str) {
    let mut it = key.splitn(2, "::");
    let a = it.next().unwrap_or_default();
    let b = it.next().unwrap_or_default();
    (a, b)
}

fn push_opt_f64(buf: &mut String, v: Option<f64>) {
    if let Some(x) = v {
        buf.push_str(&x.to_string());
    }
}

