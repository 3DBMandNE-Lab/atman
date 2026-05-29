use anyhow::{bail, Context, Result};
use atman_core::de::stability_weighted_score;
use clap::Args as ClapArgs;
use csv::ReaderBuilder;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::io::{hash_labeled_inputs, sidecar_path_for, write_run_sidecar};
use serde_json::Map as JsonMap;

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

    /// Strict BH-q threshold for the `retained_q_strict` /
    /// `n_retained_q_strict` columns. Reporting policy only — does
    /// not change which features are read or written. Must satisfy
    /// `0 < strict < relaxed < 1`.
    #[arg(long, default_value_t = 0.05)]
    report_q_strict: f64,

    /// Relaxed BH-q threshold for the `n_retained_q_relaxed` column.
    /// Same scope as `--report-q-strict`.
    #[arg(long, default_value_t = 0.10)]
    report_q_relaxed: f64,
}

#[derive(Debug, Clone)]
struct DeLite {
    gene_symbol: String,
    mean_diff: Option<f64>,
    bh_q: Option<f64>,
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    if args.top_k == 0 {
        bail!("top-k must be >= 1");
    }
    if !args.report_q_strict.is_finite() || !(0.0..1.0).contains(&args.report_q_strict) {
        bail!(
            "--report-q-strict {} must be in (0, 1)",
            args.report_q_strict
        );
    }
    if !args.report_q_relaxed.is_finite() || !(0.0..1.0).contains(&args.report_q_relaxed) {
        bail!(
            "--report-q-relaxed {} must be in (0, 1)",
            args.report_q_relaxed
        );
    }
    if args.report_q_strict >= args.report_q_relaxed {
        bail!(
            "--report-q-strict {} must be < --report-q-relaxed {}",
            args.report_q_strict,
            args.report_q_relaxed
        );
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
         baseline_q\tretained_q_strict\tn_retained_q_strict\tn_retained_q_relaxed\n",
    );
    for (cmp, base_rows) in &base {
        // `base_rows` is a HashMap; iterate in sorted-key order so the emitted
        // row order is deterministic (same class of bug as the top-k tie-break).
        let mut keys: Vec<&String> = base_rows.keys().collect();
        keys.sort();
        for key in keys {
            let b = &base_rows[key];
            let mut n_loo = 0usize;
            let mut n_sign_match = 0usize;
            let mut n_ret_strict = 0usize;
            let mut n_ret_relaxed = 0usize;
            for loo in &loos {
                if let Some(r) = loo.get(cmp).and_then(|m| m.get(key)) {
                    if let (Some(bd), Some(ld)) = (b.mean_diff, r.mean_diff) {
                        n_loo += 1;
                        if bd != 0.0 && ld != 0.0 && bd.signum() == ld.signum() {
                            n_sign_match += 1;
                        }
                    }
                    if matches!(r.bh_q, Some(q) if q < args.report_q_strict) {
                        n_ret_strict += 1;
                    }
                    if matches!(r.bh_q, Some(q) if q < args.report_q_relaxed) {
                        n_ret_relaxed += 1;
                    }
                }
            }
            let sign_rate = if n_loo == 0 {
                None
            } else {
                Some(n_sign_match as f64 / n_loo as f64)
            };
            let base_q = b.bh_q;
            let retained_q_strict = matches!(base_q, Some(q) if q < args.report_q_strict);
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
            sign_out.push_str(if retained_q_strict { "true" } else { "false" });
            sign_out.push('\t');
            sign_out.push_str(&n_ret_strict.to_string());
            sign_out.push('\t');
            sign_out.push_str(&n_ret_relaxed.to_string());
            sign_out.push('\n');
        }
    }

    let mut rank_out = String::from("comparison\tloo_index\ttop_k\toverlap_count\tjaccard_index\n");
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

    // ---- Stability-weighted ranking -------------------------------------
    // S_i = |d_i| · SSR_i · (1 - q_i). Per contrast, we also report rank_by_q
    // (1 = most significant) and rank_by_s (1 = highest stability-weighted
    // score) so downstream analyses can compare the two orderings directly.
    type StabRow = (
        String,
        String,
        String,
        Option<f64>,
        Option<f64>,
        Option<f64>,
        Option<f64>,
    );
    let mut stab_by_cmp: BTreeMap<String, Vec<StabRow>> = BTreeMap::new();
    for (cmp, base_rows) in &base {
        let mut rows: Vec<StabRow> = Vec::new();
        // Iterate the HashMap in sorted-key order so `rows` has a stable initial
        // order; the rank sorts below use a stable sort, so ties resolve to this
        // deterministic feature-id order.
        let mut keys: Vec<&String> = base_rows.keys().collect();
        keys.sort();
        for key in keys {
            let b = &base_rows[key];
            let mut n_loo = 0usize;
            let mut n_sign_match = 0usize;
            for loo in &loos {
                if let Some(r) = loo.get(cmp).and_then(|m| m.get(key)) {
                    if let (Some(bd), Some(ld)) = (b.mean_diff, r.mean_diff) {
                        n_loo += 1;
                        if bd != 0.0 && ld != 0.0 && bd.signum() == ld.signum() {
                            n_sign_match += 1;
                        }
                    }
                }
            }
            let sign_rate = if n_loo == 0 {
                None
            } else {
                Some(n_sign_match as f64 / n_loo as f64)
            };
            let score = stability_weighted_score(b.mean_diff, b.bh_q, sign_rate);
            let (panel, _assay) = split_key(key);
            rows.push((
                panel.to_string(),
                b.gene_symbol.clone(),
                key.clone(),
                b.mean_diff,
                b.bh_q,
                sign_rate,
                score,
            ));
        }
        stab_by_cmp.insert(cmp.clone(), rows);
    }

    let mut stab_out = String::from(
        "comparison\tpanel\tgene_symbol\tassay_id\t\
         baseline_mean_diff\tbaseline_bh_q\tsign_match_rate\t\
         stability_score\trank_by_q\trank_by_stability\n",
    );
    for (cmp, rows) in stab_by_cmp.iter() {
        // rank_by_q: ascending bh_q with None last → rank 1 = smallest q.
        let mut idx_q: Vec<usize> = (0..rows.len()).collect();
        idx_q.sort_by(|&a, &b| match (rows[a].4, rows[b].4) {
            (Some(qa), Some(qb)) => qa.partial_cmp(&qb).unwrap_or(std::cmp::Ordering::Equal),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => std::cmp::Ordering::Equal,
        });
        let mut rank_q = vec![0usize; rows.len()];
        for (r, &i) in idx_q.iter().enumerate() {
            rank_q[i] = r + 1;
        }
        // rank_by_s: descending stability score with None last → rank 1 = highest S.
        let mut idx_s: Vec<usize> = (0..rows.len()).collect();
        idx_s.sort_by(|&a, &b| match (rows[a].6, rows[b].6) {
            (Some(sa), Some(sb)) => sb.partial_cmp(&sa).unwrap_or(std::cmp::Ordering::Equal),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => std::cmp::Ordering::Equal,
        });
        let mut rank_s = vec![0usize; rows.len()];
        for (r, &i) in idx_s.iter().enumerate() {
            rank_s[i] = r + 1;
        }
        // Write rows sorted by stability_score descending for inspection.
        for &i in idx_s.iter() {
            let (panel, gene, key, d, q, ssr, s) = &rows[i];
            let (_p, assay) = split_key(key);
            stab_out.push_str(cmp);
            stab_out.push('\t');
            stab_out.push_str(panel);
            stab_out.push('\t');
            stab_out.push_str(gene);
            stab_out.push('\t');
            stab_out.push_str(assay);
            stab_out.push('\t');
            push_opt_f64(&mut stab_out, *d);
            stab_out.push('\t');
            push_opt_f64(&mut stab_out, *q);
            stab_out.push('\t');
            push_opt_f64(&mut stab_out, *ssr);
            stab_out.push('\t');
            push_opt_f64(&mut stab_out, *s);
            stab_out.push('\t');
            stab_out.push_str(&rank_q[i].to_string());
            stab_out.push('\t');
            stab_out.push_str(&rank_s[i].to_string());
            stab_out.push('\n');
        }
    }

    let sign_path = args.output_dir.join("loo_sign_stability.tsv");
    let rank_path = args.output_dir.join("rank_stability.tsv");
    let stab_path = args.output_dir.join("stability_ranked.tsv");
    crate::io::atomic_write(&sign_path, sign_out.as_bytes())
        .with_context(|| "writing loo_sign_stability.tsv")?;
    crate::io::atomic_write(&rank_path, rank_out.as_bytes())
        .with_context(|| "writing rank_stability.tsv")?;
    crate::io::atomic_write(&stab_path, stab_out.as_bytes())
        .with_context(|| "writing stability_ranked.tsv")?;
    eprintln!(
        "robustness: wrote {:?}, {:?}, and {:?}",
        sign_path, rank_path, stab_path,
    );

    let finished_at = SystemTime::now();
    let mut labeled: Vec<(String, PathBuf)> = vec![("baseline".into(), args.baseline.clone())];
    for (i, p) in loo_files.iter().enumerate() {
        labeled.push((format!("loo_{i}"), PathBuf::from(p)));
    }
    let labeled_refs: Vec<(&str, &Path)> = labeled
        .iter()
        .map(|(l, p)| (l.as_str(), p.as_path()))
        .collect();
    let inputs_sha256 = hash_labeled_inputs(&labeled_refs)?;
    let outputs = [sign_path.clone(), rank_path.clone(), stab_path.clone()];
    let sidecar = sidecar_path_for(&sign_path);
    let mut extras = JsonMap::new();
    extras.insert(
        "report_thresholds".into(),
        json!({
            "q_strict": args.report_q_strict,
            "q_relaxed": args.report_q_relaxed,
        }),
    );
    write_run_sidecar(
        &sidecar,
        "robustness",
        json!({
            "baseline": args.baseline.display().to_string(),
            "loo": args.loo,
            "output-dir": args.output_dir.display().to_string(),
            "top-k": args.top_k,
            "report-q-strict": args.report_q_strict,
            "report-q-relaxed": args.report_q_relaxed,
        }),
        &inputs_sha256,
        &outputs,
        started_at,
        finished_at,
        Some(extras),
    )?;
    eprintln!("robustness: sidecar={}", sidecar.display());
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
    let c_gene = headers.iter().position(|h| h == "gene_symbol");

    let mut out: BTreeMap<String, HashMap<String, DeLite>> = BTreeMap::new();
    for rec in rdr.records() {
        let row = rec?;
        let cmp = row[c_cmp].to_string();
        let key = format!("{}::{}", row[c_panel].trim(), row[c_assay].trim());
        let gene = c_gene
            .map(|i| row[i].trim().to_string())
            .unwrap_or_default();
        out.entry(cmp).or_default().insert(
            key,
            DeLite {
                gene_symbol: gene,
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
    Ok(Some(
        s.trim()
            .parse()
            .with_context(|| format!("parse f64 {:?}", s))?,
    ))
}

fn top_k_set(rows: &HashMap<String, DeLite>, k: usize) -> BTreeSet<String> {
    let mut v: Vec<(&String, f64)> = rows
        .iter()
        // Non-finite |mean_diff| (NaN/inf) is not a meaningful "top changed
        // feature": with `total_cmp`, NaN sorts ahead of every finite value in
        // descending order and would hijack the top-k set. Exclude it from the
        // ranking entirely (strict ranking semantics).
        .filter_map(|(k, r)| r.mean_diff.filter(|d| d.is_finite()).map(|d| (k, d.abs())))
        .collect();
    // Sort by |mean_diff| descending. The collection comes from a HashMap, so
    // iteration order is nondeterministic; without a tie-break, features with
    // equal magnitude at the k-boundary would be selected arbitrarily, making
    // the Jaccard/overlap output non-reproducible. The secondary key (feature
    // id, ascending lexicographic) makes the ordering a deterministic total
    // order. `total_cmp` keeps the primary comparison a total order; all
    // values here are finite after the filter above.
    v.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(b.0)));
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

#[cfg(test)]
mod tests {
    use super::*;

    fn lite(diff: Option<f64>) -> DeLite {
        DeLite {
            gene_symbol: String::new(),
            mean_diff: diff,
            bh_q: Some(0.01),
        }
    }

    #[test]
    fn top_k_set_excludes_non_finite_and_breaks_ties_by_id() {
        let mut rows = HashMap::new();
        rows.insert("ZZZ".to_string(), lite(Some(5.0)));
        rows.insert("AAA".to_string(), lite(Some(2.0))); // tie with KKK
        rows.insert("KKK".to_string(), lite(Some(2.0))); // tie with AAA
        rows.insert("NAN".to_string(), lite(Some(f64::NAN))); // must be excluded
        rows.insert("INF".to_string(), lite(Some(f64::INFINITY))); // must be excluded
        rows.insert("NONE".to_string(), lite(None)); // already excluded

        // Top-2 by |mean_diff|: ZZZ (5.0) then the tie {AAA,KKK} resolves to AAA.
        let got = top_k_set(&rows, 2);
        let expected: BTreeSet<String> = ["ZZZ", "AAA"].iter().map(|s| s.to_string()).collect();
        assert_eq!(
            got, expected,
            "NaN/inf must not enter top-k; ties resolve by ascending id"
        );
    }
}
