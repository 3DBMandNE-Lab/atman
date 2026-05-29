use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;
use csv::ReaderBuilder;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::PathBuf;
use std::time::SystemTime;

use crate::io::{hash_labeled_inputs, sidecar_path_for, write_run_sidecar};
use serde_json::Map as JsonMap;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Path to de_results.tsv produced by `atman de`.
    #[arg(long)]
    de_results: PathBuf,

    /// Output TSV path for asymmetry metrics.
    #[arg(long)]
    output: PathBuf,

    /// Comma-separated comparison pairs in `A:B` form.
    /// Example: "PT2-PT1:PR2-PR1,PT2-PR2:PT1-PR1"
    #[arg(long)]
    pairs: String,

    /// Strict BH-q threshold for the `n_sig_strict_*` columns and
    /// `sig_strict_ratio_a_over_b`. Reporting policy only — no effect
    /// on which rows are read or written. Must satisfy
    /// `0 < strict < relaxed < 1`.
    #[arg(long, default_value_t = 0.05)]
    report_q_strict: f64,

    /// Relaxed BH-q threshold for the `n_sig_relaxed_*` columns and
    /// `sig_relaxed_ratio_a_over_b`. Same scope as `--report-q-strict`.
    #[arg(long, default_value_t = 0.10)]
    report_q_relaxed: f64,
}

#[derive(Debug, Clone)]
struct DeLite {
    mean_diff: Option<f64>,
    bh_q: Option<f64>,
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
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
    let pairs = parse_pairs(&args.pairs)?;
    let by_comp = read_de_results(&args.de_results)?;

    let mut out = String::from(
        "comparison_a\tcomparison_b\tn_common\t\
         n_sig_strict_a\tn_sig_strict_b\tsig_strict_ratio_a_over_b\t\
         n_sig_relaxed_a\tn_sig_relaxed_b\tsig_relaxed_ratio_a_over_b\t\
         median_abs_effect_a\tmedian_abs_effect_b\tabs_effect_ratio_a_over_b\t\
         sign_concordance_rate\n",
    );

    for (a, b) in pairs {
        let ma = by_comp
            .get(&a)
            .with_context(|| format!("comparison {:?} not found in {:?}", a, args.de_results))?;
        let mb = by_comp
            .get(&b)
            .with_context(|| format!("comparison {:?} not found in {:?}", b, args.de_results))?;

        let keys_a: BTreeSet<&str> = ma.keys().map(|k| k.as_str()).collect();
        let keys_b: BTreeSet<&str> = mb.keys().map(|k| k.as_str()).collect();
        let common: Vec<&str> = keys_a.intersection(&keys_b).copied().collect();

        let mut sig_strict_a = 0usize;
        let mut sig_strict_b = 0usize;
        let mut sig_relaxed_a = 0usize;
        let mut sig_relaxed_b = 0usize;
        let mut abs_a: Vec<f64> = Vec::new();
        let mut abs_b: Vec<f64> = Vec::new();
        let mut sign_pairs = 0usize;
        let mut sign_same = 0usize;

        for k in common.iter().copied() {
            let ra = ma.get(k).expect("key from intersection must exist");
            let rb = mb.get(k).expect("key from intersection must exist");

            if let Some(q) = ra.bh_q {
                if q < args.report_q_strict {
                    sig_strict_a += 1;
                }
                if q < args.report_q_relaxed {
                    sig_relaxed_a += 1;
                }
            }
            if let Some(q) = rb.bh_q {
                if q < args.report_q_strict {
                    sig_strict_b += 1;
                }
                if q < args.report_q_relaxed {
                    sig_relaxed_b += 1;
                }
            }

            if let Some(d) = ra.mean_diff {
                abs_a.push(d.abs());
            }
            if let Some(d) = rb.mean_diff {
                abs_b.push(d.abs());
            }

            if let (Some(da), Some(db)) = (ra.mean_diff, rb.mean_diff) {
                if da != 0.0 && db != 0.0 {
                    sign_pairs += 1;
                    if da.signum() == db.signum() {
                        sign_same += 1;
                    }
                }
            }
        }

        let n_common = common.len();
        let sig_strict_ratio = ratio(sig_strict_a, sig_strict_b);
        let sig_relaxed_ratio = ratio(sig_relaxed_a, sig_relaxed_b);
        let med_abs_a = median(&mut abs_a);
        let med_abs_b = median(&mut abs_b);
        let abs_ratio = opt_ratio(med_abs_a, med_abs_b);
        let sign_rate = if sign_pairs == 0 {
            None
        } else {
            Some(sign_same as f64 / sign_pairs as f64)
        };

        out.push_str(&a);
        out.push('\t');
        out.push_str(&b);
        out.push('\t');
        out.push_str(&n_common.to_string());
        out.push('\t');
        out.push_str(&sig_strict_a.to_string());
        out.push('\t');
        out.push_str(&sig_strict_b.to_string());
        out.push('\t');
        push_opt_f64(&mut out, sig_strict_ratio);
        out.push('\t');
        out.push_str(&sig_relaxed_a.to_string());
        out.push('\t');
        out.push_str(&sig_relaxed_b.to_string());
        out.push('\t');
        push_opt_f64(&mut out, sig_relaxed_ratio);
        out.push('\t');
        push_opt_f64(&mut out, med_abs_a);
        out.push('\t');
        push_opt_f64(&mut out, med_abs_b);
        out.push('\t');
        push_opt_f64(&mut out, abs_ratio);
        out.push('\t');
        push_opt_f64(&mut out, sign_rate);
        out.push('\n');
    }

    crate::io::atomic_write(&args.output, out.as_bytes())
        .with_context(|| format!("writing {:?}", args.output))?;
    eprintln!("asymmetry: wrote {:?}", args.output);

    let finished_at = SystemTime::now();
    let inputs_sha256 = hash_labeled_inputs(&[("de_results", args.de_results.as_path())])?;
    let sidecar = sidecar_path_for(&args.output);
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
        "asymmetry",
        json!({
            "de-results": args.de_results.display().to_string(),
            "output": args.output.display().to_string(),
            "pairs": args.pairs,
            "report-q-strict": args.report_q_strict,
            "report-q-relaxed": args.report_q_relaxed,
        }),
        &inputs_sha256,
        std::slice::from_ref(&args.output),
        started_at,
        finished_at,
        Some(extras),
    )?;
    eprintln!("asymmetry: sidecar={}", sidecar.display());
    Ok(())
}

fn parse_pairs(raw: &str) -> Result<Vec<(String, String)>> {
    let mut pairs = Vec::new();
    for tok in raw.split(',') {
        let t = tok.trim();
        if t.is_empty() {
            bail!("invalid empty pair token in --pairs {:?}", raw);
        }
        let mut it = t.split(':').map(str::trim);
        let a = it.next().unwrap_or_default();
        let b = it.next().unwrap_or_default();
        if a.is_empty() || b.is_empty() || it.next().is_some() {
            bail!("invalid pair {:?}; expected A:B", t);
        }
        pairs.push((a.to_string(), b.to_string()));
    }
    if pairs.is_empty() {
        bail!("no pairs provided");
    }
    Ok(pairs)
}

fn read_de_results(path: &PathBuf) -> Result<BTreeMap<String, HashMap<String, DeLite>>> {
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = rdr
        .headers()
        .with_context(|| format!("reading headers from {:?}", path))?
        .clone();
    let need = |name: &str| -> Result<usize> {
        headers
            .iter()
            .position(|h| h == name)
            .with_context(|| format!("missing column {:?} in {:?}", name, path))
    };
    let c_panel = need("panel")?;
    let c_assay = need("assay_id")?;
    let c_cmp = need("comparison")?;
    let c_diff = need("mean_diff")?;
    let c_q = need("bh_q")?;

    let mut out: BTreeMap<String, HashMap<String, DeLite>> = BTreeMap::new();
    for rec in rdr.records() {
        let row = rec.with_context(|| format!("reading row from {:?}", path))?;
        let cmp = row[c_cmp].to_string();
        let key = format!("{}::{}", row[c_panel].trim(), row[c_assay].trim());
        let mean_diff = parse_opt_f64(&row[c_diff])?;
        let bh_q = parse_opt_f64(&row[c_q])?;
        out.entry(cmp)
            .or_default()
            .insert(key, DeLite { mean_diff, bh_q });
    }
    Ok(out)
}

fn parse_opt_f64(s: &str) -> Result<Option<f64>> {
    if s.trim().is_empty() {
        return Ok(None);
    }
    let v: f64 = s
        .trim()
        .parse()
        .with_context(|| format!("parsing f64 {:?}", s))?;
    Ok(Some(v))
}

fn ratio(a: usize, b: usize) -> Option<f64> {
    if b == 0 {
        None
    } else {
        Some(a as f64 / b as f64)
    }
}

fn opt_ratio(a: Option<f64>, b: Option<f64>) -> Option<f64> {
    match (a, b) {
        (Some(x), Some(y)) if y != 0.0 => Some(x / y),
        _ => None,
    }
}

fn median(v: &mut [f64]) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = v.len();
    if n % 2 == 1 {
        Some(v[n / 2])
    } else {
        Some((v[n / 2 - 1] + v[n / 2]) / 2.0)
    }
}

fn push_opt_f64(buf: &mut String, v: Option<f64>) {
    if let Some(x) = v {
        buf.push_str(&x.to_string());
    }
}

#[cfg(test)]
mod tests {
    use super::parse_pairs;

    #[test]
    fn parse_pairs_ok() {
        let got = parse_pairs("PT2-PT1:PR2-PR1,PT2-PR2:PT1-PR1").unwrap();
        assert_eq!(
            got,
            vec![
                ("PT2-PT1".to_string(), "PR2-PR1".to_string()),
                ("PT2-PR2".to_string(), "PT1-PR1".to_string())
            ]
        );
    }

    #[test]
    fn parse_pairs_bad() {
        let err = parse_pairs("PT2-PT1").unwrap_err().to_string();
        assert!(err.contains("expected A:B"));
    }
}
