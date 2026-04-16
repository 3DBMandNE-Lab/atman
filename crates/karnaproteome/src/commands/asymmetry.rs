use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;
use csv::ReaderBuilder;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::PathBuf;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Path to de_results.tsv produced by `karnaproteome de`.
    #[arg(long)]
    de_results: PathBuf,

    /// Output TSV path for asymmetry metrics.
    #[arg(long)]
    output: PathBuf,

    /// Comma-separated comparison pairs in `A:B` form.
    /// Example: "PT2-PT1:PR2-PR1,PT2-PR2:PT1-PR1"
    #[arg(long)]
    pairs: String,
}

#[derive(Debug, Clone)]
struct DeLite {
    mean_diff: Option<f64>,
    bh_q: Option<f64>,
}

pub fn run(args: Args) -> Result<()> {
    let pairs = parse_pairs(&args.pairs)?;
    let by_comp = read_de_results(&args.de_results)?;

    let mut out = String::from(
        "comparison_a\tcomparison_b\tn_common\t\
         n_sig05_a\tn_sig05_b\tsig05_ratio_a_over_b\t\
         n_sig10_a\tn_sig10_b\tsig10_ratio_a_over_b\t\
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

        let mut sig05_a = 0usize;
        let mut sig05_b = 0usize;
        let mut sig10_a = 0usize;
        let mut sig10_b = 0usize;
        let mut abs_a: Vec<f64> = Vec::new();
        let mut abs_b: Vec<f64> = Vec::new();
        let mut sign_pairs = 0usize;
        let mut sign_same = 0usize;

        for k in common.iter().copied() {
            let ra = ma.get(k).expect("key from intersection must exist");
            let rb = mb.get(k).expect("key from intersection must exist");

            if let Some(q) = ra.bh_q {
                if q < 0.05 {
                    sig05_a += 1;
                }
                if q < 0.10 {
                    sig10_a += 1;
                }
            }
            if let Some(q) = rb.bh_q {
                if q < 0.05 {
                    sig05_b += 1;
                }
                if q < 0.10 {
                    sig10_b += 1;
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
        let sig05_ratio = ratio(sig05_a, sig05_b);
        let sig10_ratio = ratio(sig10_a, sig10_b);
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
        out.push_str(&sig05_a.to_string());
        out.push('\t');
        out.push_str(&sig05_b.to_string());
        out.push('\t');
        push_opt_f64(&mut out, sig05_ratio);
        out.push('\t');
        out.push_str(&sig10_a.to_string());
        out.push('\t');
        out.push_str(&sig10_b.to_string());
        out.push('\t');
        push_opt_f64(&mut out, sig10_ratio);
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

    std::fs::write(&args.output, out.as_bytes())
        .with_context(|| format!("writing {:?}", args.output))?;
    eprintln!("asymmetry: wrote {:?}", args.output);
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
    let v: f64 = s.trim().parse().with_context(|| format!("parsing f64 {:?}", s))?;
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

