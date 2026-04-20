//! Preflight safety gates that run before any DE test path opens
//! the measurement table. Gates refuse work that the default tests
//! would silently mishandle; opting past them requires an explicit
//! CLI flag.

use anyhow::{anyhow, Context, Result};
use std::path::Path;

/// Refuse to run DE when `below_lod=1` rows exceed `max_fraction` of
/// non-dropped measurements, unless the caller explicitly opts in via
/// `--allow-censored`. Reads only the two columns the gate needs so
/// the full record parse is deferred until the chosen test path
/// loads it.
pub(super) fn check_below_lod_gate(
    path: &Path,
    max_fraction: f64,
    allow_censored: bool,
) -> Result<()> {
    if !(0.0..=1.0).contains(&max_fraction) || !max_fraction.is_finite() {
        anyhow::bail!("--max-below-lod-fraction must be in [0.0, 1.0]");
    }
    if !path.exists() {
        return Ok(());
    }
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader
        .headers()
        .with_context(|| format!("reading headers from {:?}", path))?
        .clone();
    let col_idx: std::collections::HashMap<String, usize> = headers
        .iter()
        .enumerate()
        .map(|(i, h)| (h.to_string(), i))
        .collect();
    let drop_col = col_idx
        .get("dropped_by_qc")
        .copied()
        .ok_or_else(|| anyhow!("missing column dropped_by_qc in {:?}", path))?;
    let below_col = col_idx
        .get("below_lod")
        .copied()
        .ok_or_else(|| anyhow!("missing column below_lod in {:?}", path))?;
    let mut total_usable: u64 = 0;
    let mut below: u64 = 0;
    for record in reader.records() {
        let row = record.with_context(|| format!("reading record from {:?}", path))?;
        let dropped: bool = row
            .get(drop_col)
            .and_then(|v| v.parse::<u8>().ok())
            .map(|v| v != 0)
            .unwrap_or(false);
        if dropped {
            continue;
        }
        total_usable += 1;
        let is_below: bool = row
            .get(below_col)
            .and_then(|v| v.parse::<u8>().ok())
            .map(|v| v != 0)
            .unwrap_or(false);
        if is_below {
            below += 1;
        }
    }
    if total_usable == 0 {
        return Ok(());
    }
    let fraction = below as f64 / total_usable as f64;
    if fraction > max_fraction && !allow_censored {
        anyhow::bail!(
            "below-LOD fraction {fraction:.3} exceeds --max-below-lod-fraction \
{max_fraction:.3} ({below} of {total_usable} non-dropped measurements). \
Atman's default DE tests treat missing as MCAR; left-censored data at this \
rate will bias mean_diff estimates toward zero and distort t-statistics. \
Either (a) impute or filter upstream, (b) raise --max-below-lod-fraction \
after confirming the test is appropriate, or (c) pass --allow-censored to \
bypass this gate."
        );
    }
    if fraction > max_fraction {
        eprintln!(
            "de: WARNING below-LOD fraction {fraction:.3} exceeds \
{max_fraction:.3} (allow-censored set)."
        );
    }
    Ok(())
}
