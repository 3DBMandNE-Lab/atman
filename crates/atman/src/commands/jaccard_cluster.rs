//! Shared helpers for absence-pattern clustering commands (`recover-plex`,
//! `absence-topology`).
//!
//! Both commands build per-entity absence bitsets, compute pairwise Jaccard
//! distance, cluster by connected components (union-find), summarise cluster
//! sizes (median), and format floating-point quality cells uniformly. These
//! routines were previously copy-pasted into each command; they live here as a
//! single source of truth. Pure refactor — behaviour is byte-identical to the
//! prior per-file copies.
//!
//! Note: `classify` is intentionally NOT shared. The plex and protein variants
//! have different signatures and emit different diagnostic strings, so each
//! command keeps its own.

use std::path::{Path, PathBuf};

/// Union-find `find` with path compression. The root is the minimum index in
/// the component (enforced by [`union`]), which keeps clustering deterministic.
pub fn find(parent: &mut [usize], x: usize) -> usize {
    let mut r = x;
    while parent[r] != r {
        r = parent[r];
    }
    let mut cur = x;
    while parent[cur] != r {
        let nxt = parent[cur];
        parent[cur] = r;
        cur = nxt;
    }
    r
}

/// Union-find `union`. Deterministic: the smaller index becomes the root.
pub fn union(parent: &mut [usize], a: usize, b: usize) {
    let ra = find(parent, a);
    let rb = find(parent, b);
    if ra != rb {
        if ra < rb {
            parent[rb] = ra;
        } else {
            parent[ra] = rb;
        }
    }
}

/// Jaccard distance between two equal-length bitsets: `1 - |A∩B| / |A∪B|`.
/// Two all-zero bitsets (empty absence sets) have distance 0.
pub fn jaccard_distance(a: &[u64], b: &[u64]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    let mut inter = 0u64;
    let mut uni = 0u64;
    for (x, y) in a.iter().zip(b.iter()) {
        inter += (x & y).count_ones() as u64;
        uni += (x | y).count_ones() as u64;
    }
    if uni == 0 {
        0.0
    } else {
        1.0 - (inter as f64) / (uni as f64)
    }
}

/// Median of a slice; `NaN` for an empty slice. `NaN` entries sort to the end
/// via `partial_cmp` fallback (callers pass finite cluster sizes).
pub fn median_f(v: &[f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        0.5 * (s[n / 2 - 1] + s[n / 2])
    }
}

/// Format a quality-cell float: empty string for `NaN`, `inf`/`-inf` for
/// infinities, otherwise six-decimal fixed point.
pub fn format_f(x: f64) -> String {
    if x.is_nan() {
        String::new()
    } else if x.is_infinite() {
        if x > 0.0 {
            "inf".to_string()
        } else {
            "-inf".to_string()
        }
    } else {
        format!("{:.6}", x)
    }
}

/// Derive the `<stem>_quality.<ext>` sidecar path next to `out`, falling back to
/// `default_stem` / `tsv` when `out` has no stem / extension.
pub fn quality_path_for(out: &Path, default_stem: &str) -> PathBuf {
    let stem = out
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| default_stem.to_string());
    let ext = out
        .extension()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "tsv".to_string());
    out.with_file_name(format!("{}_quality.{}", stem, ext))
}
