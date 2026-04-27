//! Compositional transforms for `atman decompose ica --transform ...`.
//!
//! Atman's canonical abundance is already on a log scale (`Log2Npx`,
//! `Log2Intensity`, `Log2Diann_pg_quantity`, etc.), so the
//! log-ratio transforms below reduce to linear operations on the
//! log values:
//!
//! - **CLR** (centered log-ratio): `x_i' = log x_i − mean_j(log x_j)`
//!   over the sample's proteins. On already-log data, this is
//!   per-sample mean-subtraction.
//! - **ALR** (additive log-ratio): `x_i' = log x_i − log x_ref`.
//!   Requires a reference protein index. On already-log data, this
//!   is per-sample subtraction of the reference column.
//! - **ILR** (isometric log-ratio): an orthonormal Helmert projection
//!   of the CLR vector onto a `(p−1)`-dimensional unconstrained
//!   subspace. Numerically better-conditioned than CLR for large K.
//! - **RatioAnchor**: divide every protein by a reference protein,
//!   then log. On already-log input this is algebraically identical
//!   to ALR; kept as a named alias because the CSF-manuscript
//!   "albumin normalization" pipeline uses this vocabulary.
//!
//! Zero-handling (`MinProb` / `KeepNa` in the feature request) is
//! only meaningful on linear-intensity input. Atman's load path
//! already emits finite values for every canonical record (QC-masked
//! cells are dropped upstream or mean-imputed by `decompose ica`),
//! so v1 implements the transforms without a separate zero-handling
//! knob and leaves that work to a follow-on once linear-scale input
//! routes exist.

/// All transforms are pure functions on borrowed slices; no IO, no
/// RNG.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Transform {
    None,
    Clr,
    /// Additive log-ratio against the protein at `reference_index`.
    Alr {
        reference_index: usize,
    },
    /// Isometric log-ratio via the Helmert orthonormal basis.
    Ilr,
    /// Per-sample ratio-against-anchor, same math as [`Transform::Alr`]
    /// on log-scale input but preserved as a distinct tag because
    /// users may reach for one term or the other depending on
    /// domain convention.
    RatioAnchor {
        reference_index: usize,
    },
}

impl Transform {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Clr => "clr",
            Self::Alr { .. } => "alr",
            Self::Ilr => "ilr",
            Self::RatioAnchor { .. } => "ratio-anchor",
        }
    }
}

/// Apply a compositional transform to a sample × protein matrix
/// (`data[i][j]` is sample `i`, protein `j`) whose values are already
/// on a log scale.
///
/// - [`Transform::None`] clones the input.
/// - [`Transform::Clr`] returns the same shape with per-row mean
///   subtracted.
/// - [`Transform::Alr`] / [`Transform::RatioAnchor`] return the same
///   shape with the reference column subtracted from every column
///   in each row. The reference column itself becomes identically
///   zero.
/// - [`Transform::Ilr`] returns an `n × (p − 1)` matrix.
///
/// Returns `Err` when the transform is malformed (e.g. ALR reference
/// out of bounds) or when the input is rectangularly invalid.
pub fn apply_transform(data: &[Vec<f64>], transform: Transform) -> Result<Vec<Vec<f64>>, String> {
    if data.is_empty() {
        return Ok(Vec::new());
    }
    let p = data[0].len();
    if data.iter().any(|row| row.len() != p) {
        return Err("non-rectangular input matrix".to_string());
    }
    match transform {
        Transform::None => Ok(data.to_vec()),
        Transform::Clr => Ok(clr(data)),
        Transform::Alr { reference_index } | Transform::RatioAnchor { reference_index } => {
            if reference_index >= p {
                return Err(format!(
                    "alr reference_index {reference_index} out of bounds (p={p})"
                ));
            }
            Ok(alr(data, reference_index))
        }
        Transform::Ilr => {
            if p < 2 {
                return Err(format!("ilr requires p ≥ 2, got p={p}"));
            }
            Ok(ilr(data))
        }
    }
}

fn clr(data: &[Vec<f64>]) -> Vec<Vec<f64>> {
    data.iter()
        .map(|row| {
            let mean: f64 = row.iter().sum::<f64>() / row.len() as f64;
            row.iter().map(|v| v - mean).collect()
        })
        .collect()
}

fn alr(data: &[Vec<f64>], reference_index: usize) -> Vec<Vec<f64>> {
    data.iter()
        .map(|row| {
            let r = row[reference_index];
            row.iter().map(|v| v - r).collect()
        })
        .collect()
}

/// ILR using a Helmert orthonormal basis. Given a CLR vector `z` of
/// length `p`, the Helmert matrix `V` of shape `(p−1) × p` has rows:
///
/// ```text
/// V[k][j] =   sqrt(k / (k + 1))                         for j == k
///         = -1 / sqrt(k * (k + 1))                      for 0 <= j < k
///         =  0                                          for j > k
/// ```
///
/// `V V^T = I_{p−1}` and `y = V z` is the ILR coordinate. Since
/// `CLR` removes one degree of freedom (row-sum = 0), the ILR is
/// full-rank over the `(p−1)`-dimensional simplex tangent space.
fn ilr(data: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let p = data[0].len();
    let z = clr(data);
    let mut basis: Vec<Vec<f64>> = vec![vec![0.0; p]; p - 1];
    for (k, row) in basis.iter_mut().enumerate() {
        let kp1 = (k + 1) as f64;
        let diag = (kp1 / (kp1 + 1.0)).sqrt();
        let off = -1.0 / (kp1 * (kp1 + 1.0)).sqrt();
        for slot in row.iter_mut().take(k + 1) {
            *slot = off;
        }
        row[k + 1] = diag;
    }
    z.iter()
        .map(|row| {
            basis
                .iter()
                .map(|basis_row| basis_row.iter().zip(row.iter()).map(|(b, v)| b * v).sum())
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transform_none_is_identity() {
        let data = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let out = apply_transform(&data, Transform::None).unwrap();
        assert_eq!(out, data);
    }

    #[test]
    fn clr_rows_sum_to_zero() {
        let data = vec![
            vec![1.0, 3.0, 5.0],
            vec![-2.0, 0.0, 2.0],
            vec![10.0, 10.0, 10.0],
        ];
        let out = apply_transform(&data, Transform::Clr).unwrap();
        for row in &out {
            let s: f64 = row.iter().sum();
            assert!(s.abs() < 1e-12, "CLR row should sum to 0, got {s}");
        }
    }

    #[test]
    fn clr_constant_row_becomes_zero() {
        let data = vec![vec![7.0, 7.0, 7.0, 7.0]];
        let out = apply_transform(&data, Transform::Clr).unwrap();
        for v in &out[0] {
            assert!(v.abs() < 1e-12);
        }
    }

    #[test]
    fn clr_invariance_protein_at_mean_maps_to_zero() {
        // If a column is always at the per-row mean, CLR sends it
        // to zero on every row. (This is the spec's "invariance
        // check" for compositional data.)
        let data = vec![
            vec![1.0, 2.0, 3.0],   // mean 2.0 → col 1 at mean
            vec![0.0, 5.0, 10.0],  // mean 5.0 → col 1 at mean
            vec![-4.0, -1.0, 2.0], // mean -1.0 → col 1 at mean
        ];
        let out = apply_transform(&data, Transform::Clr).unwrap();
        for row in &out {
            assert!(row[1].abs() < 1e-12, "mean-column should CLR to 0");
        }
    }

    #[test]
    fn alr_reference_column_becomes_zero() {
        let data = vec![vec![1.0, 5.0, 2.0], vec![3.0, 4.0, 7.0]];
        let out = apply_transform(&data, Transform::Alr { reference_index: 1 }).unwrap();
        for row in &out {
            assert!(row[1].abs() < 1e-12);
        }
    }

    #[test]
    fn alr_matches_hand_computation() {
        let data = vec![vec![1.0, 5.0, 2.0], vec![3.0, 4.0, 7.0]];
        let out = apply_transform(&data, Transform::Alr { reference_index: 0 }).unwrap();
        assert_eq!(out[0], vec![0.0, 4.0, 1.0]);
        assert_eq!(out[1], vec![0.0, 1.0, 4.0]);
    }

    #[test]
    fn ratio_anchor_equals_alr_on_log_scale_input() {
        let data = vec![vec![1.0, 5.0, 2.0], vec![3.0, 4.0, 7.0]];
        let a = apply_transform(&data, Transform::Alr { reference_index: 2 }).unwrap();
        let b = apply_transform(&data, Transform::RatioAnchor { reference_index: 2 }).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn ilr_output_has_p_minus_one_columns() {
        let data = vec![vec![1.0, 2.0, 3.0, 4.0], vec![5.0, 6.0, 7.0, 8.0]];
        let out = apply_transform(&data, Transform::Ilr).unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].len(), 3);
    }

    #[test]
    fn ilr_is_orthonormal_preserves_norm_of_clr() {
        // For any CLR-centered row z, ||y||² = ||z||² under ILR
        // orthonormal basis.
        let data = vec![vec![1.0, 3.0, -2.0, 0.0], vec![5.0, -1.0, 2.0, -6.0]];
        let clr_rows = apply_transform(&data, Transform::Clr).unwrap();
        let ilr_rows = apply_transform(&data, Transform::Ilr).unwrap();
        for (c_row, i_row) in clr_rows.iter().zip(ilr_rows.iter()) {
            let c_norm2: f64 = c_row.iter().map(|v| v * v).sum();
            let i_norm2: f64 = i_row.iter().map(|v| v * v).sum();
            assert!(
                (c_norm2 - i_norm2).abs() < 1e-9,
                "ILR should preserve CLR norm; clr²={c_norm2}, ilr²={i_norm2}"
            );
        }
    }

    #[test]
    fn alr_refuses_out_of_bounds_reference() {
        let data = vec![vec![1.0, 2.0]];
        let err = apply_transform(&data, Transform::Alr { reference_index: 5 }).unwrap_err();
        assert!(err.contains("out of bounds"), "unexpected: {err}");
    }

    #[test]
    fn ilr_refuses_p_below_two() {
        let data = vec![vec![1.0]];
        let err = apply_transform(&data, Transform::Ilr).unwrap_err();
        assert!(err.contains("p ≥ 2"), "unexpected: {err}");
    }

    #[test]
    fn transform_is_deterministic() {
        let data = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let a = apply_transform(&data, Transform::Clr).unwrap();
        let b = apply_transform(&data, Transform::Clr).unwrap();
        assert_eq!(a, b);
    }
}
