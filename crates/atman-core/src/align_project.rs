//! Project a new cohort's subject × protein matrix onto a trained
//! cross-cohort atlas of archetype loading vectors.
//!
//! The atlas is built upstream (via `atman align programs`) from N
//! per-cohort ICA loading TSVs: each multi-cohort archetype is
//! represented by the mean loading vector across its member programs,
//! aligned to the intersection of the member cohorts' protein
//! universes. Projection restricts an incoming cohort's abundance
//! matrix to the atlas universe and solves for per-subject activations
//! `β` via:
//!
//! - **Least squares:** `β = (AᵀA)⁻¹ Aᵀ x`.
//! - **Ridge:** `β = (AᵀA + λI)⁻¹ Aᵀ x` for `λ > 0` — the default.
//!
//! Both solves go through Cholesky on the normal equations; atlases
//! with `k ≤ p` (archetypes much fewer than proteins) are well-
//! conditioned for this path. NNLS is a documented follow-on.
//!
//! Transforms are applied by the caller before reaching this module,
//! so the atlas loading vectors and the incoming abundance rows must
//! already live in the same coordinate system (e.g. both CLR-centered).

#[derive(Debug, Clone)]
pub struct Atlas {
    /// Canonical ordering of protein labels; atlas loading vectors
    /// index into this.
    pub protein_labels: Vec<String>,
    /// Integer archetype id strings (e.g. `"A0001"`), in atlas order.
    pub archetype_ids: Vec<String>,
    /// `[archetype][protein]` mean loading on `protein_labels`.
    pub loadings: Vec<Vec<f64>>,
    /// Per-archetype member program count at atlas-build time.
    pub n_members: Vec<usize>,
    /// Per-archetype cohort count at atlas-build time.
    pub n_cohorts: Vec<usize>,
}

#[derive(Debug, Clone, Copy)]
pub enum ProjectionMethod {
    LeastSquares,
    Ridge(f64),
}

#[derive(Debug, Clone)]
pub struct SubjectQc {
    pub subject_id: String,
    pub residual_norm: f64,
    pub coverage_fraction: f64,
    pub n_present: usize,
    pub n_missing: usize,
}

#[derive(Debug, Clone)]
pub struct ProjectionResult {
    pub subject_ids: Vec<String>,
    pub archetype_ids: Vec<String>,
    /// `[subject][archetype]`.
    pub activations: Vec<Vec<f64>>,
    pub qc: Vec<SubjectQc>,
    /// Proteins that were in the atlas but missing from the incoming
    /// cohort (filled with 0 for the projection solve); emitted for
    /// the sidecar so the user can see what was imputed.
    pub atlas_proteins_missing_in_cohort: Vec<String>,
}

fn cholesky_lower(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    if n == 0 {
        return Some(Vec::new());
    }
    if a.iter().any(|row| row.len() != n) {
        return None;
    }
    let mut l = vec![vec![0.0_f64; n]; n];
    for j in 0..n {
        let diag_sum: f64 = l[j].iter().take(j).map(|v| v * v).sum();
        let diag = a[j][j] - diag_sum;
        if !diag.is_finite() || diag <= 0.0 {
            return None;
        }
        l[j][j] = diag.sqrt();
        for i in (j + 1)..n {
            let off_sum: f64 = l[i]
                .iter()
                .zip(l[j].iter())
                .take(j)
                .map(|(li, lj)| li * lj)
                .sum();
            l[i][j] = (a[i][j] - off_sum) / l[j][j];
        }
    }
    Some(l)
}

fn solve_cholesky(l: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
    let p = rhs.len();
    let mut z = vec![0.0; p];
    for i in 0..p {
        let mut sum = rhs[i];
        for (k, zk) in z.iter().enumerate().take(i) {
            sum -= l[i][k] * zk;
        }
        z[i] = sum / l[i][i];
    }
    let mut out = vec![0.0; p];
    for i in (0..p).rev() {
        let mut sum = z[i];
        for k in (i + 1)..p {
            sum -= l[k][i] * out[k];
        }
        out[i] = sum / l[i][i];
    }
    out
}

/// Project `abundance` (subject × protein on `cohort_protein_labels`)
/// onto `atlas`. The intersection of protein labels is used; atlas
/// proteins absent from the cohort are filled with 0 (caller's
/// responsibility to warn — `atlas_proteins_missing_in_cohort` is
/// returned for the sidecar).
pub fn project(
    atlas: &Atlas,
    cohort_protein_labels: &[String],
    subject_ids: &[String],
    abundance: &[Vec<f64>],
    method: ProjectionMethod,
) -> Result<ProjectionResult, String> {
    if atlas.protein_labels.is_empty() {
        return Err("atlas has empty protein universe".into());
    }
    if atlas.loadings.is_empty() {
        return Err("atlas has zero archetypes to project onto".into());
    }
    if abundance.len() != subject_ids.len() {
        return Err(format!(
            "subject_ids ({}) ≠ abundance rows ({})",
            subject_ids.len(),
            abundance.len()
        ));
    }
    let n_subjects = subject_ids.len();
    let k = atlas.loadings.len();
    let p_atlas = atlas.protein_labels.len();
    for (i, l) in atlas.loadings.iter().enumerate() {
        if l.len() != p_atlas {
            return Err(format!(
                "atlas archetype {} has loading length {} ≠ universe {}",
                i,
                l.len(),
                p_atlas
            ));
        }
    }
    for (i, row) in abundance.iter().enumerate() {
        if row.len() != cohort_protein_labels.len() {
            return Err(format!(
                "cohort subject {} has {} values but cohort universe is {}",
                i,
                row.len(),
                cohort_protein_labels.len()
            ));
        }
    }

    // Map cohort labels → column index so we can pull values in atlas
    // order. Atlas proteins with no cohort match contribute zero to
    // the projection (and are reported separately for QC).
    use std::collections::HashMap;
    let cohort_idx: HashMap<&str, usize> = cohort_protein_labels
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let mut atlas_to_cohort: Vec<Option<usize>> = Vec::with_capacity(p_atlas);
    let mut missing: Vec<String> = Vec::new();
    for label in &atlas.protein_labels {
        match cohort_idx.get(label.as_str()) {
            Some(&j) => atlas_to_cohort.push(Some(j)),
            None => {
                atlas_to_cohort.push(None);
                missing.push(label.clone());
            }
        }
    }

    // Build A (p_atlas × k) as column-major inside a flat Vec<Vec<f64>>
    // indexed as `a[protein][archetype]`.
    let mut a_mat: Vec<Vec<f64>> = vec![vec![0.0_f64; k]; p_atlas];
    for (ai, loading) in atlas.loadings.iter().enumerate() {
        for (pi, &v) in loading.iter().enumerate() {
            a_mat[pi][ai] = v;
        }
    }
    // AᵀA: k × k.
    let mut ata = vec![vec![0.0_f64; k]; k];
    for pi in 0..p_atlas {
        for i in 0..k {
            for j in 0..k {
                ata[i][j] += a_mat[pi][i] * a_mat[pi][j];
            }
        }
    }
    // Ridge shift on the diagonal.
    let lambda = match method {
        ProjectionMethod::LeastSquares => 0.0,
        ProjectionMethod::Ridge(l) => {
            if !(l.is_finite() && l >= 0.0) {
                return Err(format!("ridge λ must be finite and non-negative; got {l}"));
            }
            l
        }
    };
    for i in 0..k {
        ata[i][i] += lambda;
    }
    let chol = cholesky_lower(&ata).ok_or_else(|| {
        "AᵀA not positive-definite; try --projection ridge with λ > 0".to_string()
    })?;

    let mut activations = vec![vec![0.0_f64; k]; n_subjects];
    let mut qc = Vec::with_capacity(n_subjects);
    for (si, subject_id) in subject_ids.iter().enumerate() {
        // Build x on the atlas universe (missing → 0).
        let mut x = vec![0.0_f64; p_atlas];
        let mut n_present = 0usize;
        let mut n_missing = 0usize;
        for (pi, mapped) in atlas_to_cohort.iter().enumerate() {
            match mapped {
                Some(j) => {
                    let v = abundance[si][*j];
                    if v.is_finite() {
                        x[pi] = v;
                        n_present += 1;
                    } else {
                        n_missing += 1;
                    }
                }
                None => n_missing += 1,
            }
        }
        let coverage = if p_atlas == 0 {
            0.0
        } else {
            n_present as f64 / p_atlas as f64
        };

        // AᵀX.
        let mut atx = vec![0.0_f64; k];
        for pi in 0..p_atlas {
            for i in 0..k {
                atx[i] += a_mat[pi][i] * x[pi];
            }
        }
        // Solve (AᵀA + λI) β = Aᵀ x.
        let beta = solve_cholesky(&chol, &atx);
        // Residual x − A β.
        let mut resid_sq = 0.0;
        for pi in 0..p_atlas {
            let mut pred = 0.0;
            for (i, &b) in beta.iter().enumerate() {
                pred += a_mat[pi][i] * b;
            }
            let r = x[pi] - pred;
            resid_sq += r * r;
        }
        activations[si] = beta;
        qc.push(SubjectQc {
            subject_id: subject_id.clone(),
            residual_norm: resid_sq.sqrt(),
            coverage_fraction: coverage,
            n_present,
            n_missing,
        });
    }

    Ok(ProjectionResult {
        subject_ids: subject_ids.to_vec(),
        archetype_ids: atlas.archetype_ids.clone(),
        activations,
        qc,
        atlas_proteins_missing_in_cohort: missing,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_atlas_and_cohort() -> (Atlas, Vec<String>, Vec<String>, Vec<Vec<f64>>) {
        // 2 archetypes on 4 proteins, 3 subjects. Archetype 1 loads
        // strongly on P1/P2, archetype 2 on P3/P4. We synthesize
        // subject data as known linear combinations then check the
        // projection recovers those coefficients to tight tolerance.
        let atlas = Atlas {
            protein_labels: vec!["P1", "P2", "P3", "P4"]
                .into_iter()
                .map(String::from)
                .collect(),
            archetype_ids: vec!["A1".into(), "A2".into()],
            loadings: vec![
                vec![1.0, 0.9, 0.1, 0.0], // A1 loads on P1, P2
                vec![0.0, 0.1, 0.9, 1.0], // A2 loads on P3, P4
            ],
            n_members: vec![2, 2],
            n_cohorts: vec![2, 2],
        };
        let cohort_labels = atlas.protein_labels.clone();
        let subject_ids = vec!["S01".into(), "S02".into(), "S03".into()];
        // Known coefficients per subject. Subject 1: (1.0, 0.0) —
        // loads purely A1. Subject 2: (0.0, 1.0) — purely A2.
        // Subject 3: (0.5, 0.5) — mix.
        let coefs = [[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]];
        let mut abundance = vec![vec![0.0_f64; atlas.protein_labels.len()]; subject_ids.len()];
        for (si, c) in coefs.iter().enumerate() {
            for pi in 0..atlas.protein_labels.len() {
                abundance[si][pi] = c[0] * atlas.loadings[0][pi] + c[1] * atlas.loadings[1][pi];
            }
        }
        (atlas, cohort_labels, subject_ids, abundance)
    }

    #[test]
    fn least_squares_recovers_planted_coefficients_to_tight_tolerance() {
        let (atlas, cohort_labels, subject_ids, abundance) = make_atlas_and_cohort();
        let out = project(
            &atlas,
            &cohort_labels,
            &subject_ids,
            &abundance,
            ProjectionMethod::LeastSquares,
        )
        .unwrap();
        let expected = [[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]];
        for (si, row) in out.activations.iter().enumerate() {
            for (ai, &v) in row.iter().enumerate() {
                assert!(
                    (v - expected[si][ai]).abs() < 1e-10,
                    "subject {si} archetype {ai}: got {v}, expected {}",
                    expected[si][ai]
                );
            }
        }
        for q in &out.qc {
            assert!(
                q.residual_norm < 1e-10,
                "residual should be ~0 when cohort lives in atlas span: {} -> {}",
                q.subject_id,
                q.residual_norm
            );
            assert!((q.coverage_fraction - 1.0).abs() < 1e-12);
        }
        assert!(out.atlas_proteins_missing_in_cohort.is_empty());
    }

    #[test]
    fn ridge_shrinks_l2_norm_of_coefficient_vector() {
        // Ridge minimises ||A β − x||² + λ||β||², so it shrinks the
        // L2 norm of β relative to LS — though individual
        // coefficients can move either way when archetypes are
        // correlated. The correct invariant is ||β_ridge|| ≤ ||β_ls||,
        // not per-coefficient bounding.
        let (atlas, cohort_labels, subject_ids, abundance) = make_atlas_and_cohort();
        let ls = project(
            &atlas,
            &cohort_labels,
            &subject_ids,
            &abundance,
            ProjectionMethod::LeastSquares,
        )
        .unwrap();
        let ridged = project(
            &atlas,
            &cohort_labels,
            &subject_ids,
            &abundance,
            ProjectionMethod::Ridge(1.0),
        )
        .unwrap();
        for si in 0..subject_ids.len() {
            let ls_norm: f64 = ls.activations[si].iter().map(|v| v * v).sum::<f64>().sqrt();
            let ridge_norm: f64 =
                ridged.activations[si].iter().map(|v| v * v).sum::<f64>().sqrt();
            assert!(
                ridge_norm <= ls_norm + 1e-12,
                "ridge L2 norm should not exceed LS L2 norm: {ridge_norm} vs {ls_norm}"
            );
        }
    }

    #[test]
    fn partial_coverage_emits_missing_labels() {
        let (atlas, _full_labels, subject_ids, abundance) = make_atlas_and_cohort();
        // Drop P4 from the cohort universe.
        let cohort_labels: Vec<String> = atlas.protein_labels[..3].to_vec();
        let abundance_partial: Vec<Vec<f64>> =
            abundance.iter().map(|row| row[..3].to_vec()).collect();
        let out = project(
            &atlas,
            &cohort_labels,
            &subject_ids,
            &abundance_partial,
            ProjectionMethod::Ridge(0.01),
        )
        .unwrap();
        assert_eq!(out.atlas_proteins_missing_in_cohort, vec!["P4".to_string()]);
        for q in &out.qc {
            assert!(
                q.coverage_fraction < 1.0,
                "expected partial coverage, got {}",
                q.coverage_fraction
            );
            assert_eq!(q.n_missing, 1);
        }
    }

    #[test]
    fn empty_atlas_is_rejected() {
        let atlas = Atlas {
            protein_labels: vec![],
            archetype_ids: vec![],
            loadings: vec![],
            n_members: vec![],
            n_cohorts: vec![],
        };
        let err = project(&atlas, &[], &["S01".into()], &[vec![]], ProjectionMethod::Ridge(0.01))
            .unwrap_err();
        assert!(err.contains("empty protein universe"), "got: {err}");
    }

    #[test]
    fn ridge_handles_ill_conditioned_atlas_without_failing() {
        // Two perfectly collinear atlas archetypes — LS would fail
        // with a non-SPD AᵀA, ridge must still solve.
        let atlas = Atlas {
            protein_labels: vec!["P1".into(), "P2".into()],
            archetype_ids: vec!["A1".into(), "A2".into()],
            loadings: vec![vec![1.0, 0.5], vec![1.0, 0.5]],
            n_members: vec![1, 1],
            n_cohorts: vec![1, 1],
        };
        let cohort_labels = atlas.protein_labels.clone();
        let subject_ids = vec!["S01".into()];
        let abundance = vec![vec![1.0, 0.5]];
        let ls_err = project(
            &atlas,
            &cohort_labels,
            &subject_ids,
            &abundance,
            ProjectionMethod::LeastSquares,
        )
        .unwrap_err();
        assert!(ls_err.contains("positive-definite"), "got: {ls_err}");
        let ridged = project(
            &atlas,
            &cohort_labels,
            &subject_ids,
            &abundance,
            ProjectionMethod::Ridge(0.1),
        )
        .unwrap();
        assert_eq!(ridged.activations.len(), 1);
    }
}
