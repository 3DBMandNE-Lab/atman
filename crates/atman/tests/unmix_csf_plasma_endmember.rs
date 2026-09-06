//! Real-data sanity check for `atman decompose unmix` on the CSF
//! CrossDisease SIH cohort (requested by that manuscript session).
//!
//! # What it checks
//!
//! Albumin-ratio normalization (Reiber 1994) divides every protein by
//! albumin, which is itself the dominant plasma-leak marker. The claim
//! under test is that this collapses the bipolar plasma↔neuronal axis:
//!
//! - on the **primary** (log2 + median-normalized) matrix, at least one
//!   endmember is the plasma pole;
//! - on the **albnorm** matrix, no endmember is — its residual negative
//!   pole is immunoglobulins, not plasma proteins.
//!
//! # Why not cosine to the published Axis 1 vector
//!
//! The obvious criterion — cosine between an endmember and the negative
//! pole of the cross-cohort `archetype_A0002` consensus vector — does
//! not work, and was measured failing in both directions on 2026-09-06:
//! it reaches only −0.35 for the genuine plasma endmember on primary
//! (below any sensible threshold), and −0.66 for the *immunoglobulin*
//! endmember on albnorm (a false positive). The reason is that the
//! consensus vector's negative pole contains the immunoglobulins
//! alongside the plasma proteins (IGHA2 −0.230, IGHG3 −0.151, IGHA1
//! −0.096, IGHG2 −0.063, SERPINA1 −0.055), so an Ig-dominated endmember
//! scores as "plasma-like" on it. A whole-vector cosine cannot separate
//! the two poles that this test needs to tell apart.
//!
//! Instead the test scores three disjoint marker panels directly, on
//! per-protein z-scores over the cohort, and requires a plasma endmember
//! to be plasma-high, neuronal-low, *and* not merely immunoglobulin-high.
//!
//! # Running it
//!
//! The inputs are gitignored and live outside this repo, so the test is
//! `#[ignore]` and reads the canonical root from `ATMAN_CSF_CANONICAL_DIR`
//! (the directory holding `primary/sih/` and `albnorm/sih/`):
//!
//! ```text
//! ATMAN_CSF_CANONICAL_DIR=~/Cursor/CSF_CrossDisease/canonical \
//!   cargo test --test unmix_csf_plasma_endmember -- --ignored
//! ```

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Plasma-leak markers (CSF's list, minus HP which is not measured in
/// this cohort, and minus IGHA2 which is an immunoglobulin and belongs
/// to the panel below — keeping it here would blur the distinction the
/// test exists to draw).
const PLASMA: &[&str] = &["APCS", "APOB", "C4BPA", "ITIH3", "FGG", "FGB"];
/// Neuronal / brain-derived markers: the opposite pole of Axis 1.
const NEURONAL: &[&str] = &[
    "CBLN4", "SPON1", "NCAN", "SPOCK1", "CADM2", "NPTXR", "NPTX1",
];
/// Immunoglobulins and their close companions: the pole that survives
/// albumin normalization and that a naive criterion mistakes for plasma.
const IMMUNOGLOBULIN: &[&str] = &[
    "IGHA1", "IGHG1", "IGHG2", "IGHG3", "IGKC", "JCHAIN", "IGHM", "SERPINA1",
];

/// A plasma endmember must be plasma-high, neuronal-low, and not
/// explained by immunoglobulins. Threshold in per-protein SD units.
const PLASMA_Z_MIN: f64 = 0.5;

fn column_index(header: &str, name: &str, path: &Path) -> usize {
    header
        .trim_end_matches(['\n', '\r'])
        .split('\t')
        .position(|h| h == name)
        .unwrap_or_else(|| panic!("missing column {name:?} in {}", path.display()))
}

/// Per-protein (mean, sd) over observed, non-QC-dropped abundances.
/// Proteins with fewer than 3 observations or zero variance are omitted.
fn protein_moments(measurements: &Path) -> HashMap<String, (f64, f64)> {
    let text = std::fs::read_to_string(measurements)
        .unwrap_or_else(|e| panic!("read {}: {e}", measurements.display()));
    let mut lines = text.lines();
    let header = lines.next().expect("measurements.tsv is empty");
    let gene_col = column_index(header, "gene_symbol", measurements);
    let abundance_col = column_index(header, "abundance", measurements);
    let dropped_col = column_index(header, "dropped_by_qc", measurements);

    let mut values: HashMap<String, Vec<f64>> = HashMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() <= dropped_col.max(abundance_col).max(gene_col) || f[dropped_col].trim() == "1" {
            continue;
        }
        if let Ok(v) = f[abundance_col].trim().parse::<f64>() {
            if v.is_finite() {
                values.entry(f[gene_col].to_string()).or_default().push(v);
            }
        }
    }
    values
        .into_iter()
        .filter_map(|(gene, v)| {
            if v.len() < 3 {
                return None;
            }
            let n = v.len() as f64;
            let mean = v.iter().sum::<f64>() / n;
            let sd = (v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
            (sd > 0.0).then_some((gene, (mean, sd)))
        })
        .collect()
}

/// endmember_id → (protein → loading).
fn read_endmembers(path: &Path) -> Vec<(String, HashMap<String, f64>)> {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let mut lines = text.lines();
    let header = lines.next().expect("endmembers.tsv is empty");
    let id_col = column_index(header, "endmember_id", path);
    let protein_col = column_index(header, "protein", path);
    let loading_col = column_index(header, "loading", path);

    let mut by_id: Vec<(String, HashMap<String, f64>)> = Vec::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let Ok(loading) = f[loading_col].trim().parse::<f64>() else {
            continue;
        };
        let id = f[id_col].to_string();
        match by_id.iter_mut().find(|(existing, _)| *existing == id) {
            Some((_, m)) => {
                m.insert(f[protein_col].to_string(), loading);
            }
            None => {
                let mut m = HashMap::new();
                m.insert(f[protein_col].to_string(), loading);
                by_id.push((id, m));
            }
        }
    }
    by_id.sort_by(|a, b| a.0.cmp(&b.0));
    by_id
}

#[derive(Debug, Clone, Copy)]
struct PanelScores {
    plasma: f64,
    neuronal: f64,
    immunoglobulin: f64,
}

impl PanelScores {
    /// The endmember reads as a plasma-leak pole: plasma markers
    /// elevated, neuronal markers depressed, and immunoglobulins not
    /// carrying the signal on their own.
    fn is_plasma_pole(&self) -> bool {
        self.plasma > PLASMA_Z_MIN && self.neuronal < 0.0 && self.immunoglobulin < self.plasma
    }
}

fn score_panels(
    loadings: &HashMap<String, f64>,
    moments: &HashMap<String, (f64, f64)>,
) -> PanelScores {
    let panel_mean = |panel: &[&str]| -> f64 {
        let z: Vec<f64> = panel
            .iter()
            .filter_map(|g| {
                let loading = loadings.get(*g)?;
                let (mean, sd) = moments.get(*g)?;
                Some((loading - mean) / sd)
            })
            .collect();
        assert!(
            z.len() >= 4,
            "only {} of {} panel markers are measured; the panel is too thin to score",
            z.len(),
            panel.len()
        );
        z.iter().sum::<f64>() / z.len() as f64
    };
    PanelScores {
        plasma: panel_mean(PLASMA),
        neuronal: panel_mean(NEURONAL),
        immunoglobulin: panel_mean(IMMUNOGLOBULIN),
    }
}

/// Run `decompose unmix` on one canonical dir and score every endmember.
fn unmix_and_score(input_dir: &Path, output_dir: &Path) -> Vec<(String, PanelScores)> {
    assert!(
        input_dir.join("measurements.tsv").is_file(),
        "no measurements.tsv under {}",
        input_dir.display()
    );
    let out = Command::new(env!("CARGO_BIN_EXE_atman"))
        .args([
            "decompose",
            "unmix",
            "--input-dir",
            input_dir.to_str().unwrap(),
            // Fixed k rather than `--k auto`: the elbow rule selects 8
            // here and 6 on albnorm, and the finding is stable at k =
            // 4, 5 and 8 on primary, so a small explicit k keeps the
            // test's claim independent of the selection heuristic.
            "--k",
            "4",
            // ~26% of cells are missing or QC-dropped in this cohort.
            "--max-missing-fraction",
            "0.3",
            "--impute",
            "mean",
            "--output-dir",
            output_dir.to_str().unwrap(),
        ])
        .output()
        .expect("run atman decompose unmix");
    assert!(
        out.status.success(),
        "decompose unmix failed on {}:\n{}",
        input_dir.display(),
        String::from_utf8_lossy(&out.stderr)
    );

    let moments = protein_moments(&input_dir.join("measurements.tsv"));
    read_endmembers(&output_dir.join("endmembers.tsv"))
        .into_iter()
        .map(|(id, loadings)| (id, score_panels(&loadings, &moments)))
        .collect()
}

fn canonical_root() -> PathBuf {
    let raw = std::env::var("ATMAN_CSF_CANONICAL_DIR").unwrap_or_else(|_| {
        panic!(
            "set ATMAN_CSF_CANONICAL_DIR to the CSF CrossDisease canonical directory \
             (the one containing primary/sih/ and albnorm/sih/); the inputs are gitignored \
             and live outside this repo"
        )
    });
    let root = PathBuf::from(shellexpand_tilde(&raw));
    assert!(
        root.join("primary/sih").is_dir() && root.join("albnorm/sih").is_dir(),
        "ATMAN_CSF_CANONICAL_DIR={} does not contain primary/sih and albnorm/sih",
        root.display()
    );
    root
}

fn shellexpand_tilde(p: &str) -> String {
    match p.strip_prefix("~/") {
        Some(rest) => match std::env::var("HOME") {
            Ok(home) => format!("{home}/{rest}"),
            Err(_) => p.to_string(),
        },
        None => p.to_string(),
    }
}

#[test]
#[ignore = "needs the gitignored CSF CrossDisease canonical inputs; see module docs"]
fn unmix_finds_plasma_endmember_on_primary_but_not_on_albumin_normalized() {
    let root = canonical_root();
    let tmp = tempfile::tempdir().unwrap();

    let primary = unmix_and_score(&root.join("primary/sih"), &tmp.path().join("primary"));
    let albnorm = unmix_and_score(&root.join("albnorm/sih"), &tmp.path().join("albnorm"));

    let describe = |label: &str, scored: &[(String, PanelScores)]| -> String {
        let mut s = format!("{label}:\n");
        for (id, p) in scored {
            s.push_str(&format!(
                "  {id} plasma={:+.2} neuronal={:+.2} ig={:+.2} plasma_pole={}\n",
                p.plasma,
                p.neuronal,
                p.immunoglobulin,
                p.is_plasma_pole()
            ));
        }
        s
    };
    let report = format!(
        "{}{}",
        describe("primary", &primary),
        describe("albnorm", &albnorm)
    );

    let primary_poles: Vec<&String> = primary
        .iter()
        .filter(|(_, p)| p.is_plasma_pole())
        .map(|(id, _)| id)
        .collect();
    assert!(
        !primary_poles.is_empty(),
        "expected at least one plasma endmember on the primary matrix\n{report}"
    );

    let albnorm_poles: Vec<&String> = albnorm
        .iter()
        .filter(|(_, p)| p.is_plasma_pole())
        .map(|(id, _)| id)
        .collect();
    assert!(
        albnorm_poles.is_empty(),
        "albumin normalization should leave no plasma endmember, found {albnorm_poles:?}\n{report}"
    );

    // The albnorm run must still produce a usable decomposition — an
    // empty or degenerate one would satisfy the assertion above for the
    // wrong reason.
    assert!(
        albnorm.len() >= 4,
        "albnorm decomposition returned {} endmembers; expected 4\n{report}",
        albnorm.len()
    );
}
