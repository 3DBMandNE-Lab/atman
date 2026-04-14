//! Log2 fold change per spec §7.4.

use std::collections::BTreeMap;

#[derive(Debug, Clone)]
pub struct Comparison {
    pub a: String,
    pub b: String,
}

pub struct FoldChangeInput {
    /// (panel, assay) -> (exposure -> Vec<(participant, value)>)
    panels: BTreeMap<(String, String), BTreeMap<String, Vec<(String, f64)>>>,
}

impl FoldChangeInput {
    /// Construct from a sequence of cells. `S1..S4: AsRef<str>` so both
    /// `&str` literals (tests) and owned `String` (CLI) compile.
    /// Missing values (`None`) are dropped silently.
    ///
    /// Note: if an (assay, panel) pair has zero non-missing cells it will
    /// NOT appear in the output. Use `ensure_assay` to force an assay into
    /// the output universe even if no cells are contributed — this matches
    /// Dube's behavior of emitting a row per assay regardless of data.
    pub fn from_cells<I, S1, S2, S3, S4>(cells: I) -> Self
    where
        I: IntoIterator<Item = (S1, S2, S3, S4, Option<f64>)>,
        S1: AsRef<str>,
        S2: AsRef<str>,
        S3: AsRef<str>,
        S4: AsRef<str>,
    {
        let mut panels: BTreeMap<(String, String), BTreeMap<String, Vec<(String, f64)>>> =
            BTreeMap::new();
        for (panel, assay, participant, exposure, value) in cells {
            if let Some(v) = value {
                panels
                    .entry((panel.as_ref().to_string(), assay.as_ref().to_string()))
                    .or_default()
                    .entry(exposure.as_ref().to_string())
                    .or_default()
                    .push((participant.as_ref().to_string(), v));
            }
        }
        Self { panels }
    }

    /// Ensure an `(panel, assay)` pair appears in the output universe even
    /// if no cells are contributed. The resulting fold-change will be `None`
    /// for every comparison on such pairs. Matches Dube's behavior of
    /// emitting a row per assay in the panel regardless of whether the
    /// underlying data can produce a fold change.
    pub fn ensure_assay(&mut self, panel: &str, assay: &str) {
        self.panels
            .entry((panel.to_string(), assay.to_string()))
            .or_default();
    }
}

pub struct FoldChangePanel {
    pub panel: String,
    pub assays: Vec<String>,
    /// Outer index: assay. Inner index: comparisons (aligned with the
    /// `comparisons` argument passed to `compute_log2_fc`).
    pub values: Vec<Vec<Option<f64>>>,
}

pub struct FoldChangeOutput {
    pub panels: Vec<FoldChangePanel>,
}

pub fn compute_log2_fc(
    input: &FoldChangeInput,
    comparisons: &[Comparison],
) -> FoldChangeOutput {
    let mut per_panel: BTreeMap<String, Vec<(String, &BTreeMap<String, Vec<(String, f64)>>)>> =
        BTreeMap::new();
    for ((panel, assay), by_exposure) in &input.panels {
        per_panel
            .entry(panel.clone())
            .or_default()
            .push((assay.clone(), by_exposure));
    }

    let mut out_panels = Vec::new();
    for (panel_name, mut assay_rows) in per_panel {
        assay_rows.sort_by(|a, b| a.0.cmp(&b.0));
        let assays: Vec<String> = assay_rows.iter().map(|(a, _)| a.clone()).collect();
        let mut values: Vec<Vec<Option<f64>>> = Vec::with_capacity(assays.len());
        for (_assay, by_exposure) in &assay_rows {
            let mut row: Vec<Option<f64>> = Vec::with_capacity(comparisons.len());
            for c in comparisons {
                row.push(compute_one_fc(by_exposure, &c.a, &c.b));
            }
            values.push(row);
        }
        out_panels.push(FoldChangePanel {
            panel: panel_name,
            assays,
            values,
        });
    }
    FoldChangeOutput { panels: out_panels }
}

fn compute_one_fc(
    by_exposure: &BTreeMap<String, Vec<(String, f64)>>,
    a: &str,
    b: &str,
) -> Option<f64> {
    let va = by_exposure.get(a)?;
    let vb = by_exposure.get(b)?;
    if va.is_empty() || vb.is_empty() {
        return None;
    }
    let mean_a = mean_linear(va);
    let mean_b = mean_linear(vb);
    if mean_b == 0.0 {
        return None;
    }
    Some((mean_a / mean_b).log2())
}

fn mean_linear(values: &[(String, f64)]) -> f64 {
    let sum: f64 = values.iter().map(|(_, v)| 2.0_f64.powf(*v)).sum();
    sum / (values.len() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mean_linear_basic() {
        let v = vec![("a".into(), 0.0), ("b".into(), 1.0), ("c".into(), 2.0)];
        assert!((mean_linear(&v) - 7.0 / 3.0).abs() < 1e-12);
    }
}
