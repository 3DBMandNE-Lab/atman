//! CSV diff helper for reproduction integration tests.
//!
//! Aligns two CSVs by row key (the first column), compares either strictly
//! (cell-for-cell string match) or numerically (with a caller-supplied
//! tolerance, skipping the first "Assay" column). Reports the first N
//! mismatches with enough context to diagnose.

#![allow(dead_code)]

use std::collections::HashMap;
use std::path::Path;

pub struct DiffReport {
    pub file_left: String,
    pub file_right: String,
    pub header_mismatch: Option<(Vec<String>, Vec<String>)>,
    pub missing_rows: Vec<String>,
    pub extra_rows: Vec<String>,
    pub cell_mismatches: Vec<CellMismatch>,
    pub max_numeric_delta: f64,
}

pub struct CellMismatch {
    pub row_key: String,
    pub column: String,
    pub left: String,
    pub right: String,
    pub delta: Option<f64>,
}

impl DiffReport {
    pub fn is_clean(&self) -> bool {
        self.header_mismatch.is_none()
            && self.missing_rows.is_empty()
            && self.extra_rows.is_empty()
            && self.cell_mismatches.is_empty()
    }

    pub fn assert_clean(&self) {
        if self.is_clean() {
            return;
        }
        let mut msg = format!(
            "\n=== DIFF FAIL: {} vs {} ===\n",
            self.file_left, self.file_right
        );
        if let Some((l, r)) = &self.header_mismatch {
            msg.push_str(&format!(
                "header mismatch:\n  left : {:?}\n  right: {:?}\n",
                l, r
            ));
        }
        if !self.missing_rows.is_empty() {
            msg.push_str(&format!(
                "rows in left missing from right ({} total, first 5): {:?}\n",
                self.missing_rows.len(),
                &self.missing_rows[..self.missing_rows.len().min(5)]
            ));
        }
        if !self.extra_rows.is_empty() {
            msg.push_str(&format!(
                "rows in right missing from left ({} total, first 5): {:?}\n",
                self.extra_rows.len(),
                &self.extra_rows[..self.extra_rows.len().min(5)]
            ));
        }
        if !self.cell_mismatches.is_empty() {
            msg.push_str(&format!(
                "cell mismatches (first 10 of {}) — max numeric delta: {:.2e}\n",
                self.cell_mismatches.len(),
                self.max_numeric_delta,
            ));
            for m in self.cell_mismatches.iter().take(10) {
                let d = match m.delta {
                    Some(d) => format!(" delta={:+.2e}", d),
                    None => String::new(),
                };
                msg.push_str(&format!(
                    "  row={} col={} left={:?} right={:?}{}\n",
                    m.row_key, m.column, m.left, m.right, d
                ));
            }
        }
        panic!("{}", msg);
    }
}

/// Load a CSV as (header, rows_by_key). Key = first column. Remaining columns
/// include the key itself (we keep the full row indexed by key for cell access).
fn load_csv(path: &Path) -> (Vec<String>, HashMap<String, Vec<String>>) {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)
        .unwrap_or_else(|e| panic!("opening {:?}: {}", path, e));
    let header: Vec<String> = reader
        .headers()
        .unwrap_or_else(|e| panic!("reading headers from {:?}: {}", path, e))
        .iter()
        .map(|s| s.to_string())
        .collect();
    let mut rows: HashMap<String, Vec<String>> = HashMap::new();
    for r in reader.records() {
        let r = r.unwrap_or_else(|e| panic!("reading record from {:?}: {}", path, e));
        let cells: Vec<String> = r.iter().map(|s| s.to_string()).collect();
        let key = cells[0].clone();
        rows.insert(key, cells);
    }
    (header, rows)
}

/// Strict cell-for-cell string diff. Use for filtered NPX files where
/// reproduction must be byte-exact.
pub fn diff_strict(left: &Path, right: &Path) -> DiffReport {
    let (lh, lr) = load_csv(left);
    let (rh, rr) = load_csv(right);
    let mut report = DiffReport {
        file_left: left.display().to_string(),
        file_right: right.display().to_string(),
        header_mismatch: None,
        missing_rows: vec![],
        extra_rows: vec![],
        cell_mismatches: vec![],
        max_numeric_delta: 0.0,
    };
    if lh != rh {
        report.header_mismatch = Some((lh.clone(), rh.clone()));
        return report;
    }
    for (key, lrow) in &lr {
        match rr.get(key) {
            Some(rrow) => {
                for (i, col) in lh.iter().enumerate() {
                    let lv = lrow.get(i).cloned().unwrap_or_default();
                    let rv = rrow.get(i).cloned().unwrap_or_default();
                    if lv != rv {
                        report.cell_mismatches.push(CellMismatch {
                            row_key: key.clone(),
                            column: col.clone(),
                            left: lv,
                            right: rv,
                            delta: None,
                        });
                    }
                }
            }
            None => report.missing_rows.push(key.clone()),
        }
    }
    for key in rr.keys() {
        if !lr.contains_key(key) {
            report.extra_rows.push(key.clone());
        }
    }
    report
}

/// Numeric cell diff with tolerance. First column must match exactly (the row
/// key, typically `Assay`). Remaining columns are parsed as f64 and compared
/// with `|delta| <= tol`. Empty cells must match empty cells.
pub fn diff_numeric(left: &Path, right: &Path, tol: f64) -> DiffReport {
    let (lh, lr) = load_csv(left);
    let (rh, rr) = load_csv(right);
    let mut report = DiffReport {
        file_left: left.display().to_string(),
        file_right: right.display().to_string(),
        header_mismatch: None,
        missing_rows: vec![],
        extra_rows: vec![],
        cell_mismatches: vec![],
        max_numeric_delta: 0.0,
    };
    if lh != rh {
        report.header_mismatch = Some((lh.clone(), rh.clone()));
        return report;
    }
    for (key, lrow) in &lr {
        match rr.get(key) {
            Some(rrow) => {
                for (i, col) in lh.iter().enumerate() {
                    let lv = lrow.get(i).cloned().unwrap_or_default();
                    let rv = rrow.get(i).cloned().unwrap_or_default();
                    if i == 0 {
                        if lv != rv {
                            report.cell_mismatches.push(CellMismatch {
                                row_key: key.clone(),
                                column: col.clone(),
                                left: lv,
                                right: rv,
                                delta: None,
                            });
                        }
                        continue;
                    }
                    let lf = lv.parse::<f64>().ok();
                    let rf = rv.parse::<f64>().ok();
                    match (lf, rf) {
                        (Some(a), Some(b)) => {
                            let d = (a - b).abs();
                            if d > report.max_numeric_delta {
                                report.max_numeric_delta = d;
                            }
                            if d > tol {
                                report.cell_mismatches.push(CellMismatch {
                                    row_key: key.clone(),
                                    column: col.clone(),
                                    left: lv,
                                    right: rv,
                                    delta: Some(a - b),
                                });
                            }
                        }
                        (None, None) if lv == rv => {}
                        _ => {
                            report.cell_mismatches.push(CellMismatch {
                                row_key: key.clone(),
                                column: col.clone(),
                                left: lv,
                                right: rv,
                                delta: None,
                            });
                        }
                    }
                }
            }
            None => report.missing_rows.push(key.clone()),
        }
    }
    for key in rr.keys() {
        if !lr.contains_key(key) {
            report.extra_rows.push(key.clone());
        }
    }
    report
}
