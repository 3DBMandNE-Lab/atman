//! Contrast manifest shared by `axes contrast` and `axes displacement`.
//!
//! Columns: `label`, `cohort`, `case`, `control`, `family`; optional
//! `condition_col` (default `condition`) and `subset` (predicates such as
//! `diagnosis_group!=MS`). `case`/`control` accept `a|b` lists; `control=*`
//! means every other level inside the subset.

use anyhow::{anyhow, bail, Context, Result};
use std::path::Path;

use super::Context as RowContext;
use crate::design::{parse_predicates, predicates_match, Predicate};

#[derive(Debug, Clone)]
pub struct Contrast {
    pub label: String,
    pub cohort: String,
    pub case: Vec<String>,
    pub control: Vec<String>,
    pub control_is_rest: bool,
    pub family: String,
    pub condition_col: String,
    pub subset: Vec<Predicate>,
}

fn split_levels(text: &str) -> Vec<String> {
    text.split('|')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect()
}

pub fn read_contrast_manifest(path: &Path) -> Result<Vec<Contrast>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = reader.headers()?.clone();
    let col = |name: &str| headers.iter().position(|h| h.trim() == name);
    let need =
        |name: &str| col(name).ok_or_else(|| anyhow!("missing column {:?} in {:?}", name, path));
    let (c_label, c_cohort, c_case, c_control, c_family) = (
        need("label")?,
        need("cohort")?,
        need("case")?,
        need("control")?,
        need("family")?,
    );
    let c_cond = col("condition_col");
    let c_subset = col("subset");
    let mut out = Vec::new();
    for row in reader.records() {
        let row = row?;
        let get = |i: usize| row.get(i).unwrap_or("").trim().to_string();
        let label = get(c_label);
        if label.is_empty() {
            continue;
        }
        if out.iter().any(|c: &Contrast| c.label == label) {
            bail!("duplicate contrast label {:?} in {:?}", label, path);
        }
        let control_text = get(c_control);
        let control_is_rest = control_text == "*";
        let case = split_levels(&get(c_case));
        if case.is_empty() {
            bail!("contrast {:?}: empty case", label);
        }
        let control = if control_is_rest {
            Vec::new()
        } else {
            split_levels(&control_text)
        };
        if !control_is_rest && control.is_empty() {
            bail!(
                "contrast {:?}: empty control (use * for all other levels)",
                label
            );
        }
        let condition_col = c_cond
            .map(get)
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| "condition".to_string());
        let subset_text = c_subset.map(get).unwrap_or_default();
        let subset = if subset_text.is_empty() {
            Vec::new()
        } else {
            parse_predicates(&subset_text)?
        };
        out.push(Contrast {
            label,
            cohort: get(c_cohort),
            case,
            control,
            control_is_rest,
            family: get(c_family),
            condition_col,
            subset,
        });
    }
    if out.is_empty() {
        bail!("no contrasts in {:?}", path);
    }
    Ok(out)
}

#[derive(Debug, Clone, Default)]
pub struct Selection {
    pub case: Vec<usize>,
    pub control: Vec<usize>,
}

/// Rows of the contrast's cohort passing `subset`, split into case and
/// control by the value of `condition_col`.
pub fn select_rows(c: &Contrast, ctx: &RowContext) -> Result<Selection> {
    let mut sel = Selection::default();
    for i in 0..ctx.table.len() {
        if ctx.table.cohort[i] != c.cohort {
            continue;
        }
        if !c.subset.is_empty() && !predicates_match(&c.subset, &|col| ctx.text(i, col)) {
            continue;
        }
        let Some(level) = ctx.text(i, &c.condition_col) else {
            continue;
        };
        if c.case.contains(&level) {
            sel.case.push(i);
        } else if c.control_is_rest || c.control.contains(&level) {
            sel.control.push(i);
        }
    }
    if sel.case.is_empty() {
        bail!(
            "contrast {:?}: no case rows (cohort {:?}, {} = {})",
            c.label,
            c.cohort,
            c.condition_col,
            c.case.join("|")
        );
    }
    if sel.control.is_empty() {
        bail!(
            "contrast {:?}: no control rows (cohort {:?})",
            c.label,
            c.cohort
        );
    }
    Ok(sel)
}
