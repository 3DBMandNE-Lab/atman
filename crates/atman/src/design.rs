//! Sample-metadata frames and formula → design-matrix construction shared
//! by the `axes` command group and `residuals`.
//!
//! A [`CovariateFrame`] is the union of every `sample_id`-keyed TSV a
//! command was pointed at (cohort `samples.tsv` files, `--covariates-tsv`
//! files). A [`DesignSpec`] is a parsed `~ case + z(age) + sex` formula and
//! [`build_design`] turns it into a complete-case numeric design matrix
//! following the covariate-typing rule used by `de --test ols`: a column
//! whose trimmed non-empty values all parse as numbers is numeric;
//! otherwise it is categorical and one-hot encoded with the alphabetically
//! first level as reference.

use anyhow::{anyhow, bail, Context, Result};
use atman_core::expr::{evaluate_column, split_top_level, Expr};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

#[derive(Debug, Default, Clone)]
pub struct CovariateFrame {
    columns: BTreeSet<String>,
    by_sample: HashMap<String, BTreeMap<String, String>>,
}

impl CovariateFrame {
    pub fn new() -> Self {
        Self::default()
    }

    /// Load every column of a `sample_id`-keyed TSV. Values are trimmed and
    /// empty cells skipped. Returns the number of data rows read. A
    /// (sample, column) pair that already holds a different value is an
    /// error so two files cannot silently disagree.
    pub fn load_tsv(&mut self, path: &Path) -> Result<usize> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .flexible(true)
            .from_path(path)
            .with_context(|| format!("opening {:?}", path))?;
        let headers = reader.headers()?.clone();
        let sample_col = headers
            .iter()
            .position(|h| h.trim() == "sample_id")
            .ok_or_else(|| anyhow!("missing column \"sample_id\" in {:?}", path))?;
        let names: Vec<String> = headers.iter().map(|h| h.trim().to_string()).collect();
        let mut n = 0usize;
        for row in reader.records() {
            let row = row.with_context(|| format!("reading {:?}", path))?;
            let sid = row.get(sample_col).unwrap_or("").trim().to_string();
            if sid.is_empty() {
                continue;
            }
            n += 1;
            let entry = self.by_sample.entry(sid.clone()).or_default();
            for (i, name) in names.iter().enumerate() {
                if i == sample_col || name.is_empty() {
                    continue;
                }
                let value = row.get(i).unwrap_or("").trim();
                if value.is_empty() {
                    continue;
                }
                if let Some(existing) = entry.get(name) {
                    if existing != value {
                        bail!(
                            "sample {:?} column {:?} is {:?} in an earlier file but {:?} in {:?}",
                            sid,
                            name,
                            existing,
                            value,
                            path
                        );
                    }
                } else {
                    entry.insert(name.clone(), value.to_string());
                }
                self.columns.insert(name.clone());
            }
        }
        Ok(n)
    }

    pub fn columns(&self) -> &BTreeSet<String> {
        &self.columns
    }

    pub fn has_column(&self, col: &str) -> bool {
        self.columns.contains(col)
    }

    pub fn get(&self, sample_id: &str, col: &str) -> Option<&str> {
        self.by_sample
            .get(sample_id)
            .and_then(|m| m.get(col))
            .map(String::as_str)
    }

    pub fn contains_sample(&self, sample_id: &str) -> bool {
        self.by_sample.contains_key(sample_id)
    }

    pub fn sample_ids(&self) -> impl Iterator<Item = &String> {
        self.by_sample.keys()
    }
}

/// Parse `label=path,label=path`.
pub fn parse_labeled_paths(spec: &str) -> Result<Vec<(String, PathBuf)>> {
    let mut out = Vec::new();
    for chunk in spec.split(',') {
        let chunk = chunk.trim();
        if chunk.is_empty() {
            continue;
        }
        let (label, path) = chunk
            .split_once('=')
            .ok_or_else(|| anyhow!("expected label=path, got {:?}", chunk))?;
        let label = label.trim();
        let path = path.trim();
        if label.is_empty() || path.is_empty() {
            bail!("expected label=path, got {:?}", chunk);
        }
        out.push((label.to_string(), PathBuf::from(path)));
    }
    if out.is_empty() {
        bail!("no label=path entries in {:?}", spec);
    }
    Ok(out)
}

/// A categorical term whose fitted rows carry fewer than two levels.
/// Raised as a typed error so bootstrap loops can skip the replicate.
#[derive(Debug, Clone)]
pub struct DegenerateFactor {
    pub column: String,
}

impl std::fmt::Display for DegenerateFactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "categorical covariate {:?} has fewer than two levels among the fitted rows",
            self.column
        )
    }
}

impl std::error::Error for DegenerateFactor {}

/// How protein groups (assays) that share a gene symbol are reduced to one
/// value per sample. Shared by `de`, `residuals`, and `score weighted`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum CollapseGenes {
    /// Keep one assay per gene: the lexically first assay id.
    None,
    /// Per sample, the mean of the non-missing assays sharing the symbol.
    Mean,
    /// Keep the single assay observed in the most samples (ties: lexically first).
    MaxObserved,
}

impl CollapseGenes {
    pub fn as_str(&self) -> &'static str {
        match self {
            CollapseGenes::None => "none",
            CollapseGenes::Mean => "mean",
            CollapseGenes::MaxObserved => "max-observed",
        }
    }

    /// Choose the representative assay id for a gene from
    /// `(assay_id, n_observed)` pairs.
    pub fn representative<'a>(&self, counts: &'a BTreeMap<String, usize>) -> Option<&'a String> {
        match self {
            CollapseGenes::MaxObserved => counts
                .iter()
                .max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0)))
                .map(|(a, _)| a),
            _ => counts.keys().next(),
        }
    }

    /// Reduce one sample's `(assay_id, value)` list to a single value.
    pub fn collapse(&self, values: &[(String, f64)], representative: &str) -> Option<f64> {
        match self {
            CollapseGenes::Mean => {
                if values.is_empty() {
                    None
                } else {
                    Some(values.iter().map(|(_, v)| *v).sum::<f64>() / values.len() as f64)
                }
            }
            _ => values
                .iter()
                .find(|(a, _)| a == representative)
                .map(|(_, v)| *v),
        }
    }
}

#[derive(Debug, Clone)]
pub enum Term {
    /// The contrast indicator (e.g. `case`), supplied by the caller per row.
    Indicator(String),
    /// A raw metadata column, typed numeric or categorical at build time.
    Column(String),
    /// A numeric expression over metadata columns.
    Expression { text: String, expr: Expr },
}

#[derive(Debug, Clone)]
pub struct DesignSpec {
    pub terms: Vec<Term>,
    pub formula: String,
}

/// Parse `~ term + term + ...`. When `indicator` is `Some(name)` the formula
/// must contain that bare term exactly once.
pub fn parse_design(formula: &str, indicator: Option<&str>) -> Result<DesignSpec> {
    let trimmed = formula.trim();
    let rhs = trimmed
        .strip_prefix('~')
        .ok_or_else(|| anyhow!("design {:?} must start with `~`", formula))?
        .trim();
    if rhs.is_empty() {
        bail!("design {:?} has no terms", formula);
    }
    let mut terms = Vec::new();
    let mut seen = BTreeSet::new();
    for raw in split_top_level(rhs, '+') {
        if raw == "1" {
            continue;
        }
        if raw == "0" || raw == "-1" {
            bail!("intercept removal is not supported; Atman includes an intercept");
        }
        if !seen.insert(raw.clone()) {
            bail!("duplicate term {:?} in design {:?}", raw, formula);
        }
        if Some(raw.as_str()) == indicator {
            terms.push(Term::Indicator(raw));
            continue;
        }
        let expr = Expr::parse(&raw)
            .map_err(|e| anyhow!("design {:?}: term {:?}: {}", formula, raw, e))?;
        if let Some(name) = expr.bare_identifier() {
            terms.push(Term::Column(name.to_string()));
        } else {
            terms.push(Term::Expression { text: raw, expr });
        }
    }
    if let Some(name) = indicator {
        if !terms.iter().any(|t| matches!(t, Term::Indicator(_))) {
            bail!(
                "design {:?} must include the contrast term `{}`",
                formula,
                name
            );
        }
    }
    Ok(DesignSpec {
        terms,
        formula: trimmed.to_string(),
    })
}

#[derive(Debug, Clone)]
pub struct DesignMatrix {
    /// One row per kept sample: `[1, columns...]`.
    pub rows: Vec<Vec<f64>>,
    /// Column labels aligned with `rows[i]`.
    pub labels: Vec<String>,
    /// Indices (into the caller's `0..n`) of the kept, complete-case rows.
    pub kept: Vec<usize>,
    /// Column index of the indicator term, when the spec has one.
    pub indicator_col: Option<usize>,
    /// For every column, the index of the `spec.terms` entry that produced
    /// it (`None` for the intercept). A multi-level categorical term owns
    /// several consecutive columns.
    pub term_index: Vec<Option<usize>>,
}

enum ColumnKind {
    Numeric,
    Categorical(Vec<String>),
}

/// Build the complete-case design matrix. `text(i, col)` returns the raw
/// metadata value of row `i`; `indicator(i)` returns the indicator value for
/// row `i` (`None` excludes the row).
pub fn build_design(
    spec: &DesignSpec,
    n: usize,
    text: &dyn Fn(usize, &str) -> Option<String>,
    indicator: &dyn Fn(usize) -> Option<f64>,
) -> Result<DesignMatrix> {
    let num = |i: usize, name: &str| -> Option<f64> {
        text(i, name)
            .and_then(|s| s.trim().parse::<f64>().ok())
            .filter(|v| v.is_finite())
    };
    let mut kept = Vec::new();
    for i in 0..n {
        let complete = spec.terms.iter().all(|term| match term {
            Term::Indicator(_) => indicator(i).is_some(),
            Term::Column(c) => text(i, c).map(|s| !s.trim().is_empty()).unwrap_or(false),
            Term::Expression { expr, .. } => expr.inner().eval(&|name| num(i, name)).is_some(),
        });
        if complete {
            kept.push(i);
        }
    }
    if kept.is_empty() {
        bail!(
            "no rows with complete covariates for design {:?}",
            spec.formula
        );
    }

    let mut kinds: Vec<Option<ColumnKind>> = Vec::with_capacity(spec.terms.len());
    for term in &spec.terms {
        match term {
            Term::Column(c) => {
                let values: Vec<String> = kept
                    .iter()
                    .map(|&i| text(i, c).unwrap_or_default().trim().to_string())
                    .collect();
                if values
                    .iter()
                    .all(|v| v.parse::<f64>().map(|x| x.is_finite()).unwrap_or(false))
                {
                    kinds.push(Some(ColumnKind::Numeric));
                } else {
                    let levels: Vec<String> = values
                        .iter()
                        .cloned()
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .collect();
                    if levels.len() < 2 {
                        return Err(anyhow::Error::new(DegenerateFactor { column: c.clone() }));
                    }
                    kinds.push(Some(ColumnKind::Categorical(levels)));
                }
            }
            _ => kinds.push(None),
        }
    }

    let mut expr_values: BTreeMap<String, Vec<Option<f64>>> = BTreeMap::new();
    for term in &spec.terms {
        if let Term::Expression { text: label, expr } = term {
            let col = evaluate_column(expr, kept.len(), &|k, name| num(kept[k], name));
            if col.iter().any(|v| v.is_none()) {
                bail!(
                    "expression {:?} is constant or undefined over the fitted rows",
                    label
                );
            }
            expr_values.insert(label.clone(), col);
        }
    }

    let mut labels = vec!["(Intercept)".to_string()];
    let mut term_index: Vec<Option<usize>> = vec![None];
    let mut indicator_col = None;
    for (t, term) in spec.terms.iter().enumerate() {
        match term {
            Term::Indicator(name) => {
                indicator_col = Some(labels.len());
                labels.push(name.clone());
                term_index.push(Some(t));
            }
            Term::Column(c) => match kinds[t].as_ref().expect("column kind") {
                ColumnKind::Numeric => {
                    labels.push(c.clone());
                    term_index.push(Some(t));
                }
                ColumnKind::Categorical(levels) => {
                    for level in levels.iter().skip(1) {
                        labels.push(format!("{c}{level}"));
                        term_index.push(Some(t));
                    }
                }
            },
            Term::Expression { text: label, .. } => {
                labels.push(label.clone());
                term_index.push(Some(t));
            }
        }
    }

    let mut rows = Vec::with_capacity(kept.len());
    for (k, &i) in kept.iter().enumerate() {
        let mut row = Vec::with_capacity(labels.len());
        row.push(1.0);
        for (t, term) in spec.terms.iter().enumerate() {
            match term {
                Term::Indicator(_) => {
                    row.push(indicator(i).expect("indicator present for kept row"))
                }
                Term::Column(c) => {
                    let value = text(i, c).unwrap_or_default().trim().to_string();
                    match kinds[t].as_ref().expect("column kind") {
                        ColumnKind::Numeric => {
                            row.push(value.parse::<f64>().expect("numeric column"))
                        }
                        ColumnKind::Categorical(levels) => {
                            for level in levels.iter().skip(1) {
                                row.push(if &value == level { 1.0 } else { 0.0 });
                            }
                        }
                    }
                }
                Term::Expression { text: label, .. } => {
                    row.push(expr_values[label][k].expect("expression evaluated"));
                }
            }
        }
        rows.push(row);
    }
    Ok(DesignMatrix {
        rows,
        labels,
        kept,
        indicator_col,
        term_index,
    })
}

#[derive(Debug, Clone, PartialEq)]
pub enum PredicateOp {
    Eq,
    Ne,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Predicate {
    pub column: String,
    pub op: PredicateOp,
    pub value: String,
}

/// Parse `col!=value; col==value; col=value` (semicolon-separated).
pub fn parse_predicates(text: &str) -> Result<Vec<Predicate>> {
    let mut out = Vec::new();
    for raw in text.split(';') {
        let raw = raw.trim();
        if raw.is_empty() {
            continue;
        }
        let (column, op, value) = if let Some((c, v)) = raw.split_once("!=") {
            (c, PredicateOp::Ne, v)
        } else if let Some((c, v)) = raw.split_once("==") {
            (c, PredicateOp::Eq, v)
        } else if let Some((c, v)) = raw.split_once('=') {
            (c, PredicateOp::Eq, v)
        } else {
            bail!(
                "predicate {:?} must be col==value, col=value, or col!=value",
                raw
            );
        };
        let column = column.trim();
        let value = value.trim();
        if column.is_empty() || value.is_empty() {
            bail!("predicate {:?} has an empty column or value", raw);
        }
        out.push(Predicate {
            column: column.to_string(),
            op,
            value: value.to_string(),
        });
    }
    if out.is_empty() {
        bail!("no predicates in {:?}", text);
    }
    Ok(out)
}

/// All predicates must hold. A missing value fails `==` and passes `!=`.
pub fn predicates_match(preds: &[Predicate], get: &dyn Fn(&str) -> Option<String>) -> bool {
    preds.iter().all(|p| {
        let actual = get(&p.column)
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty());
        match (&p.op, actual) {
            (PredicateOp::Eq, Some(v)) => v == p.value,
            (PredicateOp::Eq, None) => false,
            (PredicateOp::Ne, Some(v)) => v != p.value,
            (PredicateOp::Ne, None) => true,
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame_from(text: &str) -> CovariateFrame {
        let tmp = tempfile::tempdir().unwrap();
        let p = tmp.path().join("samples.tsv");
        std::fs::write(&p, text).unwrap();
        let mut f = CovariateFrame::new();
        f.load_tsv(&p).unwrap();
        f
    }

    #[test]
    fn frame_trims_and_skips_empty_cells() {
        let f = frame_from("sample_id\tage\tsex\nS1\t50\tM \nS2\t\tF\n");
        assert_eq!(f.get("S1", "sex"), Some("M"));
        assert_eq!(f.get("S1", "age"), Some("50"));
        assert_eq!(f.get("S2", "age"), None);
        assert!(f.has_column("sex"));
        assert!(!f.has_column("qalb"));
        assert!(f.contains_sample("S2"));
    }

    #[test]
    fn frame_rejects_conflicting_values_across_files() {
        let tmp = tempfile::tempdir().unwrap();
        let a = tmp.path().join("a.tsv");
        let b = tmp.path().join("b.tsv");
        std::fs::write(&a, "sample_id\tage\nS1\t50\n").unwrap();
        std::fs::write(&b, "sample_id\tage\nS1\t51\n").unwrap();
        let mut f = CovariateFrame::new();
        f.load_tsv(&a).unwrap();
        assert!(f.load_tsv(&b).is_err());
    }

    #[test]
    fn labeled_paths_parse() {
        let v = parse_labeled_paths("sih=/a/sih, ms=/b/ms").unwrap();
        assert_eq!(v[0].0, "sih");
        assert_eq!(v[1].1, PathBuf::from("/b/ms"));
        assert!(parse_labeled_paths("nolabel").is_err());
    }

    #[test]
    fn design_parse_classifies_terms() {
        let spec = parse_design("~ case + z(age) + sex + log10(QAlb + 1)", Some("case")).unwrap();
        assert_eq!(spec.terms.len(), 4);
        assert!(matches!(&spec.terms[0], Term::Indicator(n) if n == "case"));
        assert!(matches!(&spec.terms[1], Term::Expression { text, .. } if text == "z(age)"));
        assert!(matches!(&spec.terms[2], Term::Column(c) if c == "sex"));
        assert!(parse_design("~ age", Some("case")).is_err());
        assert!(parse_design("age", None).is_err());
        assert!(parse_design("~ age + age", None).is_err());
        assert!(parse_design("~ 0 + age", None).is_err());
        assert_eq!(parse_design("~ 1 + age", None).unwrap().terms.len(), 1);
    }

    #[test]
    fn design_build_encodes_numeric_categorical_expression_and_indicator() {
        let f = frame_from(
            "sample_id\tage\tsex\tQAlb\nS1\t10\tF\t10\nS2\t20\tM\t100\nS3\t30\tF\t1000\nS4\t\tM\t10\nS5\t40\tM\t100\n",
        );
        let ids = ["S1", "S2", "S3", "S4", "S5"];
        let spec = parse_design("~ case + z(age) + sex + log10(QAlb)", Some("case")).unwrap();
        let dm = build_design(
            &spec,
            5,
            &|i, col| f.get(ids[i], col).map(str::to_string),
            &|i| Some(if i % 2 == 0 { 1.0 } else { 0.0 }),
        )
        .unwrap();
        assert_eq!(dm.kept, vec![0, 1, 2, 4]);
        assert_eq!(
            dm.labels,
            vec!["(Intercept)", "case", "z(age)", "sexM", "log10(QAlb)"]
        );
        assert_eq!(dm.indicator_col, Some(1));
        assert_eq!(
            dm.term_index,
            vec![None, Some(0), Some(1), Some(2), Some(3)]
        );
        assert!((dm.rows[0][2] + 1.161895003862225).abs() < 1e-9);
        assert!((dm.rows[3][2] - 1.161895003862225).abs() < 1e-9);
        assert_eq!(dm.rows[1][3], 1.0);
        assert_eq!(dm.rows[0][3], 0.0);
        assert!((dm.rows[2][4] - 3.0).abs() < 1e-12);
        assert_eq!(dm.rows[0][1], 1.0);
        assert_eq!(dm.rows[1][1], 0.0);
    }

    #[test]
    fn design_build_rejects_single_level_categorical() {
        let f = frame_from("sample_id\tsex\nS1\tF\nS2\tF\n");
        let ids = ["S1", "S2"];
        let spec = parse_design("~ sex", None).unwrap();
        let err = build_design(
            &spec,
            2,
            &|i, col| f.get(ids[i], col).map(str::to_string),
            &|_| None,
        )
        .unwrap_err();
        assert!(err.to_string().contains("fewer than two levels"));
    }

    #[test]
    fn collapse_genes_rules() {
        let mut counts = BTreeMap::new();
        counts.insert("B2".to_string(), 4usize);
        counts.insert("A1".to_string(), 2usize);
        counts.insert("C3".to_string(), 4usize);
        assert_eq!(CollapseGenes::None.representative(&counts).unwrap(), "A1");
        assert_eq!(
            CollapseGenes::MaxObserved.representative(&counts).unwrap(),
            "B2"
        );
        let values = vec![("A1".to_string(), 1.0), ("B2".to_string(), 3.0)];
        assert_eq!(CollapseGenes::Mean.collapse(&values, "A1"), Some(2.0));
        assert_eq!(CollapseGenes::None.collapse(&values, "A1"), Some(1.0));
        assert_eq!(
            CollapseGenes::MaxObserved.collapse(&values, "B2"),
            Some(3.0)
        );
        assert_eq!(CollapseGenes::None.collapse(&values, "Z9"), None);
    }

    #[test]
    fn predicates_parse_and_match() {
        let preds = parse_predicates("diagnosis_group!=MS; cohort == ms").unwrap();
        assert_eq!(preds.len(), 2);
        let get = |col: &str| match col {
            "diagnosis_group" => Some("Headache".to_string()),
            "cohort" => Some("ms".to_string()),
            _ => None,
        };
        assert!(predicates_match(&preds, &get));
        let get_ms = |col: &str| match col {
            "diagnosis_group" => Some("MS".to_string()),
            "cohort" => Some("ms".to_string()),
            _ => None,
        };
        assert!(!predicates_match(&preds, &get_ms));
        assert!(predicates_match(
            &parse_predicates("x!=1").unwrap(),
            &|_| None
        ));
        assert!(!predicates_match(
            &parse_predicates("x=1").unwrap(),
            &|_| None
        ));
        assert!(parse_predicates("nonsense").is_err());
    }
}
