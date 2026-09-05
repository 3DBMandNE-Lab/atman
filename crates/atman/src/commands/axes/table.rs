//! Wide per-subject score table: `sample_id`, `cohort`, `condition`,
//! `is_control`, then any number of numeric score columns. Non-numeric
//! extra columns are kept as text so they can drive `--group-by` and
//! subset predicates.

use anyhow::{anyhow, bail, Context, Result};
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Debug, Clone, Default)]
pub struct ScoreTable {
    pub sample_id: Vec<String>,
    pub cohort: Vec<String>,
    pub condition: Vec<String>,
    pub is_control: Vec<bool>,
    /// Numeric columns in file order.
    pub numeric_order: Vec<String>,
    pub numeric: BTreeMap<String, Vec<Option<f64>>>,
    pub text: BTreeMap<String, Vec<Option<String>>>,
}

fn parse_flag(value: &str, path: &Path) -> Result<bool> {
    match value.trim().to_ascii_lowercase().as_str() {
        "1" | "true" => Ok(true),
        "0" | "false" => Ok(false),
        other => bail!("is_control value {:?} in {:?} is not 0/1", other, path),
    }
}

impl ScoreTable {
    pub fn read(path: &Path) -> Result<ScoreTable> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .flexible(true)
            .from_path(path)
            .with_context(|| format!("opening {:?}", path))?;
        let headers: Vec<String> = reader
            .headers()?
            .iter()
            .map(|h| h.trim().to_string())
            .collect();
        let col = |name: &str| -> Result<usize> {
            headers
                .iter()
                .position(|h| h == name)
                .ok_or_else(|| anyhow!("missing column {:?} in {:?}", name, path))
        };
        let c_sample = col("sample_id")?;
        let c_cohort = col("cohort")?;
        let c_condition = col("condition")?;
        let c_control = col("is_control")?;
        let key_cols = [c_sample, c_cohort, c_condition, c_control];

        let mut table = ScoreTable::default();
        let mut raw: Vec<Vec<String>> = vec![Vec::new(); headers.len()];
        for row in reader.records() {
            let row = row.with_context(|| format!("reading {:?}", path))?;
            let sid = row.get(c_sample).unwrap_or("").trim().to_string();
            if sid.is_empty() {
                continue;
            }
            table.sample_id.push(sid);
            table
                .cohort
                .push(row.get(c_cohort).unwrap_or("").trim().to_string());
            table
                .condition
                .push(row.get(c_condition).unwrap_or("").trim().to_string());
            table
                .is_control
                .push(parse_flag(row.get(c_control).unwrap_or(""), path)?);
            for (i, bucket) in raw.iter_mut().enumerate() {
                if key_cols.contains(&i) {
                    continue;
                }
                bucket.push(row.get(i).unwrap_or("").trim().to_string());
            }
        }
        for (i, name) in headers.iter().enumerate() {
            if key_cols.contains(&i) || name.is_empty() {
                continue;
            }
            let values = &raw[i];
            let numeric = values
                .iter()
                .filter(|v| !v.is_empty())
                .all(|v| v.parse::<f64>().is_ok());
            if numeric {
                table.numeric_order.push(name.clone());
                table.numeric.insert(
                    name.clone(),
                    values
                        .iter()
                        .map(|v| {
                            if v.is_empty() {
                                None
                            } else {
                                v.parse::<f64>().ok().filter(|x| x.is_finite())
                            }
                        })
                        .collect(),
                );
            } else {
                table.text.insert(
                    name.clone(),
                    values
                        .iter()
                        .map(|v| if v.is_empty() { None } else { Some(v.clone()) })
                        .collect(),
                );
            }
        }
        if table.sample_id.is_empty() {
            bail!("no rows in {:?}", path);
        }
        Ok(table)
    }

    pub fn len(&self) -> usize {
        self.sample_id.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sample_id.is_empty()
    }

    pub fn numeric_column(&self, name: &str) -> Result<&[Option<f64>]> {
        self.numeric
            .get(name)
            .map(Vec::as_slice)
            .ok_or_else(|| anyhow!("score column {:?} not found or not numeric", name))
    }

    /// Text view of any column for row `i` (key columns, text columns, and
    /// numeric columns rendered with `{}`).
    pub fn text_value(&self, i: usize, col: &str) -> Option<String> {
        match col {
            "sample_id" => Some(self.sample_id[i].clone()),
            "cohort" => Some(self.cohort[i].clone()),
            "condition" => Some(self.condition[i].clone()),
            "is_control" => Some(if self.is_control[i] { "1" } else { "0" }.to_string()),
            _ => {
                if let Some(v) = self.text.get(col) {
                    return v[i].clone();
                }
                self.numeric
                    .get(col)
                    .and_then(|v| v[i].map(|x| format!("{x}")))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn score_table_classifies_columns() {
        let tmp = tempfile::tempdir().unwrap();
        let p = tmp.path().join("scores.tsv");
        std::fs::write(
            &p,
            "sample_id\tcohort\tcondition\tis_control\taxis1_z\tnote\tumap_1\nS1\tsih\tLeak\t0\t0.5\tfoo\t1.0\nS2\tsih\tNoLeak\t1\t\tbar\t2.0\n",
        )
        .unwrap();
        let t = ScoreTable::read(&p).unwrap();
        assert_eq!(t.len(), 2);
        assert!(!t.is_empty());
        assert_eq!(t.numeric_order, vec!["axis1_z", "umap_1"]);
        assert_eq!(t.numeric_column("axis1_z").unwrap(), &[Some(0.5), None]);
        assert!(t.numeric_column("note").is_err());
        assert_eq!(t.text_value(0, "note"), Some("foo".to_string()));
        assert_eq!(t.text_value(1, "is_control"), Some("1".to_string()));
        assert_eq!(t.text_value(1, "condition"), Some("NoLeak".to_string()));
        assert!(t.is_control[1]);
        assert!(!t.is_control[0]);
    }

    #[test]
    fn score_table_requires_key_columns() {
        let tmp = tempfile::tempdir().unwrap();
        let p = tmp.path().join("scores.tsv");
        std::fs::write(&p, "sample_id\tcohort\taxis1_z\nS1\tsih\t0.5\n").unwrap();
        assert!(ScoreTable::read(&p).is_err());
    }
}
