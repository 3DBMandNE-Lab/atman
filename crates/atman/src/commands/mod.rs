use anyhow::{bail, Result};

pub mod align;
pub mod asymmetry;
pub mod bench;
pub mod bootstrap;
pub mod coupling;
pub mod de;
pub mod decompose;
pub mod enrich;
pub mod enrich_gprofiler;
pub mod fold_change;
pub mod ingest;
pub mod ingest_matrix;
pub mod matrix;
pub mod meta;
pub mod module_de;
pub mod module_trajectory;
pub mod modules;
pub mod network;
pub mod null;
pub mod programs;
pub mod qc;
pub mod ratio;
pub mod report;
pub mod residuals;
pub mod robustness;
pub mod run;
pub mod score;
pub mod validate;
pub mod within_cohort_rank;

/// Parse comma-separated comparisons in `A-B` form with strict validation.
pub fn parse_comparisons(groups: &str) -> Result<Vec<(String, String)>> {
    let mut comparisons = Vec::new();
    for raw in groups.split(',') {
        let token = raw.trim();
        if token.is_empty() {
            bail!("invalid empty comparison in --groups {:?}", groups);
        }

        let mut parts = token.split('-').map(str::trim);
        let a = parts.next().unwrap_or_default();
        let b = parts.next().unwrap_or_default();
        if a.is_empty() || b.is_empty() || parts.next().is_some() {
            bail!(
                "invalid comparison {:?}; expected exactly one A-B pair",
                token
            );
        }
        comparisons.push((a.to_string(), b.to_string()));
    }

    if comparisons.is_empty() {
        bail!("no comparisons given");
    }
    Ok(comparisons)
}

#[cfg(test)]
mod tests {
    use super::parse_comparisons;

    #[test]
    fn parse_comparisons_accepts_valid_list() {
        let parsed = parse_comparisons("PT2-PT1,PR2-PR1").expect("parse comparisons");
        assert_eq!(
            parsed,
            vec![
                ("PT2".to_string(), "PT1".to_string()),
                ("PR2".to_string(), "PR1".to_string())
            ]
        );
    }

    #[test]
    fn parse_comparisons_rejects_empty_tokens() {
        let err = parse_comparisons("PT2-PT1,").expect_err("must fail on trailing comma");
        assert!(err.to_string().contains("invalid empty comparison"));
    }

    #[test]
    fn parse_comparisons_rejects_malformed_pair() {
        let err = parse_comparisons("PT2").expect_err("must fail without separator");
        assert!(err.to_string().contains("expected exactly one A-B pair"));
    }
}
