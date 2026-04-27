use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn write_loadings(path: &std::path::Path, rows: &[(&str, &str, f64)]) {
    let mut buf = String::from("program\tgene_symbol\tloading\n");
    for (program, gene, loading) in rows {
        buf.push_str(&format!("{program}\t{gene}\t{loading}\n"));
    }
    std::fs::write(path, buf).unwrap();
}

fn write_annotations(path: &std::path::Path, rows: &[(&str, &str)]) {
    let mut buf = String::from("program\tcategory\n");
    for (program, category) in rows {
        buf.push_str(&format!("{program}\t{category}\n"));
    }
    std::fs::write(path, buf).unwrap();
}

#[test]
fn align_programs_builds_multi_cohort_archetypes() {
    let tmp = tempfile::tempdir().unwrap();
    let cohort_a = tmp.path().join("a_loadings.tsv");
    let cohort_b = tmp.path().join("b_loadings.tsv");
    let cohort_c = tmp.path().join("c_loadings.tsv");
    let out = tmp.path().join("archetypes.tsv");

    // Shared "astro" program expressed similarly across A/B/C; each has one
    // idiosyncratic program that shouldn't match across cohorts.
    write_loadings(
        &cohort_a,
        &[
            ("program_01", "AQP4", 0.90),
            ("program_01", "GFAP", 0.80),
            ("program_01", "S100B", 0.70),
            ("program_01", "GENE_X", 0.05),
            ("program_02", "SYN1", 0.85),
            ("program_02", "SYP", 0.80),
            ("program_02", "GENE_X", 0.10),
        ],
    );
    write_loadings(
        &cohort_b,
        &[
            ("program_01", "AQP4", 0.88),
            ("program_01", "GFAP", 0.82),
            ("program_01", "S100B", 0.65),
            ("program_01", "GENE_X", 0.02),
            ("program_02", "IGHG1", 0.90),
            ("program_02", "IGHM", 0.80),
            ("program_02", "GENE_X", 0.05),
        ],
    );
    write_loadings(
        &cohort_c,
        &[
            ("program_01", "AQP4", 0.85),
            ("program_01", "GFAP", 0.78),
            ("program_01", "S100B", 0.68),
            ("program_01", "GENE_X", 0.01),
        ],
    );

    let output = run_atman(&[
        "align",
        "programs",
        "--loadings",
        &format!(
            "{},{},{}",
            cohort_a.to_str().unwrap(),
            cohort_b.to_str().unwrap(),
            cohort_c.to_str().unwrap()
        ),
        "--cohorts",
        "A,B,C",
        "--metric",
        "cosine",
        "--top-n",
        "10",
        "--tau",
        "0.5",
        "--reciprocal-best",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.starts_with(
        "archetype_id\tcohort\tprogram\tcategory\tn_members\tn_cohorts\tis_singleton\n"
    ));
    let data_rows: Vec<Vec<&str>> = text
        .lines()
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect();
    let universal_rows: Vec<_> = data_rows.iter().filter(|row| row[5] == "3").collect();
    // The three astro programs must land in one archetype spanning all 3 cohorts.
    assert_eq!(universal_rows.len(), 3);
    let cohorts: std::collections::BTreeSet<&str> = universal_rows.iter().map(|r| r[1]).collect();
    assert_eq!(cohorts.len(), 3);
}

#[test]
fn align_programs_sweep_produces_grid() {
    let tmp = tempfile::tempdir().unwrap();
    let cohort_a = tmp.path().join("a_loadings.tsv");
    let cohort_b = tmp.path().join("b_loadings.tsv");
    let ann_a = tmp.path().join("a_annotations.tsv");
    let ann_b = tmp.path().join("b_annotations.tsv");
    let out_matrix = tmp.path().join("sensitivity.tsv");

    write_loadings(
        &cohort_a,
        &[
            ("program_01", "AQP4", 0.9),
            ("program_01", "GFAP", 0.8),
            ("program_02", "SYN1", 0.85),
            ("program_02", "SYP", 0.70),
        ],
    );
    write_loadings(
        &cohort_b,
        &[
            ("program_01", "AQP4", 0.88),
            ("program_01", "GFAP", 0.78),
            ("program_02", "IGHG1", 0.9),
            ("program_02", "IGHM", 0.8),
        ],
    );
    write_annotations(
        &ann_a,
        &[("program_01", "astrocyte"), ("program_02", "synaptic")],
    );
    write_annotations(
        &ann_b,
        &[("program_01", "astrocyte"), ("program_02", "humoral")],
    );

    let output = run_atman(&[
        "align",
        "programs",
        "--loadings",
        &format!(
            "{},{}",
            cohort_a.to_str().unwrap(),
            cohort_b.to_str().unwrap()
        ),
        "--annotations",
        &format!("{},{}", ann_a.to_str().unwrap(), ann_b.to_str().unwrap()),
        "--cohorts",
        "A,B",
        "--sweep",
        "--metrics",
        "jaccard:top_n=[2,3]:tau=[0.1,0.5],cosine:tau=[0.3,0.6]",
        "--reciprocal-best",
        "--compare-constrained-vs-unconstrained",
        "--output-matrix",
        out_matrix.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let text = std::fs::read_to_string(&out_matrix).unwrap();
    assert!(text.starts_with(
        "metric\ttop_n\ttau\treciprocal_best\tcategory_constraint\tn_archetypes\tn_multi_cohort\tn_universal\tcategory_recovery\n"
    ));
    let rows: Vec<&str> = text.lines().skip(1).collect();
    // Jaccard: 2 top_n * 2 tau = 4 cells; Cosine: 1 top_n (default) * 2 tau = 2 cells.
    // Each cell emits constrained + unconstrained = 2 rows. Total = (4+2)*2 = 12.
    assert_eq!(rows.len(), 12);
    assert!(rows.iter().any(|r| r.starts_with("jaccard\t2\t0.1\t1\t1")));
    assert!(rows.iter().any(|r| r.starts_with("cosine\t40\t0.3\t1\t0")));
}
