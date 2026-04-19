use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

#[test]
fn enrich_ora_uses_measured_universe_and_reports_overlap() {
    let tmp = tempfile::tempdir().unwrap();
    let de = tmp.path().join("de_results.tsv");
    let sets = tmp.path().join("gene_sets.tsv");
    let out = tmp.path().join("ora.tsv");
    std::fs::write(
        &de,
        "\
panel\tassay_id\tgene_symbol\tuniprot\tcomparison\tn_pairs\tmean_a\tmean_b\tmean_diff\tt\tdf\tp_value\tbh_q\tskip_reason\n\
p\tA\tG1\t\tCase-Control\t4\t\t\t\t\t\t0.001\t0.01\t\n\
p\tB\tG2\t\tCase-Control\t4\t\t\t\t\t\t0.002\t0.02\t\n\
p\tC\tG3\t\tCase-Control\t4\t\t\t\t\t\t0.5\t0.5\t\n\
p\tD\tG4\t\tCase-Control\t4\t\t\t\t\t\t0.6\t0.6\t\n\
p\tE\tG5\t\tOther-Control\t4\t\t\t\t\t\t0.001\t0.01\t\n",
    )
    .unwrap();
    std::fs::write(
        &sets,
        "\
set_name\tgene_symbol\n\
signal\tG1\n\
signal\tG2\n\
signal\tG3\n\
outside\tG1\n\
outside\tNOPE\n",
    )
    .unwrap();

    let output = run_atman(&[
        "enrich",
        "ora",
        "--de-results",
        de.to_str().unwrap(),
        "--gene-sets",
        sets.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
        "--comparison",
        "Case-Control",
        "--q",
        "0.05",
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.starts_with("set_name\tuniverse_size\thit_count"));
    assert!(text.contains("signal\t4\t2\t3\t2\t"));
    assert!(text.contains("\tG1,G2\n"));
    assert!(text.contains("outside\t4\t2\t1\t1\t"));
    assert!(!text.contains("G5"));
}

#[test]
fn enrich_ora_accepts_explicit_universe_file() {
    let tmp = tempfile::tempdir().unwrap();
    let de = tmp.path().join("de_results.tsv");
    let sets = tmp.path().join("gene_sets.tsv");
    let universe = tmp.path().join("universe.tsv");
    let out = tmp.path().join("ora.tsv");
    std::fs::write(
        &de,
        "\
panel\tassay_id\tgene_symbol\tuniprot\tcomparison\tn_pairs\tmean_a\tmean_b\tmean_diff\tt\tdf\tp_value\tbh_q\tskip_reason\n\
p\tA\tG1\t\tCase-Control\t4\t\t\t\t\t\t0.001\t0.01\t\n\
p\tB\tG2\t\tCase-Control\t4\t\t\t\t\t\t0.5\t0.5\t\n",
    )
    .unwrap();
    std::fs::write(
        &sets,
        "\
set_name\tgene_symbol\n\
set_a\tG1\n\
set_a\tG3\n",
    )
    .unwrap();
    std::fs::write(&universe, "gene_symbol\nG1\nG2\nG3\nG4\n").unwrap();

    let output = run_atman(&[
        "enrich",
        "ora",
        "--de-results",
        de.to_str().unwrap(),
        "--gene-sets",
        sets.to_str().unwrap(),
        "--universe",
        universe.to_str().unwrap(),
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(out).unwrap();
    assert!(text.contains("set_a\t4\t1\t2\t1\t"));
}
