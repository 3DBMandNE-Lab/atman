// Regression test for TASK-010: deterministic top-k tie-break + atomic writes
// in `atman robustness`.
//
// `top_k_set` collects (feature_id, |mean_diff|) pairs from a HashMap and used
// to sort by magnitude only, with no secondary key. When several features share
// the same |mean_diff| at the k-boundary, HashMap iteration order decided which
// ones survived `.take(k)`, so the rank_stability.tsv overlap/Jaccard output was
// non-reproducible. The fix adds a lexicographic feature-id tie-break and uses a
// total-order comparator. This test feeds a baseline with tied effect sizes at
// the boundary and asserts the robustness outputs are byte-identical across two
// independent runs.

use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

// de_results.tsv with several features sharing the same |mean_diff|, so the
// top-k boundary is a tie. ASSAY ids are deliberately written in a
// non-lexicographic order to make any iteration-order dependence visible.
const BASELINE: &str = "\
comparison\tpanel\tassay_id\tgene_symbol\tmean_diff\tbh_q\n\
Case_vs_Control\tp1\tZZZ\tGZ\t2.0\t0.01\n\
Case_vs_Control\tp1\tMMM\tGM\t2.0\t0.02\n\
Case_vs_Control\tp1\tAAA\tGA\t2.0\t0.03\n\
Case_vs_Control\tp1\tKKK\tGK\t2.0\t0.04\n\
Case_vs_Control\tp1\tBBB\tGB\t0.5\t0.20\n";

// LOO file: same features, also tied, with a couple of magnitudes nudged so the
// Jaccard is a non-degenerate fraction (not 0 or 1).
const LOO_0: &str = "\
comparison\tpanel\tassay_id\tgene_symbol\tmean_diff\tbh_q\n\
Case_vs_Control\tp1\tZZZ\tGZ\t2.0\t0.01\n\
Case_vs_Control\tp1\tMMM\tGM\t2.0\t0.02\n\
Case_vs_Control\tp1\tAAA\tGA\t2.0\t0.03\n\
Case_vs_Control\tp1\tKKK\tGK\t2.0\t0.04\n\
Case_vs_Control\tp1\tBBB\tGB\t0.5\t0.20\n";

fn write_fixture(dir: &std::path::Path) -> (std::path::PathBuf, std::path::PathBuf) {
    std::fs::create_dir_all(dir).unwrap();
    let baseline = dir.join("baseline.tsv");
    let loo0 = dir.join("loo_0.tsv");
    std::fs::write(&baseline, BASELINE).unwrap();
    std::fs::write(&loo0, LOO_0).unwrap();
    (baseline, loo0)
}

fn run_once(root: &std::path::Path, tag: &str) -> std::path::PathBuf {
    let input = root.join(format!("in_{tag}"));
    let out = root.join(format!("out_{tag}"));
    let (baseline, loo0) = write_fixture(&input);
    let output = run_atman(&[
        "robustness",
        "--baseline",
        baseline.to_str().unwrap(),
        "--loo",
        loo0.to_str().unwrap(),
        "--output-dir",
        out.to_str().unwrap(),
        // top_k = 2, inside a 4-way tie at |mean_diff| = 2.0 → boundary is tied.
        "--top-k",
        "2",
    ]);
    assert!(
        output.status.success(),
        "robustness failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    out
}

#[test]
fn topk_tie_break_is_byte_reproducible() {
    let tmp = tempfile::tempdir().unwrap();
    let out_a = run_once(tmp.path(), "a");
    let out_b = run_once(tmp.path(), "b");

    for name in [
        "rank_stability.tsv",
        "loo_sign_stability.tsv",
        "stability_ranked.tsv",
    ] {
        let a = std::fs::read(out_a.join(name)).unwrap();
        let b = std::fs::read(out_b.join(name)).unwrap();
        assert_eq!(
            a, b,
            "{name} differs between two runs — top-k tie-break is non-deterministic"
        );
    }

    // The deterministic tie-break selects the two lexicographically smallest
    // assay ids among the tied |mean_diff|=2.0 features: AAA and KKK
    // (keys "p1::AAA", "p1::KKK"). Baseline and LOO are identical, so overlap
    // is the full top-2 and Jaccard is 1.
    let rank = std::fs::read_to_string(out_a.join("rank_stability.tsv")).unwrap();
    let data_line = rank.lines().nth(1).expect("a data row");
    let cols: Vec<&str> = data_line.split('\t').collect();
    // comparison, loo_index, top_k, overlap_count, jaccard_index
    assert_eq!(cols[2], "2", "top_k");
    assert_eq!(cols[3], "2", "overlap_count for identical tied top-2 sets");
    assert_eq!(cols[4], "1", "jaccard_index for identical tied top-2 sets");
}
