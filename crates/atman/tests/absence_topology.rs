use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

/// Build a fixture with `n_clusters` protein cliques, each containing
/// `proteins_per_clique` proteins that all fail in the same set of
/// `samples_per_clique` samples.  In addition, `n_unique_fragile` proteins
/// each fail in their own random subset of samples (so they should not
/// cluster), and `n_complete` proteins are observed in every sample.
fn write_fixture(
    dir: &std::path::Path,
    n_samples: usize,
    n_clusters: usize,
    proteins_per_clique: usize,
    samples_per_clique: usize,
    n_unique_fragile: usize,
    n_complete: usize,
) {
    std::fs::create_dir_all(dir).unwrap();
    let mut samples_buf =
        String::from("sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n");
    for s in 0..n_samples {
        samples_buf.push_str(&format!(
            "S{:04}\tsubj_{}\tcase\t0\ttumor\t{}\n",
            s,
            s,
            s + 1
        ));
    }
    std::fs::write(dir.join("samples.tsv"), samples_buf).unwrap();

    let total_proteins = n_clusters * proteins_per_clique + n_unique_fragile + n_complete;
    let mut proteins_buf =
        String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for i in 0..total_proteins {
        proteins_buf.push_str(&format!("diann_report\tA{:05}\t\tG{:05}\tms\t\n", i, i));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins_buf).unwrap();

    // Build per-protein absent sample sets.
    let mut absent_samples_for: Vec<std::collections::HashSet<usize>> =
        vec![std::collections::HashSet::new(); total_proteins];
    let mut prot = 0usize;
    for c in 0..n_clusters {
        let lo = c * samples_per_clique;
        let hi = lo + samples_per_clique;
        for _ in 0..proteins_per_clique {
            for s in lo..hi {
                absent_samples_for[prot].insert(s);
            }
            prot += 1;
        }
    }
    // Unique-fragile: each gets a random-looking but deterministic subset of samples.
    let mut rng_state: u64 = 0xCAFEBABEDEADBEEF;
    for _ in 0..n_unique_fragile {
        for _ in 0..10 {
            rng_state = rng_state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let s = (rng_state >> 32) as usize % n_samples;
            absent_samples_for[prot].insert(s);
        }
        prot += 1;
    }
    // Complete proteins: empty absence set; nothing to do.

    let mut buf = String::from("platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\tabundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\tdropped_by_qc\tplate_id\tpanel_lot\tingest_order\n");
    let mut row_order = 1u64;
    for s in 0..n_samples {
        for p in 0..total_proteins {
            if absent_samples_for[p].contains(&s) {
                continue;
            }
            buf.push_str(&format!(
                "diann_report\tS{:04}\tA{:05}\tG{:05}\tms\t\t10.0\t10.0\tlog2_intensity\tPASS\tPASS\t\t0\t0\t\t\t{}\n",
                s, p, p, row_order
            ));
            row_order += 1;
        }
    }
    std::fs::write(dir.join("measurements.tsv"), buf).unwrap();
}

fn parse_quality(path: &std::path::Path) -> std::collections::HashMap<String, String> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    let row = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    header
        .iter()
        .zip(row.iter())
        .map(|(h, v)| (h.to_string(), v.to_string()))
        .collect()
}

#[test]
fn absence_topology_recovers_fragility_cliques() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output = tmp.path().join("topology.tsv");
    // 60 samples, 4 cliques of 5 proteins each failing in 10 samples,
    // 10 unique-fragile proteins, 50 complete proteins.
    write_fixture(&input, 60, 4, 5, 10, 10, 50);

    let res = run_atman(&[
        "absence-topology",
        "--input-dir",
        input.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
        "--min-absent",
        "5",
    ]);
    assert!(
        res.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&res.stderr)
    );

    let q = parse_quality(&output.with_file_name("topology_quality.tsv"));
    let n_complete: usize = q["n_complete"].parse().unwrap();
    let n_clusterable: usize = q["n_clusterable"].parse().unwrap();
    // 4 cliques × 5 + 10 unique-fragile = 30 clusterable; 50 complete proteins.
    assert_eq!(n_complete, 50);
    assert_eq!(n_clusterable, 30);

    // 4 distinct cliques, each ending up as its own cluster of 5 proteins.
    let body = std::fs::read_to_string(&output).unwrap();
    let mut clique_sizes: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();
    for line in body.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols[2] != "complete" {
            *clique_sizes.entry(cols[2].to_string()).or_default() += 1;
        }
    }
    let big_clusters: Vec<usize> = clique_sizes.values().copied().filter(|&s| s == 5).collect();
    assert!(
        big_clusters.len() >= 4,
        "expected at least 4 cliques of size 5, got sizes: {:?}",
        clique_sizes
    );

    let sidecar: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(output.with_extension("tsv.run.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(sidecar["command"], "absence-topology");
    assert!(sidecar["reinvoke"]
        .as_str()
        .unwrap()
        .starts_with("atman absence-topology "));
}

#[test]
fn absence_topology_flags_random_absence() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input");
    let output = tmp.path().join("topology.tsv");
    // No cliques; every "fragile" protein gets a random absence set.
    write_fixture(&input, 60, 0, 0, 0, 50, 50);

    let res = run_atman(&[
        "absence-topology",
        "--input-dir",
        input.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
        "--min-absent",
        "3",
    ]);
    assert!(
        res.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&res.stderr)
    );

    let q = parse_quality(&output.with_file_name("topology_quality.tsv"));
    let diagnostic = &q["diagnostic"];
    assert_ne!(
        diagnostic, "protein_coherent",
        "random absence should not pass the protein_coherent gate (got {diagnostic})"
    );
}
