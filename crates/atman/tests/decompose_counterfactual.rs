use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

#[test]
fn counterfactual_delta_equals_loading_times_activation() {
    let tmp = tempfile::tempdir().unwrap();

    // 2 programs, 3 proteins, 2 samples.
    // Loadings:
    //   program_01: P1=0.5, P2=-0.3, P3=0.0
    //   program_02: P1=0.0, P2=0.4,  P3=0.8
    std::fs::write(
        tmp.path().join("loadings.tsv"),
        "program\tassay_id\tgene_symbol\tloading\n\
         program_01\tP1\tGENE1\t0.5\n\
         program_01\tP2\tGENE2\t-0.3\n\
         program_01\tP3\tGENE3\t0.0\n\
         program_02\tP1\tGENE1\t0.0\n\
         program_02\tP2\tGENE2\t0.4\n\
         program_02\tP3\tGENE3\t0.8\n",
    )
    .unwrap();

    // Activations:
    //   sample_A: program_01=2.0, program_02=1.0
    //   sample_B: program_01=-1.0, program_02=3.0
    std::fs::write(
        tmp.path().join("activations.tsv"),
        "cohort\tsubject_id\tsample_id\tprogram\tactivation\n\
         c\tA\tsample_A\tprogram_01\t2.0\n\
         c\tA\tsample_A\tprogram_02\t1.0\n\
         c\tB\tsample_B\tprogram_01\t-1.0\n\
         c\tB\tsample_B\tprogram_02\t3.0\n",
    )
    .unwrap();

    let output = tmp.path().join("delta.tsv");
    let result = run_atman(&[
        "decompose",
        "counterfactual",
        "--loadings",
        tmp.path().join("loadings.tsv").to_str().unwrap(),
        "--activations",
        tmp.path().join("activations.tsv").to_str().unwrap(),
        "--set-program",
        "program_01",
        "--to",
        "0",
        "--output",
        output.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let text = std::fs::read_to_string(&output).unwrap();
    let lines: Vec<&str> = text.lines().collect();

    // Header: assay_id, gene_symbol, sample_A, sample_B
    assert_eq!(lines[0], "assay_id\tgene_symbol\tsample_A\tsample_B");

    // Delta = loading × (original_activation - 0) = loading × activation
    // P1/GENE1: loading=0.5, sample_A act=2.0 → delta=1.0, sample_B act=-1.0 → delta=-0.5
    let p1: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(p1[0], "P1");
    assert_eq!(p1[1], "GENE1");
    assert!((p1[2].parse::<f64>().unwrap() - 1.0).abs() < 1e-9);
    assert!((p1[3].parse::<f64>().unwrap() - (-0.5)).abs() < 1e-9);

    // P2/GENE2: loading=-0.3, sample_A → delta=-0.6, sample_B → delta=0.3
    let p2: Vec<&str> = lines[2].split('\t').collect();
    assert_eq!(p2[0], "P2");
    assert!((p2[2].parse::<f64>().unwrap() - (-0.6)).abs() < 1e-9);
    assert!((p2[3].parse::<f64>().unwrap() - 0.3).abs() < 1e-9);

    // P3/GENE3: loading=0.0 on program_01 → delta=0 for both samples
    let p3: Vec<&str> = lines[3].split('\t').collect();
    assert_eq!(p3[0], "P3");
    assert_eq!(p3[2], "0");
    assert_eq!(p3[3], "0");

    // Sidecar exists
    assert!(tmp.path().join("delta.tsv.run.json").exists());
}

#[test]
fn counterfactual_rejects_unknown_program() {
    let tmp = tempfile::tempdir().unwrap();
    std::fs::write(
        tmp.path().join("loadings.tsv"),
        "program\tassay_id\tgene_symbol\tloading\n\
         program_01\tP1\tG1\t1.0\n",
    )
    .unwrap();
    std::fs::write(
        tmp.path().join("activations.tsv"),
        "cohort\tsubject_id\tsample_id\tprogram\tactivation\n\
         c\tA\tS1\tprogram_01\t1.0\n",
    )
    .unwrap();

    let result = run_atman(&[
        "decompose",
        "counterfactual",
        "--loadings",
        tmp.path().join("loadings.tsv").to_str().unwrap(),
        "--activations",
        tmp.path().join("activations.tsv").to_str().unwrap(),
        "--set-program",
        "program_99",
        "--to",
        "0",
        "--output",
        tmp.path().join("out.tsv").to_str().unwrap(),
    ]);
    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("program_99") && stderr.contains("not found"),
        "should name the missing program: {stderr}"
    );
}
