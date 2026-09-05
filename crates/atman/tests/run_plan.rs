use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

#[test]
fn run_plan_executes_stages_and_writes_manifest() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = tmp.path().join("work");
    std::fs::create_dir(&input_dir).unwrap();
    std::fs::write(input_dir.join("seed.txt"), "hello\n").unwrap();

    let plan = tmp.path().join("plan.yaml");
    std::fs::write(
        &plan,
        r#"
name: two_stage_demo
plan_commit: sha-test-001
stages:
  - id: copy_seed
    command: cp seed.txt staged.txt
    inputs:
      - seed.txt
    outputs:
      - staged.txt
  - id: append_tag
    command: printf 'tag\n' >> staged.txt
    inputs:
      - staged.txt
    outputs:
      - staged.txt
"#,
    )
    .unwrap();

    let output_dir = tmp.path().join("manifest");
    let result = run_atman(&[
        "run",
        "--plan",
        plan.to_str().unwrap(),
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let manifest_path = output_dir.join("plan_manifest.tsv");
    let manifest = std::fs::read_to_string(&manifest_path).unwrap();
    assert!(manifest.starts_with(
        "plan_name\tplan_commit\tplan_hash\tstage_id\tcommand\tinput_hash\toutput_hash\truntime_s\texit_code\tatman_version\tsystem\tstarted_at_unix_s\tsidecar_hash\n"
    ));
    let rows: Vec<Vec<&str>> = manifest
        .lines()
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0][0], "two_stage_demo");
    assert_eq!(rows[0][1], "sha-test-001");
    assert_eq!(rows[0][3], "copy_seed");
    assert_eq!(rows[0][8], "0"); // exit_code

    // Plan hash matches byte-for-byte across both rows (same plan).
    assert_eq!(rows[0][2], rows[1][2]);

    // Verify the staged file was actually produced.
    let staged = std::fs::read_to_string(input_dir.join("staged.txt")).unwrap();
    assert_eq!(staged, "hello\ntag\n");
}

#[test]
fn run_plan_detects_drift_without_plan_commit_bump() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = tmp.path().join("work");
    std::fs::create_dir(&input_dir).unwrap();

    let plan = tmp.path().join("plan.yaml");
    let v1 = r#"
name: drift_demo
plan_commit: sha-001
stages:
  - id: noop
    command: "true"
"#;
    std::fs::write(&plan, v1).unwrap();

    let output_dir = tmp.path().join("manifest");
    let first = run_atman(&[
        "run",
        "--plan",
        plan.to_str().unwrap(),
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
    ]);
    assert!(first.status.success());

    // Edit the plan but keep plan_commit the same.
    let v2 = r#"
name: drift_demo
plan_commit: sha-001
stages:
  - id: noop
    command: "false"
"#;
    std::fs::write(&plan, v2).unwrap();

    let drift = run_atman(&[
        "run",
        "--plan",
        plan.to_str().unwrap(),
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
    ]);
    assert!(!drift.status.success());
    let stderr = String::from_utf8_lossy(&drift.stderr);
    assert!(
        stderr.contains("plan content drift"),
        "expected drift warning in stderr: {stderr}"
    );

    // Passing --allow-drift should permit the re-run.
    let allowed = run_atman(&[
        "run",
        "--plan",
        plan.to_str().unwrap(),
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
        "--allow-drift",
    ]);
    let allowed_stderr = String::from_utf8_lossy(&allowed.stderr);
    // --allow-drift must suppress the drift guard: the run gets past
    // consistency checking and into stage execution.
    assert!(
        !allowed_stderr.contains("plan content drift"),
        "--allow-drift should suppress the drift error, got: {allowed_stderr}"
    );
    // The (only) stage runs `false`, which exits 1; without
    // --continue-on-error the run must surface that stage failure.
    assert!(
        !allowed.status.success(),
        "run should fail because stage `noop` (command `false`) exits non-zero"
    );
    assert!(
        allowed_stderr.contains("stage noop exited")
            || allowed_stderr.contains("one or more stages failed"),
        "expected a stage-failure message, got: {allowed_stderr}"
    );
}

#[test]
fn run_plan_propagates_stage_failure() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = tmp.path().join("work");
    std::fs::create_dir(&input_dir).unwrap();
    let plan = tmp.path().join("plan.yaml");
    std::fs::write(
        &plan,
        r#"
name: fail_demo
stages:
  - id: boom
    command: "exit 3"
"#,
    )
    .unwrap();
    let output_dir = tmp.path().join("manifest");

    let result = run_atman(&[
        "run",
        "--plan",
        plan.to_str().unwrap(),
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
    ]);
    assert!(!result.status.success());

    let manifest = std::fs::read_to_string(output_dir.join("plan_manifest.tsv")).unwrap();
    let row = manifest.lines().nth(1).unwrap();
    let cols: Vec<&str> = row.split('\t').collect();
    assert_eq!(cols[3], "boom");
    assert_eq!(cols[8], "3");
}

#[test]
fn run_plan_dry_run_validates_inputs_without_executing() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = tmp.path().join("work");
    std::fs::create_dir(&input_dir).unwrap();
    std::fs::write(input_dir.join("seed.txt"), "hello\n").unwrap();
    let plan = tmp.path().join("plan.yaml");
    std::fs::write(
        &plan,
        r#"
name: dry
vars:
  out: results
stages:
  - id: a
    command: mkdir -p ${out} && cp seed.txt ${out}/staged.txt
    inputs: [seed.txt]
    outputs: ["${out}/staged.txt"]
  - id: b
    command: cp ${out}/staged.txt ${out}/final.txt
    inputs: ["${out}/staged.txt"]
    outputs: ["${out}/final.txt"]
"#,
    )
    .unwrap();
    let output_dir = tmp.path().join("manifest");
    let r = run_atman(&[
        "run",
        "--plan",
        plan.to_str().unwrap(),
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
        "--dry-run",
    ]);
    let stderr = String::from_utf8_lossy(&r.stderr);
    assert!(r.status.success(), "stderr:\n{stderr}");
    assert!(
        stderr.contains("results/staged.txt  (produced by an earlier stage)"),
        "{stderr}"
    );
    assert!(stderr.contains("dry run ok"));
    assert!(
        !input_dir.join("results").exists(),
        "dry run must not execute"
    );
    assert!(!output_dir.join("plan_manifest.tsv").exists());

    // A missing input that no earlier stage produces is rejected.
    std::fs::write(
        &plan,
        r#"
name: dry
stages:
  - id: a
    command: cat nope.txt
    inputs: [nope.txt]
"#,
    )
    .unwrap();
    let r = run_atman(&[
        "run",
        "--plan",
        plan.to_str().unwrap(),
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
        "--dry-run",
    ]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("nope.txt"));

    // An undefined variable is rejected at parse time.
    std::fs::write(
        &plan,
        "name: dry\nstages:\n  - id: a\n    command: echo ${missing}\n",
    )
    .unwrap();
    let r = run_atman(&[
        "run",
        "--plan",
        plan.to_str().unwrap(),
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
        "--dry-run",
    ]);
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("undefined variable"));
}

#[test]
fn run_plan_hashes_sidecars_and_enforces_strict_outputs() {
    let tmp = tempfile::tempdir().unwrap();
    let input_dir = tmp.path().join("work");
    std::fs::create_dir(&input_dir).unwrap();
    let plan = tmp.path().join("plan.yaml");
    std::fs::write(
        &plan,
        r#"
name: sidecars
stages:
  - id: with_sidecar
    command: printf 'x\n' > out.tsv && printf '{}\n' > out.tsv.run.json
    outputs: [out.tsv]
  - id: forgets_output
    command: "true"
    outputs: [never.tsv]
"#,
    )
    .unwrap();
    let output_dir = tmp.path().join("manifest");
    let r = run_atman(&[
        "run",
        "--plan",
        plan.to_str().unwrap(),
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--output-dir",
        output_dir.to_str().unwrap(),
    ]);
    assert!(
        !r.status.success(),
        "strict outputs must fail the second stage"
    );
    assert!(String::from_utf8_lossy(&r.stderr).contains("never.tsv"));
    let manifest = std::fs::read_to_string(output_dir.join("plan_manifest.tsv")).unwrap();
    let rows: Vec<Vec<&str>> = manifest
        .lines()
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0][3], "with_sidecar");
    assert_eq!(rows[0][8], "0");
    assert_eq!(
        rows[0][12].len(),
        64,
        "sidecar_hash populated: {:?}",
        rows[0][12]
    );
    assert_eq!(rows[1][3], "forgets_output");
    assert_eq!(rows[1][8], "2");
    assert!(rows[1][12].is_empty());

    // --strict-outputs false restores the permissive behaviour.
    let output_dir2 = tmp.path().join("manifest2");
    let r = run_atman(&[
        "run",
        "--plan",
        plan.to_str().unwrap(),
        "--input-dir",
        input_dir.to_str().unwrap(),
        "--output-dir",
        output_dir2.to_str().unwrap(),
        "--strict-outputs",
        "false",
    ]);
    assert!(
        r.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&r.stderr)
    );
}
