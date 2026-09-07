//! The `harmonize fit` / `harmonize apply` contract, exercised through
//! the CLI.
//!
//! The core algebra is tested in `atman_core::harmonize`. What is tested
//! here is the layer between it and the user: the model file. That file
//! is the thing that travels — "a replay can prove the held-out cohort
//! was absent" is the reason it exists — and until now nothing tested
//! it end to end.
//!
//! Two of its fields carry the module's safety claims. `fit_cohorts` is
//! what `apply` refuses a training cohort with, and `permuted_labels` is
//! what marks the negative-control arm. Both were read with
//! `unwrap_or_default`, so a truncated or hand-edited file made the
//! refusal permit every cohort and made a permuted model report itself
//! as a real one. Both failures are silent, which is why they are
//! tested by mutating a real model file rather than by inspection.

use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

/// A cohort with a planted case/control contrast on the first half of
/// its features, so a fitted direction has something real to carry.
fn write_cohort(dir: &Path, tag: &str, n: usize, p: usize, effect: f64) {
    std::fs::create_dir_all(dir).unwrap();

    let subjects: Vec<String> = (0..n).map(|i| format!("{tag}_S{i:03}")).collect();
    let features: Vec<String> = (0..p).map(|j| format!("P{j:03}")).collect();

    let mut samples = String::from(
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\n",
    );
    for (i, s) in subjects.iter().enumerate() {
        let cond = if i % 2 == 0 { "Control" } else { "Case" };
        samples.push_str(&format!("{s}\t{s}\t{cond}\t0\tcsf\t{}\n", i + 1));
    }
    std::fs::write(dir.join("samples.tsv"), samples).unwrap();

    let mut proteins = String::from("platform\tassay_id\tuniprot\tgene_symbol\tpanel\tpanel_lot\n");
    for f in &features {
        proteins.push_str(&format!(
            "spectronaut_report\t{f}\t{f}\tG{f}\tspectronaut\t\n"
        ));
    }
    std::fs::write(dir.join("proteins.tsv"), proteins).unwrap();

    // Deterministic pseudo-noise: no rand dependency in the test tree.
    let mut state: u64 = 0x9E3779B97F4A7C15;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        ((state >> 11) as f64 / (1u64 << 53) as f64) - 0.5
    };

    let hdr = "platform\tsample_id\tassay_id\tgene_symbol\tpanel\tnpx_source_str\tabundance\t\
               abundance_raw\tabundance_unit\tqc_sample\tqc_assay\tdetection_limit\tbelow_lod\t\
               dropped_by_qc\tplate_id\tpanel_lot\tingest_order\n";
    let mut m = String::from(hdr);
    for (i, s) in subjects.iter().enumerate() {
        let is_case = i % 2 == 1;
        for (j, f) in features.iter().enumerate() {
            let signal = if is_case && j < p / 2 { effect } else { 0.0 };
            let a = 10.0 + signal + next();
            m.push_str(&format!(
                "spectronaut_report\t{s}\t{f}\tG{f}\tspectronaut\t{a:.4}\t{a:.4}\t{a:.4}\t\
                 log2_intensity\tPASS\tPASS\t\t0\t0\t\t\t1\n"
            ));
        }
    }
    std::fs::write(dir.join("measurements.tsv"), m).unwrap();
}

struct Fixture {
    _tmp: tempfile::TempDir,
    root: std::path::PathBuf,
}

impl Fixture {
    /// Two training cohorts and one held-out cohort.
    fn new() -> Self {
        let tmp = tempfile::tempdir().unwrap();
        let root = tmp.path().to_path_buf();
        write_cohort(&root.join("A"), "A", 24, 20, 1.5);
        write_cohort(&root.join("B"), "B", 24, 20, 1.5);
        write_cohort(&root.join("H"), "H", 24, 20, 1.5);
        Self { _tmp: tmp, root }
    }

    fn p(&self, rel: &str) -> String {
        self.root.join(rel).display().to_string()
    }

    fn fit(&self, method: &str, extra: &[&str]) -> Output {
        let cohorts = format!("A={},B={}", self.p("A"), self.p("B"));
        let model = self.p("model.json");
        let mut args = vec![
            "harmonize",
            "fit",
            "--cohorts",
            &cohorts,
            "--method",
            method,
            "--case-value",
            "Case",
            "--output-model",
            &model,
        ];
        args.extend_from_slice(extra);
        run_atman(&args)
    }

    fn apply(&self, cohort_spec: &str, extra: &[&str]) -> Output {
        let model = self.p("model.json");
        let out = self.p("scores.tsv");
        let mut args = vec![
            "harmonize",
            "apply",
            "--model",
            &model,
            "--cohort",
            cohort_spec,
            "--case-value",
            "Case",
            "--output",
            &out,
        ];
        args.extend_from_slice(extra);
        run_atman(&args)
    }

    fn model_json(&self) -> serde_json::Value {
        let text = std::fs::read_to_string(self.root.join("model.json")).unwrap();
        serde_json::from_str(&text).unwrap()
    }

    /// Rewrite the model file with `mutate` applied, to simulate a
    /// truncated, hand-edited, or foreign-produced artefact.
    fn mutate_model(&self, mutate: impl FnOnce(&mut serde_json::Map<String, serde_json::Value>)) {
        let mut doc = self.model_json();
        mutate(doc.as_object_mut().unwrap());
        std::fs::write(
            self.root.join("model.json"),
            serde_json::to_string_pretty(&doc).unwrap(),
        )
        .unwrap();
    }
}

fn stderr(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).to_string()
}

/// The contract works at all: fit writes a model, apply reads it and
/// scores a cohort that was not in the fit.
#[test]
fn fit_then_apply_round_trips_through_the_model_file() {
    let f = Fixture::new();
    let out = f.fit("zscore", &[]);
    assert!(out.status.success(), "fit failed: {}", stderr(&out));

    let doc = f.model_json();
    assert_eq!(doc["schema_version"], 1);
    assert_eq!(doc["method"], "zscore");
    assert_eq!(doc["permuted_labels"], false);
    assert_eq!(doc["fit_cohorts"], serde_json::json!(["A", "B"]));
    assert!(doc["fit_inputs_sha256"].is_object());

    let held = format!("H={}", f.p("H"));
    let out = f.apply(&held, &[]);
    assert!(out.status.success(), "apply failed: {}", stderr(&out));

    let scores = std::fs::read_to_string(f.root.join("scores.tsv")).unwrap();
    assert_eq!(
        scores.lines().count(),
        25,
        "expected a header plus 24 subjects"
    );
}

/// The leakage refusal fires on a cohort the model was fit on.
#[test]
fn apply_refuses_a_training_cohort_through_the_cli() {
    let f = Fixture::new();
    assert!(f.fit("zscore", &[]).status.success());

    let training = format!("A={}", f.p("A"));
    let out = f.apply(&training, &[]);
    assert!(
        !out.status.success(),
        "applying to a training cohort must fail"
    );
    let e = stderr(&out);
    assert!(
        e.contains("training cohorts") || e.contains("held-out"),
        "refusal must say why: {e}"
    );
}

/// An absent `fit_cohorts` must be refused, not read as an empty list.
///
/// Read as empty, the refusal in `apply_with` matches nothing and every
/// cohort is permitted — including the ones the model was fit on. The
/// run then succeeds and reports a held-out evaluation that is not one.
#[test]
fn a_model_without_fit_cohorts_is_refused_rather_than_failing_open() {
    let f = Fixture::new();
    assert!(f.fit("zscore", &[]).status.success());
    f.mutate_model(|m| {
        m.remove("fit_cohorts");
    });

    // The dangerous case: the cohort this model was actually fit on.
    let training = format!("A={}", f.p("A"));
    let out = f.apply(&training, &[]);
    assert!(
        !out.status.success(),
        "a model with no fit_cohorts silently permitted its own \
         training cohort; the leakage guard failed open"
    );
    assert!(
        stderr(&out).contains("fit_cohorts"),
        "the error must name the missing field: {}",
        stderr(&out)
    );
}

/// An empty `fit_cohorts` list fails open the same way as an absent one.
#[test]
fn a_model_with_an_empty_fit_cohorts_list_is_refused() {
    let f = Fixture::new();
    assert!(f.fit("zscore", &[]).status.success());
    f.mutate_model(|m| {
        m.insert("fit_cohorts".into(), serde_json::json!([]));
    });

    let training = format!("A={}", f.p("A"));
    let out = f.apply(&training, &[]);
    assert!(!out.status.success());
    assert!(stderr(&out).contains("no training cohorts"));
}

/// An absent `permuted_labels` must be refused, not read as `false`.
///
/// Read as false, a negative-control model reports itself as a
/// real-labels one in both the log line and the sidecar.
#[test]
fn a_model_without_permuted_labels_is_refused_rather_than_read_as_real() {
    let f = Fixture::new();
    assert!(f.fit("zscore", &["--permute-labels"]).status.success());
    assert_eq!(f.model_json()["permuted_labels"], true);

    f.mutate_model(|m| {
        m.remove("permuted_labels");
    });
    let held = format!("H={}", f.p("H"));
    let out = f.apply(&held, &[]);
    assert!(
        !out.status.success(),
        "a model with no permuted_labels was read as a real-labels model"
    );
    assert!(stderr(&out).contains("permuted_labels"));
}

/// The permuted flag survives the round trip and is reported on apply.
#[test]
fn the_permuted_arm_is_marked_through_the_model_file() {
    let f = Fixture::new();
    assert!(f.fit("zscore", &["--permute-labels"]).status.success());
    let held = format!("H={}", f.p("H"));
    let out = f.apply(&held, &[]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(
        stderr(&out).contains("permuted=true"),
        "apply must report that this is the negative-control arm: {}",
        stderr(&out)
    );
}

/// A per-sample method warns when `--center-direction` is off, because
/// the anchor term otherwise lets any direction separate the arms.
#[test]
fn apply_warns_when_center_direction_is_off_for_a_per_sample_method() {
    let f = Fixture::new();
    let out = f.fit("reference-protein", &["--reference-k", "4"]);
    assert!(out.status.success(), "{}", stderr(&out));

    let held = format!("H={}", f.p("H"));
    let out = f.apply(&held, &[]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(
        stderr(&out).contains("center-direction"),
        "a per-sample method must warn when centring is off: {}",
        stderr(&out)
    );

    // And it must not warn when the flag is on.
    let out = f.apply(&held, &["--center-direction"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(!stderr(&out).contains("warning:"));
}

/// A model naming a method whose fitted state is absent is refused
/// rather than scoring every subject as missing.
#[test]
fn a_reference_protein_model_without_its_reference_set_is_refused() {
    let f = Fixture::new();
    assert!(f
        .fit("reference-protein", &["--reference-k", "4"])
        .status
        .success());
    f.mutate_model(|m| {
        m.insert("reference_features".into(), serde_json::json!([]));
    });

    let held = format!("H={}", f.p("H"));
    let out = f.apply(&held, &[]);
    assert!(!out.status.success());
    assert!(stderr(&out).contains("empty reference set"));
}

/// Structural corruption is named rather than coerced to a default.
#[test]
fn a_malformed_direction_entry_is_refused() {
    let f = Fixture::new();
    assert!(f.fit("zscore", &[]).status.success());
    f.mutate_model(|m| {
        let mut d = m["direction"].as_array().unwrap().clone();
        d[0] = serde_json::json!("not-a-number");
        m.insert("direction".into(), serde_json::Value::Array(d));
    });

    let held = format!("H={}", f.p("H"));
    let out = f.apply(&held, &[]);
    assert!(!out.status.success());
    assert!(stderr(&out).contains("non-numeric"));
}

/// An unknown schema version is refused rather than read as version 1.
#[test]
fn a_future_schema_version_is_refused() {
    let f = Fixture::new();
    assert!(f.fit("zscore", &[]).status.success());
    f.mutate_model(|m| {
        m.insert("schema_version".into(), serde_json::json!(2));
    });

    let held = format!("H={}", f.p("H"));
    let out = f.apply(&held, &[]);
    assert!(!out.status.success());
    assert!(stderr(&out).contains("schema_version 2"));
}
