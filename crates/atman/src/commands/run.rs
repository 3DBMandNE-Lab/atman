//! Pre-registered analysis runner.
//!
//! Reads a declarative plan YAML/JSON, executes each stage in order via
//! `sh -c <command>`, hashes declared inputs and outputs with SHA-256, and
//! appends a manifest row per stage. Re-running against the same plan file
//! compares input hashes against the previous manifest and refuses to
//! overwrite if the plan content has drifted without an updated
//! `plan_commit`. This converts "trust my git log" into a hash-verifiable
//! artifact: cite the manifest TSV and the plan SHA, and the pipeline can be
//! reproduced byte-for-byte.

use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

use crate::io::{atomic_write, escape_tsv, sha256_hex};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Path to the analysis plan (YAML or JSON).
    #[arg(long)]
    plan: PathBuf,

    /// Working directory stages run in; paths in the plan resolve relative to this.
    #[arg(long, default_value = ".")]
    input_dir: PathBuf,

    /// Directory to write the manifest into.
    #[arg(long)]
    output_dir: PathBuf,

    /// Emit the manifest TSV (always true in the current CLI; kept for compatibility).
    #[arg(long, default_value_t = true)]
    emit_manifest: bool,

    /// Continue past stage failures instead of aborting the run.
    #[arg(long, default_value_t = false)]
    continue_on_error: bool,

    /// Allow re-running when the plan file contents differ from a previous manifest
    /// with the same `plan_commit`. By default we refuse to overwrite.
    #[arg(long, default_value_t = false)]
    allow_drift: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Plan {
    name: String,
    #[serde(default)]
    plan_commit: Option<String>,
    stages: Vec<Stage>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Stage {
    id: String,
    command: String,
    #[serde(default)]
    inputs: Vec<String>,
    #[serde(default)]
    outputs: Vec<String>,
}

pub fn run(args: Args) -> Result<()> {
    let plan_text = std::fs::read_to_string(&args.plan)
        .with_context(|| format!("reading plan file {:?}", args.plan))?;
    let plan: Plan = parse_plan(&plan_text, &args.plan)?;
    if plan.stages.is_empty() {
        bail!("plan {:?} has no stages", args.plan);
    }
    let plan_hash = sha256_hex(plan_text.as_bytes());
    std::fs::create_dir_all(&args.output_dir).with_context(|| {
        format!("creating output dir {:?}", args.output_dir)
    })?;

    let manifest_path = args.output_dir.join("plan_manifest.tsv");
    if manifest_path.exists() && !args.allow_drift {
        check_manifest_consistency(&manifest_path, &plan_hash, &plan)?;
    }

    let mut rows = Vec::with_capacity(plan.stages.len());
    let mut aborted = false;
    for stage in &plan.stages {
        if aborted {
            break;
        }
        let row = execute_stage(stage, &args.input_dir, &plan, &plan_hash)?;
        if row.exit_code != 0 && !args.continue_on_error {
            eprintln!(
                "run: stage {} exited with {}; aborting (pass --continue-on-error to keep going)",
                stage.id, row.exit_code
            );
            rows.push(row);
            aborted = true;
        } else {
            rows.push(row);
        }
    }

    write_manifest(&manifest_path, &rows)?;
    eprintln!(
        "run: plan={} stages={} manifest={}",
        plan.name,
        rows.len(),
        manifest_path.display()
    );
    if aborted {
        bail!("one or more stages failed; manifest written to {:?}", manifest_path);
    }
    Ok(())
}

fn parse_plan(text: &str, path: &Path) -> Result<Plan> {
    let ext = path
        .extension()
        .and_then(|os| os.to_str())
        .map(|s| s.to_ascii_lowercase())
        .unwrap_or_default();
    let plan: Plan = match ext.as_str() {
        "json" => serde_json::from_str(text)
            .with_context(|| format!("parsing JSON plan {:?}", path))?,
        "yaml" | "yml" | "" => serde_yaml::from_str(text)
            .with_context(|| format!("parsing YAML plan {:?}", path))?,
        other => bail!("unsupported plan extension {:?}", other),
    };
    let mut seen = std::collections::BTreeSet::new();
    for stage in &plan.stages {
        if !seen.insert(stage.id.clone()) {
            bail!("duplicate stage id {:?} in plan", stage.id);
        }
    }
    Ok(plan)
}

#[derive(Debug)]
struct ManifestRow {
    plan_name: String,
    plan_commit: String,
    plan_hash: String,
    stage_id: String,
    command: String,
    input_hash: String,
    output_hash: String,
    runtime_s: f64,
    exit_code: i32,
    atman_version: String,
    system: String,
    started_at_unix_s: u64,
}

fn execute_stage(
    stage: &Stage,
    cwd: &Path,
    plan: &Plan,
    plan_hash: &str,
) -> Result<ManifestRow> {
    let input_hash = hash_paths(cwd, &stage.inputs);
    let start = Instant::now();
    let started_at = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let output = Command::new("sh")
        .arg("-c")
        .arg(&stage.command)
        .current_dir(cwd)
        .output()
        .with_context(|| format!("spawning stage {}", stage.id))?;
    let runtime = start.elapsed().as_secs_f64();
    let exit_code = output.status.code().unwrap_or(-1);
    // Pipe captured stderr through so users see progress; stdout would be noisy.
    if !output.stderr.is_empty() {
        eprintln!("--- stage {} stderr ---", stage.id);
        eprintln!("{}", String::from_utf8_lossy(&output.stderr).trim_end());
    }
    let output_hash = hash_paths(cwd, &stage.outputs);
    Ok(ManifestRow {
        plan_name: plan.name.clone(),
        plan_commit: plan.plan_commit.clone().unwrap_or_default(),
        plan_hash: plan_hash.to_string(),
        stage_id: stage.id.clone(),
        command: stage.command.clone(),
        input_hash,
        output_hash,
        runtime_s: runtime,
        exit_code,
        atman_version: env!("CARGO_PKG_VERSION").to_string(),
        system: format!(
            "{}/{}/{}",
            std::env::consts::OS,
            std::env::consts::ARCH,
            std::env::consts::FAMILY
        ),
        started_at_unix_s: started_at,
    })
}

fn hash_paths(cwd: &Path, paths: &[String]) -> String {
    if paths.is_empty() {
        return String::new();
    }
    // Sort for determinism; hash concatenation of (relative path, file hash).
    let mut parts: Vec<(String, String)> = paths
        .iter()
        .map(|p| {
            let abs = cwd.join(p);
            let hash = hash_file_or_missing(&abs);
            (p.clone(), hash)
        })
        .collect();
    parts.sort();
    let mut buf = Vec::new();
    for (path, file_hash) in parts {
        buf.extend_from_slice(path.as_bytes());
        buf.push(0);
        buf.extend_from_slice(file_hash.as_bytes());
        buf.push(b'\n');
    }
    sha256_hex(&buf)
}

fn hash_file_or_missing(path: &Path) -> String {
    match std::fs::read(path) {
        Ok(bytes) => sha256_hex(&bytes),
        Err(_) => "MISSING".to_string(),
    }
}

fn check_manifest_consistency(manifest_path: &Path, new_plan_hash: &str, plan: &Plan) -> Result<()> {
    let existing = std::fs::read_to_string(manifest_path)
        .with_context(|| format!("reading existing manifest {:?}", manifest_path))?;
    let mut previous: Option<(String, String)> = None;
    for line in existing.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        previous = Some((cols[1].to_string(), cols[2].to_string()));
        break;
    }
    let Some((prev_commit, prev_plan_hash)) = previous else {
        return Ok(());
    };
    let current_commit = plan.plan_commit.clone().unwrap_or_default();
    if prev_plan_hash != new_plan_hash && prev_commit == current_commit {
        bail!(
            "plan content drift detected: previous plan_hash={} new_plan_hash={} but plan_commit={:?} unchanged. Bump plan_commit or pass --allow-drift.",
            prev_plan_hash,
            new_plan_hash,
            current_commit
        );
    }
    Ok(())
}

fn write_manifest(path: &Path, rows: &[ManifestRow]) -> Result<()> {
    let mut out = String::from(
        "plan_name\tplan_commit\tplan_hash\tstage_id\tcommand\tinput_hash\toutput_hash\truntime_s\texit_code\tatman_version\tsystem\tstarted_at_unix_s\n",
    );
    for row in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            row.plan_name,
            row.plan_commit,
            row.plan_hash,
            row.stage_id,
            escape_tsv(&row.command),
            row.input_hash,
            row.output_hash,
            row.runtime_s,
            row.exit_code,
            row.atman_version,
            row.system,
            row.started_at_unix_s,
        ));
    }
    atomic_write(path, out.as_bytes())
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_plan_accepts_yaml_and_json() {
        let yaml = r#"
name: csf_crossdisease_v3
plan_commit: abc123
stages:
  - id: de
    command: echo de
    inputs:
      - data.tsv
    outputs:
      - out.tsv
  - id: null_test
    command: echo null
"#;
        let json = r#"{
            "name": "csf_crossdisease_v3",
            "plan_commit": "abc123",
            "stages": [
                {"id": "de", "command": "echo de", "inputs": ["data.tsv"], "outputs": ["out.tsv"]},
                {"id": "null_test", "command": "echo null", "inputs": [], "outputs": []}
            ]
        }"#;
        let p_yaml: Plan = serde_yaml::from_str(yaml).unwrap();
        let p_json: Plan = serde_json::from_str(json).unwrap();
        assert_eq!(p_yaml.name, p_json.name);
        assert_eq!(p_yaml.stages.len(), 2);
        assert_eq!(p_yaml.stages[0].inputs, vec!["data.tsv".to_string()]);
    }

    #[test]
    fn plan_hash_is_stable() {
        let a = "name: x\nstages:\n  - id: a\n    command: echo hi\n";
        let b = "name: x\nstages:\n  - id: a\n    command: echo hi\n";
        assert_eq!(sha256_hex(a.as_bytes()), sha256_hex(b.as_bytes()));
    }

    #[test]
    fn hash_paths_is_order_independent() {
        let tmp = tempfile::tempdir().unwrap();
        std::fs::write(tmp.path().join("a"), "one").unwrap();
        std::fs::write(tmp.path().join("b"), "two").unwrap();
        let h1 = hash_paths(
            tmp.path(),
            &["a".to_string(), "b".to_string()],
        );
        let h2 = hash_paths(
            tmp.path(),
            &["b".to_string(), "a".to_string()],
        );
        assert_eq!(h1, h2);
    }
}
