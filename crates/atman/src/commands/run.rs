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
//!
//! Plans may declare `vars:` (name → string) substituted as `${name}` in
//! every stage's `command`, `inputs`, and `outputs`. `--dry-run` validates a
//! plan without executing it: every input must either exist on disk or be
//! declared as an output of an earlier stage. After each stage, a declared
//! output's sibling `<output>.run.json` sidecar (if present) is hashed into
//! the `sidecar_hash` column, and `--strict-outputs` (default on) marks a
//! stage failed when a declared output is missing afterwards.

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

    /// Validate the plan (schema, `${vars}`, input provenance) and list the
    /// resolved stages without executing anything or writing a manifest.
    #[arg(long, default_value_t = false)]
    dry_run: bool,

    /// Treat a stage whose declared outputs are missing after it runs as
    /// failed (manifest exit_code 2). Pass `--strict-outputs false` to keep
    /// the historical behaviour of recording `MISSING` and continuing.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set, num_args = 1)]
    strict_outputs: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Plan {
    name: String,
    #[serde(default)]
    plan_commit: Option<String>,
    /// `${name}` substitutions applied to every stage's command, inputs, and outputs.
    #[serde(default)]
    vars: std::collections::BTreeMap<String, String>,
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
    let plan = substitute_vars(plan)?;
    let plan_hash = sha256_hex(plan_text.as_bytes());
    if args.dry_run {
        return dry_run(&plan, &args.input_dir);
    }
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating output dir {:?}", args.output_dir))?;

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
        let row = execute_stage(
            stage,
            &args.input_dir,
            &plan,
            &plan_hash,
            args.strict_outputs,
        )?;
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
        bail!(
            "one or more stages failed; manifest written to {:?}",
            manifest_path
        );
    }
    Ok(())
}

/// Replace every `${name}` in stage commands, inputs, and outputs with the
/// plan's `vars`. An undefined name is an error; `$` not followed by `{` is
/// left alone (shell variables in commands still work).
fn substitute_vars(mut plan: Plan) -> Result<Plan> {
    let vars = plan.vars.clone();
    let apply = |text: &str, stage_id: &str| -> Result<String> {
        let mut out = String::with_capacity(text.len());
        let mut rest = text;
        while let Some(start) = rest.find("${") {
            out.push_str(&rest[..start]);
            let after = &rest[start + 2..];
            let Some(end) = after.find('}') else {
                bail!("stage {:?}: unterminated ${{ in {:?}", stage_id, text);
            };
            let name = &after[..end];
            match vars.get(name) {
                Some(value) => out.push_str(value),
                None => bail!(
                    "stage {:?}: undefined variable ${{{}}} (declare it under `vars:`)",
                    stage_id,
                    name
                ),
            }
            rest = &after[end + 1..];
        }
        out.push_str(rest);
        Ok(out)
    };
    for stage in &mut plan.stages {
        stage.command = apply(&stage.command, &stage.id)?;
        stage.inputs = stage
            .inputs
            .iter()
            .map(|p| apply(p, &stage.id))
            .collect::<Result<_>>()?;
        stage.outputs = stage
            .outputs
            .iter()
            .map(|p| apply(p, &stage.id))
            .collect::<Result<_>>()?;
    }
    Ok(plan)
}

/// Validate a plan without executing it. Every input must exist under
/// `cwd` or be declared as an output of an earlier stage; the resolved
/// stages are listed on stderr.
fn dry_run(plan: &Plan, cwd: &Path) -> Result<()> {
    let mut produced: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    let mut problems: Vec<String> = Vec::new();
    for (i, stage) in plan.stages.iter().enumerate() {
        eprintln!("stage {} [{}]: {}", i + 1, stage.id, stage.command);
        for input in &stage.inputs {
            let on_disk = cwd.join(input).exists();
            let from_earlier = produced.contains(input);
            eprintln!(
                "  input  {}{}",
                input,
                if from_earlier {
                    "  (produced by an earlier stage)"
                } else if on_disk {
                    ""
                } else {
                    "  (MISSING)"
                }
            );
            if !on_disk && !from_earlier {
                problems.push(format!(
                    "stage {:?}: input {:?} does not exist and no earlier stage declares it as an output",
                    stage.id, input
                ));
            }
        }
        for output in &stage.outputs {
            eprintln!("  output {}", output);
            produced.insert(output.clone());
        }
        if stage.outputs.is_empty() {
            eprintln!("  (no declared outputs)");
        }
    }
    if !problems.is_empty() {
        bail!(
            "dry run found {} problem(s):\n{}",
            problems.len(),
            problems.join("\n")
        );
    }
    eprintln!(
        "run: dry run ok — plan={} stages={} (nothing executed, no manifest written)",
        plan.name,
        plan.stages.len()
    );
    Ok(())
}

fn parse_plan(text: &str, path: &Path) -> Result<Plan> {
    let ext = path
        .extension()
        .and_then(|os| os.to_str())
        .map(|s| s.to_ascii_lowercase())
        .unwrap_or_default();
    let plan: Plan =
        match ext.as_str() {
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
    /// Hash over the `<output>.run.json` sidecars that exist next to the
    /// declared outputs after the stage ran (empty when there are none).
    sidecar_hash: String,
    /// Declared outputs that exist but carry no data row.
    ///
    /// Existence and non-emptiness are different audits, and a check
    /// that conflates them gives false assurance about exactly the
    /// stages worth inspecting: a header-only table is consumed
    /// downstream as a run of absences rather than as a missing input.
    n_empty_outputs: usize,
}

/// Whether a declared output exists but carries no data.
///
/// Two cases: zero bytes, or a `.tsv` holding nothing but its header.
///
/// The header rule is deliberately narrow — a single line containing a
/// tab. Every atman writer emits a tab-separated header and then rows,
/// so a one-line `.tsv` with a tab is a table with no rows; a one-line
/// file *without* a tab is more likely a legitimate single-value
/// output, and flagging it would fail a stage that did its job. Erring
/// toward not flagging is the right direction for a check that can fail
/// a pipeline.
fn output_is_empty(path: &Path) -> bool {
    let Ok(meta) = std::fs::metadata(path) else {
        // Missing is a different finding, reported separately.
        return false;
    };
    if meta.len() == 0 {
        return true;
    }
    if path.extension().and_then(|e| e.to_str()) != Some("tsv") {
        return false;
    }
    let Ok(text) = std::fs::read_to_string(path) else {
        return false;
    };
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    match (lines.next(), lines.next()) {
        (Some(header), None) => header.contains('\t'),
        _ => false,
    }
}

fn execute_stage(
    stage: &Stage,
    cwd: &Path,
    plan: &Plan,
    plan_hash: &str,
    strict_outputs: bool,
) -> Result<ManifestRow> {
    let input_hash = hash_paths(cwd, &stage.inputs)?;
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
    let mut exit_code = output.status.code().unwrap_or(-1);
    // Pipe captured stderr through so users see progress; stdout would be noisy.
    if !output.stderr.is_empty() {
        eprintln!("--- stage {} stderr ---", stage.id);
        eprintln!("{}", String::from_utf8_lossy(&output.stderr).trim_end());
    }
    let output_hash = hash_paths(cwd, &stage.outputs)?;
    let missing: Vec<&String> = stage
        .outputs
        .iter()
        .filter(|p| !cwd.join(p).exists())
        .collect();
    if strict_outputs && exit_code == 0 && !missing.is_empty() {
        eprintln!(
            "run: stage {} declared outputs missing after it ran: {}",
            stage.id,
            missing
                .iter()
                .map(|p| p.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        );
        exit_code = 2;
    }
    // Existence was checked above; content is a separate question.
    let empty: Vec<&String> = stage
        .outputs
        .iter()
        .filter(|p| {
            let full = cwd.join(p);
            full.exists() && output_is_empty(&full)
        })
        .collect();
    if !empty.is_empty() {
        eprintln!(
            "run: stage {} produced {} declared output(s) with no data rows: {}. A header-only \
             table is not a missing file and will be consumed downstream as a run of absences \
             rather than as an absent input.",
            stage.id,
            empty.len(),
            empty
                .iter()
                .map(|p| p.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        );
        if strict_outputs && exit_code == 0 {
            exit_code = 2;
        }
    }
    let n_empty_outputs = empty.len();

    let sidecars: Vec<String> = stage
        .outputs
        .iter()
        .map(|p| format!("{p}.run.json"))
        .filter(|p| cwd.join(p).exists())
        .collect();
    let sidecar_hash = hash_paths(cwd, &sidecars)?;
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
        sidecar_hash,
        n_empty_outputs,
    })
}

fn hash_paths(cwd: &Path, paths: &[String]) -> Result<String> {
    if paths.is_empty() {
        return Ok(String::new());
    }
    // Sort for determinism; hash concatenation of (relative path, file hash).
    let mut parts: Vec<(String, String)> = paths
        .iter()
        .map(|p| {
            let abs = cwd.join(p);
            let hash = hash_file_or_missing(&abs)?;
            Ok((p.clone(), hash))
        })
        .collect::<Result<_>>()?;
    parts.sort();
    let mut buf = Vec::new();
    for (path, file_hash) in parts {
        buf.extend_from_slice(path.as_bytes());
        buf.push(0);
        buf.extend_from_slice(file_hash.as_bytes());
        buf.push(b'\n');
    }
    Ok(sha256_hex(&buf))
}

/// Hash a declared input/output. A genuinely absent path hashes to the
/// sentinel `"MISSING"`; any other I/O error (permission denied, corrupted
/// read) is propagated rather than masquerading as an absent file, so a
/// manifest never silently records `"MISSING"` for a file that exists but
/// could not be read.
fn hash_file_or_missing(path: &Path) -> Result<String> {
    match std::fs::read(path) {
        Ok(bytes) => Ok(sha256_hex(&bytes)),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok("MISSING".to_string()),
        Err(e) => Err(e).with_context(|| format!("hashing declared path {:?}", path)),
    }
}

fn check_manifest_consistency(
    manifest_path: &Path,
    new_plan_hash: &str,
    plan: &Plan,
) -> Result<()> {
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
        "plan_name\tplan_commit\tplan_hash\tstage_id\tcommand\tinput_hash\toutput_hash\truntime_s\texit_code\tatman_version\tsystem\tstarted_at_unix_s\tsidecar_hash\tn_empty_outputs\n",
    );
    for row in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
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
            row.sidecar_hash,
            row.n_empty_outputs,
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
    fn vars_substitute_and_reject_undefined() {
        let plan: Plan = serde_yaml::from_str(
            "name: x\nvars:\n  out: results\n  seed: '7'\nstages:\n  - id: a\n    command: atman de --seed ${seed} --output-dir ${out}/de\n    inputs: [\"${out}/in.tsv\"]\n    outputs: [\"${out}/de/de_results.tsv\"]\n",
        )
        .unwrap();
        let plan = substitute_vars(plan).unwrap();
        assert_eq!(
            plan.stages[0].command,
            "atman de --seed 7 --output-dir results/de"
        );
        assert_eq!(plan.stages[0].inputs, vec!["results/in.tsv".to_string()]);
        assert_eq!(
            plan.stages[0].outputs,
            vec!["results/de/de_results.tsv".to_string()]
        );
        let bad: Plan =
            serde_yaml::from_str("name: x\nstages:\n  - id: a\n    command: echo ${nope}\n")
                .unwrap();
        let err = substitute_vars(bad).unwrap_err().to_string();
        assert!(err.contains("undefined variable ${nope}"), "{err}");
        // A bare `$` (shell variable) is left untouched.
        let shell: Plan =
            serde_yaml::from_str("name: x\nstages:\n  - id: a\n    command: echo $HOME\n").unwrap();
        assert_eq!(
            substitute_vars(shell).unwrap().stages[0].command,
            "echo $HOME"
        );
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
        let h1 = hash_paths(tmp.path(), &["a".to_string(), "b".to_string()]).unwrap();
        let h2 = hash_paths(tmp.path(), &["b".to_string(), "a".to_string()]).unwrap();
        assert_eq!(h1, h2);
    }
}
