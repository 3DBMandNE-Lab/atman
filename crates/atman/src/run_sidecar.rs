//! Invocation-provenance sidecar shared by every atman subcommand.
//!
//! Writes a JSON file next to the command's primary output capturing enough
//! detail to reproduce the run from the sidecar alone: atman version + git
//! SHA, build environment (rustc version, Cargo.lock hash, profile, target
//! triple), fully resolved argument dict, per-file SHA-256 of every input
//! the command read, per-file SHA-256 of every output it wrote, UTC
//! start/finish timestamps, per-run UUID, cwd anchor, and a reconstructed
//! `reinvoke` command-line string.
//!
//! ## Schema version 1
//!
//! The current on-disk format is `schema_version: 1`. Readers that pre-date
//! the bump see no `schema_version` key and should treat the file as v0:
//! in v0, inputs are a single roll-up `input_dir_sha256: string`; in v1,
//! inputs are a per-file dict `inputs_sha256: {relative_path: sha256}`.
//! v1 is otherwise a superset of v0.
//!
//! Re-exported from `crate::io` so callers can write
//! `use crate::io::write_run_sidecar;` in the same spirit as the other TSV
//! helpers.

use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

use anyhow::{Context, Result};
use serde_json::{json, Map, Value};

use crate::io::{atomic_write, sha256_hex};

/// On-disk sidecar schema version. Bumped whenever an incompatible
/// change lands (field rename, semantic shift); additive changes
/// alone do not bump this.
pub const SIDECAR_SCHEMA_VERSION: u32 = 1;

/// Returns the sidecar path for a primary output file: `<primary>.run.json`.
pub fn sidecar_path_for(primary_output: &Path) -> PathBuf {
    let mut s = primary_output.as_os_str().to_os_string();
    s.push(".run.json");
    PathBuf::from(s)
}

/// Per-file input hashes keyed by a stable label (basename for canonical
/// TSV inputs, caller-supplied label for `hash_labeled_inputs`). Ordered
/// so the JSON output is deterministic under a fixed input set.
pub type InputHashes = BTreeMap<String, String>;

/// Per-file SHA-256 for the canonical input TSVs in a directory.
///
/// Keys are the basenames passed in `filenames`; values are hex SHA-256
/// of each file's bytes. Missing files are skipped so optional inputs
/// (e.g. `qc_measurements.tsv`) don't appear in the output at all when
/// absent, rather than silently altering a roll-up hash.
pub fn hash_canonical_inputs(input_dir: &Path, filenames: &[&str]) -> Result<InputHashes> {
    let entries: Vec<(&str, PathBuf)> = filenames.iter().map(|n| (*n, input_dir.join(n))).collect();
    let refs: Vec<(&str, &Path)> = entries.iter().map(|(n, p)| (*n, p.as_path())).collect();
    hash_labeled_inputs(&refs)
}

/// Per-file SHA-256 for arbitrary labeled inputs. Use when inputs are
/// not a single canonical directory (e.g. `align programs` takes
/// per-cohort loadings TSVs with caller-chosen labels).
pub fn hash_labeled_inputs(entries: &[(&str, &Path)]) -> Result<InputHashes> {
    let mut out: InputHashes = BTreeMap::new();
    for (label, path) in entries {
        if !path.exists() {
            continue;
        }
        let bytes = std::fs::read(path).with_context(|| format!("reading {:?}", path))?;
        out.insert((*label).to_string(), sha256_hex(&bytes));
    }
    Ok(out)
}

/// Write the one-per-invocation provenance sidecar.
///
/// The caller owns the `args` dict (typically built with
/// `serde_json::json!` using kebab-case keys that match CLI flag names).
///
/// `inputs_sha256` is a per-file hash dict from [`hash_canonical_inputs`]
/// or [`hash_labeled_inputs`]; `output_files` are hashed here and
/// emitted under `output_files` as `{path: "sha256:<hex>"}`.
///
/// `extras` optionally merges additional fields into the top-level
/// object — intended for recording *resolved values* of rule-form
/// arguments (e.g. `k_resolved: 7` when `--k cumulative-variance=0.8`
/// was requested). Keys here must not collide with the standard
/// fields; `schema_version`, `atman_version`, `build_env`,
/// `output_files`, etc. are enforced reserved and will be shadowed by
/// the canonical fields.
#[allow(clippy::too_many_arguments)]
pub fn write_run_sidecar(
    sidecar_path: &Path,
    command: &str,
    args: Value,
    inputs_sha256: &InputHashes,
    output_files: &[PathBuf],
    started_at: SystemTime,
    finished_at: SystemTime,
    extras: Option<Map<String, Value>>,
) -> Result<()> {
    let mut output_map = serde_json::Map::new();
    for p in output_files {
        let bytes = std::fs::read(p).with_context(|| format!("reading output {:?}", p))?;
        let hash = sha256_hex(&bytes);
        output_map.insert(
            p.to_string_lossy().into_owned(),
            Value::String(format!("sha256:{hash}")),
        );
    }

    let mut inputs_obj = Map::new();
    for (k, v) in inputs_sha256 {
        inputs_obj.insert(k.clone(), Value::String(format!("sha256:{v}")));
    }

    let cwd = std::env::current_dir()
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|_| String::new());

    let build_env = json!({
        "rustc_version":       env!("ATMAN_RUSTC_VERSION"),
        "cargo_lock_sha256":   env!("ATMAN_CARGO_LOCK_SHA256"),
        "profile":             env!("ATMAN_PROFILE"),
        "target_triple":       env!("ATMAN_TARGET"),
    });

    let reinvoke = reinvoke_string(command, &args);
    let run_uuid = uuid::Uuid::new_v4().to_string();

    let mut payload = Map::new();
    payload.insert("schema_version".into(), json!(SIDECAR_SCHEMA_VERSION));
    payload.insert("run_uuid".into(), json!(run_uuid));
    payload.insert("command".into(), json!(command));
    payload.insert("reinvoke".into(), json!(reinvoke));
    payload.insert("args".into(), args);
    payload.insert("atman_version".into(), json!(env!("CARGO_PKG_VERSION")));
    payload.insert("atman_git_sha".into(), json!(env!("ATMAN_GIT_SHA")));
    payload.insert("build_env".into(), build_env);
    payload.insert("cwd_at_start".into(), json!(cwd));
    payload.insert("inputs_sha256".into(), Value::Object(inputs_obj));
    payload.insert("output_files".into(), Value::Object(output_map));
    payload.insert("started_at".into(), json!(format_iso8601_utc(started_at)));
    payload.insert("finished_at".into(), json!(format_iso8601_utc(finished_at)));
    payload.insert("exit_code".into(), json!(0));

    // `os_arch` is retained for compatibility with v0 consumers that
    // read it directly. Identical to `build_env.target_triple`.
    payload.insert("os_arch".into(), json!(env!("ATMAN_TARGET")));

    if let Some(extra) = extras {
        // Reserved keys cannot be overwritten — protects the schema
        // against commands that might accidentally emit a field name
        // that collides with a canonical one.
        const RESERVED: &[&str] = &[
            "schema_version",
            "run_uuid",
            "command",
            "reinvoke",
            "args",
            "atman_version",
            "atman_git_sha",
            "build_env",
            "cwd_at_start",
            "inputs_sha256",
            "output_files",
            "started_at",
            "finished_at",
            "exit_code",
            "os_arch",
        ];
        for (k, v) in extra {
            if RESERVED.contains(&k.as_str()) {
                continue;
            }
            payload.insert(k, v);
        }
    }

    let mut text = serde_json::to_string_pretty(&Value::Object(payload))?;
    text.push('\n');
    atomic_write(sidecar_path, text.as_bytes())
}

/// Best-effort reconstruction of the command line.
///
/// Produces `atman <command> [--key value]*` by walking the resolved
/// argument dict in sort order. Handles the common types:
///
/// - bool: emits `--key true` / `--key false` (clap accepts either).
/// - string / number: emits `--key <value>`, shell-quoting if the
///   value contains whitespace, single quotes, or `$`.
/// - null: skipped (argument was unset).
/// - array / object: serialized as a single JSON token; the reviewer
///   may have to hand-edit these for complex nested args.
///
/// The raw `args` dict is also in the sidecar, so `reinvoke` is a
/// convenience for paste-and-run, not a source of truth.
fn reinvoke_string(command: &str, args: &Value) -> String {
    let mut parts: Vec<String> = vec!["atman".into()];
    for tok in command.split_whitespace() {
        parts.push(tok.to_string());
    }
    if let Value::Object(map) = args {
        // Deterministic iteration order: BTreeMap-sorted keys.
        let mut keys: Vec<&String> = map.keys().collect();
        keys.sort();
        for k in keys {
            let v = &map[k];
            let flag = format!("--{k}");
            match v {
                Value::Null => continue,
                Value::Bool(b) => {
                    parts.push(flag);
                    parts.push(b.to_string());
                }
                Value::Number(n) => {
                    parts.push(flag);
                    parts.push(n.to_string());
                }
                Value::String(s) => {
                    parts.push(flag);
                    parts.push(shell_quote(s));
                }
                other => {
                    parts.push(flag);
                    parts.push(shell_quote(&other.to_string()));
                }
            }
        }
    }
    parts.join(" ")
}

/// POSIX single-quote escaping sufficient for pasting into sh / bash / zsh.
fn shell_quote(s: &str) -> String {
    if !s.is_empty()
        && s.chars()
            .all(|c| c.is_ascii_alphanumeric() || "-_.,/=@:+".contains(c))
    {
        return s.to_string();
    }
    let mut out = String::with_capacity(s.len() + 2);
    out.push('\'');
    for c in s.chars() {
        if c == '\'' {
            out.push_str(r"'\''");
        } else {
            out.push(c);
        }
    }
    out.push('\'');
    out
}

/// Format a `SystemTime` as ISO-8601 UTC with second precision:
/// `YYYY-MM-DDTHH:MM:SSZ`. Uses Howard Hinnant's civil-from-days algorithm to
/// avoid pulling a time-formatting dependency into the workspace.
pub fn format_iso8601_utc(t: SystemTime) -> String {
    let secs = match t.duration_since(UNIX_EPOCH) {
        Ok(d) => d.as_secs() as i64,
        Err(e) => -(e.duration().as_secs() as i64),
    };
    let days = secs.div_euclid(86_400);
    let sod = secs.rem_euclid(86_400);
    let (year, month, day) = civil_from_days(days);
    let hour = sod / 3600;
    let minute = (sod % 3600) / 60;
    let second = sod % 60;
    format!("{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}Z")
}

fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 }.div_euclid(146_097);
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    (y, m as u32, d as u32)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;
    use tempfile::TempDir;

    #[test]
    fn sidecar_path_appends_run_json_suffix() {
        let p = Path::new("/out/loadings.tsv");
        assert_eq!(
            sidecar_path_for(p),
            PathBuf::from("/out/loadings.tsv.run.json")
        );
    }

    #[test]
    fn format_iso8601_utc_epoch_is_zero() {
        assert_eq!(format_iso8601_utc(UNIX_EPOCH), "1970-01-01T00:00:00Z");
    }

    #[test]
    fn format_iso8601_utc_handles_leap_day() {
        let t = UNIX_EPOCH + Duration::from_secs(1_709_210_096);
        assert_eq!(format_iso8601_utc(t), "2024-02-29T12:34:56Z");
    }

    #[test]
    fn format_iso8601_utc_handles_manuscript_timestamp() {
        let t = UNIX_EPOCH + Duration::from_secs(1_776_572_207);
        assert_eq!(format_iso8601_utc(t), "2026-04-19T04:16:47Z");
    }

    #[test]
    fn hash_canonical_inputs_skips_missing_files() {
        let dir = TempDir::new().unwrap();
        std::fs::write(dir.path().join("a.tsv"), b"alpha").unwrap();
        let out = hash_canonical_inputs(dir.path(), &["a.tsv", "missing.tsv"]).unwrap();
        assert_eq!(out.len(), 1);
        assert!(out.contains_key("a.tsv"));
        assert!(!out.contains_key("missing.tsv"));
    }

    #[test]
    fn hash_canonical_inputs_is_per_file_and_order_independent() {
        let dir = TempDir::new().unwrap();
        std::fs::write(dir.path().join("a.tsv"), b"alpha").unwrap();
        std::fs::write(dir.path().join("b.tsv"), b"beta").unwrap();
        let h1 = hash_canonical_inputs(dir.path(), &["a.tsv", "b.tsv"]).unwrap();
        let h2 = hash_canonical_inputs(dir.path(), &["b.tsv", "a.tsv"]).unwrap();
        assert_eq!(h1, h2);
        // Per-file, not a single roll-up: changing one file changes
        // only that file's value, not the whole map.
        std::fs::write(dir.path().join("a.tsv"), b"ALPHA").unwrap();
        let h3 = hash_canonical_inputs(dir.path(), &["a.tsv", "b.tsv"]).unwrap();
        assert_ne!(h1["a.tsv"], h3["a.tsv"]);
        assert_eq!(h1["b.tsv"], h3["b.tsv"]);
    }

    #[test]
    fn write_run_sidecar_v1_has_all_required_fields() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("loadings.tsv");
        std::fs::write(&out, b"program\tassay\n").unwrap();
        let sidecar = sidecar_path_for(&out);
        let t0 = UNIX_EPOCH + Duration::from_secs(1_700_000_000);
        let t1 = UNIX_EPOCH + Duration::from_secs(1_700_000_042);
        let mut inputs = InputHashes::new();
        inputs.insert("samples.tsv".into(), "deadbeef".into());
        inputs.insert("proteins.tsv".into(), "cafef00d".into());
        write_run_sidecar(
            &sidecar,
            "decompose ica",
            json!({"k": 2, "seed": 20260418, "max-iter": 200, "robust": true, "tol": null}),
            &inputs,
            std::slice::from_ref(&out),
            t0,
            t1,
            None,
        )
        .unwrap();
        let text = std::fs::read_to_string(&sidecar).unwrap();
        let parsed: Value = serde_json::from_str(&text).unwrap();

        // v1 contract fields
        assert_eq!(parsed["schema_version"], 1);
        assert!(parsed["run_uuid"].as_str().unwrap().len() == 36);
        assert_eq!(parsed["command"], "decompose ica");
        let reinvoke = parsed["reinvoke"].as_str().unwrap();
        assert!(reinvoke.starts_with("atman decompose ica"));
        assert!(reinvoke.contains("--k 2"));
        assert!(reinvoke.contains("--seed 20260418"));
        assert!(reinvoke.contains("--robust true"));
        // `tol` was null — should be skipped in reinvoke
        assert!(!reinvoke.contains("--tol"));

        // Per-file inputs; v0 roll-up is gone
        assert!(parsed.get("input_dir_sha256").is_none());
        assert_eq!(parsed["inputs_sha256"]["samples.tsv"], "sha256:deadbeef");
        assert_eq!(parsed["inputs_sha256"]["proteins.tsv"], "sha256:cafef00d");

        // Build env
        assert!(parsed["build_env"]["rustc_version"].is_string());
        assert!(parsed["build_env"]["cargo_lock_sha256"].is_string());
        assert!(parsed["build_env"]["profile"].is_string());
        assert!(parsed["build_env"]["target_triple"].is_string());

        // Existing v0 fields kept
        assert!(parsed["atman_version"].is_string());
        assert!(parsed["atman_git_sha"].is_string());
        assert!(parsed["os_arch"].is_string());
        assert_eq!(parsed["started_at"], "2023-11-14T22:13:20Z");
        assert_eq!(parsed["finished_at"], "2023-11-14T22:14:02Z");
        assert_eq!(parsed["exit_code"], 0);
        assert!(parsed["cwd_at_start"].is_string());
        let out_key = out.to_string_lossy().into_owned();
        assert!(parsed["output_files"][&out_key]
            .as_str()
            .unwrap()
            .starts_with("sha256:"));
    }

    #[test]
    fn write_run_sidecar_merges_extras_and_rejects_reserved_keys() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("result.tsv");
        std::fs::write(&out, b"x\n").unwrap();
        let sidecar = sidecar_path_for(&out);
        let t0 = UNIX_EPOCH + Duration::from_secs(1_700_000_000);
        let t1 = t0 + Duration::from_secs(1);
        let mut extras = Map::new();
        extras.insert("k_resolved".into(), json!(7));
        extras.insert("selection_rule".into(), json!("cumulative-variance=0.80"));
        // Attempt to overwrite a reserved field — should be ignored.
        extras.insert("schema_version".into(), json!(999));
        extras.insert("command".into(), json!("evil-command"));
        write_run_sidecar(
            &sidecar,
            "decompose ica",
            json!({"k": "auto"}),
            &InputHashes::new(),
            std::slice::from_ref(&out),
            t0,
            t1,
            Some(extras),
        )
        .unwrap();
        let parsed: Value =
            serde_json::from_str(&std::fs::read_to_string(&sidecar).unwrap()).unwrap();
        assert_eq!(parsed["k_resolved"], 7);
        assert_eq!(parsed["selection_rule"], "cumulative-variance=0.80");
        // Reserved keys survived the merge attempt
        assert_eq!(parsed["schema_version"], 1);
        assert_eq!(parsed["command"], "decompose ica");
    }

    #[test]
    fn reinvoke_shell_quotes_values_with_spaces() {
        let s = reinvoke_string(
            "de",
            &json!({
                "input-dir": "/tmp/my data",
                "groups": "A-B,C-D",
                "min-pairs": 5,
            }),
        );
        assert!(s.contains("--input-dir '/tmp/my data'"));
        // Values with commas/dashes that don't need quoting stay bare
        assert!(s.contains("--groups A-B,C-D"));
        assert!(s.contains("--min-pairs 5"));
    }
}
