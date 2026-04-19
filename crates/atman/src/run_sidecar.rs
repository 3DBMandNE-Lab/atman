//! Invocation-provenance sidecar shared by every atman subcommand.
//!
//! Writes a JSON file next to the command's primary output capturing enough
//! detail to reproduce the run from the sidecar alone: atman version + git SHA,
//! the full resolved argument dict, a hash of the canonical input TSVs, hashes
//! of each output file, UTC start/finish timestamps, and the OS/arch triple.
//!
//! Re-exported from `crate::io` so callers can write
//! `use crate::io::write_run_sidecar;` in the same spirit as the other TSV
//! helpers.

use std::{
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

use anyhow::{Context, Result};
use serde_json::{json, Value};

use crate::io::{atomic_write, sha256_hex};

/// Returns the sidecar path for a primary output file: `<primary>.run.json`.
pub fn sidecar_path_for(primary_output: &Path) -> PathBuf {
    let mut s = primary_output.as_os_str().to_os_string();
    s.push(".run.json");
    PathBuf::from(s)
}

/// SHA-256 of a sorted `<basename>\t<sha256(contents)>\n` concatenation over
/// the canonical input files the command actually reads.
///
/// Missing files are skipped (optional canonical TSVs like `qc_measurements.tsv`
/// shouldn't flip the hash when absent). Basenames — not full paths — are
/// hashed so the same canonical directory hashes identically regardless of
/// where on disk it lives.
pub fn hash_canonical_inputs(input_dir: &Path, filenames: &[&str]) -> Result<String> {
    let mut entries: Vec<(&str, PathBuf)> = Vec::with_capacity(filenames.len());
    for name in filenames {
        entries.push((*name, input_dir.join(name)));
    }
    let refs: Vec<(&str, &Path)> = entries.iter().map(|(n, p)| (*n, p.as_path())).collect();
    hash_labeled_inputs(&refs)
}

/// SHA-256 of a sorted `<label>\t<sha256(contents)>\n` concatenation over
/// arbitrary input files. Use when inputs are not a single canonical
/// directory (e.g. `align programs` takes per-cohort loadings TSVs).
/// Missing files are skipped.
pub fn hash_labeled_inputs(entries: &[(&str, &Path)]) -> Result<String> {
    let mut pairs: Vec<(String, String)> = Vec::new();
    for (label, path) in entries {
        if !path.exists() {
            continue;
        }
        let bytes = std::fs::read(path).with_context(|| format!("reading {:?}", path))?;
        pairs.push(((*label).to_string(), sha256_hex(&bytes)));
    }
    pairs.sort();
    let mut concat = String::new();
    for (label, hash) in &pairs {
        concat.push_str(label);
        concat.push('\t');
        concat.push_str(hash);
        concat.push('\n');
    }
    Ok(sha256_hex(concat.as_bytes()))
}

/// Write the one-per-invocation provenance sidecar.
///
/// The caller owns the full argument dict (typically built with `serde_json::json!`
/// and kebab-case keys that match the CLI flag names, so the sidecar reads the
/// way a reviewer would invoke the command on the shell).
///
/// `output_files` are hashed in place — list every artifact this invocation
/// wrote except the sidecar itself.
pub fn write_run_sidecar(
    sidecar_path: &Path,
    command: &str,
    args: Value,
    input_dir_sha256: &str,
    output_files: &[PathBuf],
    started_at: SystemTime,
    finished_at: SystemTime,
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
    let payload = json!({
        "atman_version": env!("CARGO_PKG_VERSION"),
        "atman_git_sha": env!("ATMAN_GIT_SHA"),
        "command": command,
        "args": args,
        "input_dir_sha256": input_dir_sha256,
        "output_files": Value::Object(output_map),
        "started_at": format_iso8601_utc(started_at),
        "finished_at": format_iso8601_utc(finished_at),
        "os_arch": env!("ATMAN_TARGET"),
        "exit_code": 0,
    });
    let mut text = serde_json::to_string_pretty(&payload)?;
    text.push('\n');
    atomic_write(sidecar_path, text.as_bytes())
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
        // 2024-02-29T12:34:56Z → seconds since epoch:
        // days = 54*365 + 13 leaps + 31 + 28 = 19782; 19782*86400 + 45296 = 1_709_210_096
        let t = UNIX_EPOCH + Duration::from_secs(1_709_210_096);
        assert_eq!(format_iso8601_utc(t), "2024-02-29T12:34:56Z");
    }

    #[test]
    fn format_iso8601_utc_handles_manuscript_timestamp() {
        // 2026-04-19T04:16:47Z matches the reviewer-feedback example.
        // Days to 2026-01-01 = 56*365 + 14 leaps = 20454; +31+28+31+18 = 108 →
        // day 20562. 20562*86400 = 1_776_556_800. +4*3600 + 16*60 + 47 = 15_407.
        let t = UNIX_EPOCH + Duration::from_secs(1_776_572_207);
        assert_eq!(format_iso8601_utc(t), "2026-04-19T04:16:47Z");
    }

    #[test]
    fn hash_canonical_inputs_skips_missing_files() {
        let dir = TempDir::new().unwrap();
        std::fs::write(dir.path().join("a.tsv"), b"alpha").unwrap();
        let with_missing =
            hash_canonical_inputs(dir.path(), &["a.tsv", "missing.tsv"]).unwrap();
        let only_present = hash_canonical_inputs(dir.path(), &["a.tsv"]).unwrap();
        assert_eq!(with_missing, only_present);
    }

    #[test]
    fn hash_canonical_inputs_is_order_independent() {
        let dir = TempDir::new().unwrap();
        std::fs::write(dir.path().join("a.tsv"), b"alpha").unwrap();
        std::fs::write(dir.path().join("b.tsv"), b"beta").unwrap();
        let h1 = hash_canonical_inputs(dir.path(), &["a.tsv", "b.tsv"]).unwrap();
        let h2 = hash_canonical_inputs(dir.path(), &["b.tsv", "a.tsv"]).unwrap();
        assert_eq!(h1, h2);
    }

    #[test]
    fn hash_canonical_inputs_changes_with_content() {
        let dir = TempDir::new().unwrap();
        std::fs::write(dir.path().join("a.tsv"), b"alpha").unwrap();
        let h1 = hash_canonical_inputs(dir.path(), &["a.tsv"]).unwrap();
        std::fs::write(dir.path().join("a.tsv"), b"ALPHA").unwrap();
        let h2 = hash_canonical_inputs(dir.path(), &["a.tsv"]).unwrap();
        assert_ne!(h1, h2);
    }

    #[test]
    fn write_run_sidecar_round_trip_parses_as_json_with_expected_keys() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("loadings.tsv");
        std::fs::write(&out, b"program\tassay\n").unwrap();
        let sidecar = sidecar_path_for(&out);
        // 2023-11-14T22:13:20Z and +42s
        let t0 = UNIX_EPOCH + Duration::from_secs(1_700_000_000);
        let t1 = UNIX_EPOCH + Duration::from_secs(1_700_000_042);
        write_run_sidecar(
            &sidecar,
            "decompose ica",
            json!({"k": 2, "seed": 20260418}),
            "deadbeef",
            &[out.clone()],
            t0,
            t1,
        )
        .unwrap();
        let text = std::fs::read_to_string(&sidecar).unwrap();
        let parsed: Value = serde_json::from_str(&text).unwrap();
        assert_eq!(parsed["command"], "decompose ica");
        assert_eq!(parsed["args"]["k"], 2);
        assert_eq!(parsed["args"]["seed"], 20260418);
        assert_eq!(parsed["input_dir_sha256"], "deadbeef");
        assert_eq!(parsed["started_at"], "2023-11-14T22:13:20Z");
        assert_eq!(parsed["finished_at"], "2023-11-14T22:14:02Z");
        assert_eq!(parsed["exit_code"], 0);
        assert!(parsed["atman_version"].is_string());
        assert!(parsed["atman_git_sha"].is_string());
        assert!(parsed["os_arch"].is_string());
        let out_key = out.to_string_lossy().into_owned();
        assert!(parsed["output_files"][&out_key]
            .as_str()
            .unwrap()
            .starts_with("sha256:"));
    }
}
