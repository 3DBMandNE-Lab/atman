//! Capture the git commit SHA and target triple at build time so every
//! subcommand's run-sidecar can report the exact binary that wrote the output.
//!
//! Sets:
//! - `ATMAN_GIT_SHA` — output of `git rev-parse HEAD`, or `"unknown"` when the
//!   binary is built outside a git checkout (e.g. from a crates.io tarball).
//! - `ATMAN_TARGET`  — the Cargo `TARGET` triple (e.g. `aarch64-apple-darwin`).

use std::{path::Path, process::Command};

fn main() {
    let sha = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()
        .and_then(|o| {
            if o.status.success() {
                String::from_utf8(o.stdout).ok()
            } else {
                None
            }
        })
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "unknown".to_string());
    println!("cargo:rustc-env=ATMAN_GIT_SHA={sha}");

    let target = std::env::var("TARGET").unwrap_or_default();
    println!("cargo:rustc-env=ATMAN_TARGET={target}");

    // Rebuild when HEAD moves so the baked-in SHA stays current across commits.
    println!("cargo:rerun-if-changed=build.rs");
    let head = Path::new("../../.git/HEAD");
    if head.exists() {
        println!("cargo:rerun-if-changed=../../.git/HEAD");
        if let Ok(text) = std::fs::read_to_string(head) {
            if let Some(rest) = text.trim().strip_prefix("ref: ") {
                println!("cargo:rerun-if-changed=../../.git/{rest}");
            }
        }
    }
}
