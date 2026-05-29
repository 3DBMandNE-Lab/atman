//! Capture build-time provenance into compile-time env vars so every
//! subcommand's run sidecar can report the exact binary that wrote the
//! output.
//!
//! Sets:
//! - `ATMAN_GIT_SHA` — `git rev-parse HEAD`, or `"unknown"` when built
//!   outside a git checkout. A `-dirty` suffix is appended when the
//!   working tree has uncommitted tracked changes at build time, so a
//!   binary built from a dirty tree never claims a clean commit.
//! - `ATMAN_TARGET` — Cargo `TARGET` triple (e.g. `aarch64-apple-darwin`).
//! - `ATMAN_RUSTC_VERSION` — `rustc --version` verbatim, or `"unknown"`.
//! - `ATMAN_CARGO_LOCK_SHA256` — SHA-256 of the workspace `Cargo.lock`,
//!   or `"unknown"` when the lock file is absent (shouldn't happen in a
//!   normal build, but `cargo install` from a git dep can).
//! - `ATMAN_PROFILE` — `debug` / `release`, from Cargo.

use std::{
    path::{Path, PathBuf},
    process::Command,
};

fn main() {
    let sha = git_sha_with_dirty();
    println!("cargo:rustc-env=ATMAN_GIT_SHA={sha}");

    let target = std::env::var("TARGET").unwrap_or_default();
    println!("cargo:rustc-env=ATMAN_TARGET={target}");

    let rustc_version = std::env::var("RUSTC")
        .ok()
        .and_then(|rustc| {
            Command::new(rustc)
                .arg("--version")
                .output()
                .ok()
                .filter(|o| o.status.success())
                .and_then(|o| String::from_utf8(o.stdout).ok())
        })
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "unknown".to_string());
    println!("cargo:rustc-env=ATMAN_RUSTC_VERSION={rustc_version}");

    let lock_sha = ["Cargo.lock", "../../Cargo.lock", "../Cargo.lock"]
        .iter()
        .map(Path::new)
        .find(|p| p.exists())
        .and_then(|p| std::fs::read(p).ok())
        .map(|bytes| sha256_hex(&bytes))
        .unwrap_or_else(|| "unknown".to_string());
    println!("cargo:rustc-env=ATMAN_CARGO_LOCK_SHA256={lock_sha}");

    let profile = std::env::var("PROFILE").unwrap_or_else(|_| "unknown".to_string());
    println!("cargo:rustc-env=ATMAN_PROFILE={profile}");

    // Rebuild when HEAD moves so the baked-in SHA stays current across
    // commits. Also re-run if Cargo.lock changes (so ATMAN_CARGO_LOCK_SHA256
    // tracks dependency pin changes), and when this crate's own sources change
    // so the `-dirty` flag refreshes on uncommitted edits during incremental
    // dev builds (emitting explicit rerun-if-changed opts out of cargo's
    // default "any package file" trigger, so `src` must be listed explicitly).
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=src");
    println!("cargo:rerun-if-changed=../../Cargo.lock");
    // Resolve the git directory. In a normal checkout `../../.git` is a
    // directory; in a LINKED git worktree it is a FILE containing
    // `gitdir: <path-to-per-worktree-gitdir>`. The kanban build system runs
    // entirely in linked worktrees, so we must handle the file form — otherwise
    // the HEAD rerun trigger is silently skipped and the baked SHA goes stale
    // across commits in a worktree.
    if let Some(git_dir) = resolve_git_dir(Path::new("../../.git")) {
        let head = git_dir.join("HEAD");
        if head.exists() {
            println!("cargo:rerun-if-changed={}", head.display());
            if let Ok(text) = std::fs::read_to_string(&head) {
                if let Some(rest) = text.trim().strip_prefix("ref: ") {
                    // Branch refs are shared via the common dir (worktrees
                    // record it in a `commondir` file); resolve against it so
                    // the watched ref path is correct in both layouts.
                    let base = read_commondir(&git_dir).unwrap_or_else(|| git_dir.clone());
                    println!("cargo:rerun-if-changed={}", base.join(rest).display());
                }
            }
        }
    }
}

/// Resolve the real git directory for `marker` (`<repo>/.git`). Returns the
/// directory itself for a normal checkout, or the `gitdir:` target when
/// `marker` is the file form used by linked worktrees. `None` if neither.
fn resolve_git_dir(marker: &Path) -> Option<PathBuf> {
    if marker.is_dir() {
        return Some(marker.to_path_buf());
    }
    if marker.is_file() {
        let text = std::fs::read_to_string(marker).ok()?;
        let dir = text
            .lines()
            .find_map(|l| l.trim().strip_prefix("gitdir: ").map(str::trim))?;
        return Some(PathBuf::from(dir));
    }
    None
}

/// In a linked worktree the shared refs live in the common git dir, recorded
/// in `<git_dir>/commondir` as a path relative to `git_dir`. Returns that
/// resolved common dir, or `None` for a normal checkout (no `commondir` file).
fn read_commondir(git_dir: &Path) -> Option<PathBuf> {
    let rel = std::fs::read_to_string(git_dir.join("commondir")).ok()?;
    Some(git_dir.join(rel.trim()))
}

/// Resolve the build-time git SHA, appending a `-dirty` suffix when the
/// working tree has uncommitted tracked changes.
///
/// `git status --porcelain` emits one line per changed-or-untracked path and
/// nothing at all for a clean tree, so a non-empty (successful) output marks
/// the tree dirty. Untracked (non-ignored) files are INCLUDED: a new,
/// uncommitted source file means the build is not reproducible from the
/// recorded commit, which is exactly the provenance gap `-dirty` exists to
/// flag. (`.gitignore`d paths like `target/` are not reported, so routine
/// build artifacts do not spuriously trip it.)
///
/// When `git` is unavailable the base SHA is `"unknown"` (preserving the
/// existing fallback) and no suffix is added — we never claim dirtiness we
/// cannot observe. Side-effect-free and deterministic.
///
/// LIMITATION: the flag is captured when this build script runs. Cargo re-runs
/// it on changes to `build.rs`, this crate's `src`, `Cargo.lock`, and
/// `.git/HEAD`/the branch ref — so it refreshes on commits and on edits to
/// this crate. An uncommitted edit to a DIFFERENT workspace crate that does
/// not also rebuild this crate can leave a stale clean SHA. For guaranteed
/// provenance, build from a clean checkout (as CI/release does); for
/// incremental dev builds `-dirty` is best-effort.
fn git_sha_with_dirty() -> String {
    let sha = capture_command("git", &["rev-parse", "HEAD"]);
    if sha == "unknown" {
        return sha;
    }
    let status = Command::new("git")
        .args(["status", "--porcelain"])
        .output();
    let dirty = matches!(
        status,
        Ok(o) if o.status.success() && !o.stdout.iter().all(u8::is_ascii_whitespace)
    );
    if dirty {
        format!("{sha}-dirty")
    } else {
        sha
    }
}

fn capture_command(program: &str, args: &[&str]) -> String {
    Command::new(program)
        .args(args)
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "unknown".to_string())
}

/// Tiny SHA-256 without pulling a heavy dep into build-dependencies.
/// Correctness matters more than speed — this runs once per rebuild.
fn sha256_hex(bytes: &[u8]) -> String {
    let h = Sha256Builder::new().update(bytes).finalize();
    let mut out = String::with_capacity(64);
    for b in h.iter() {
        out.push_str(&format!("{b:02x}"));
    }
    out
}

/// Minimal SHA-256 implementation for build.rs usage only. Produces
/// identical output to the `sha2` crate on any RFC-6234 test vector.
struct Sha256Builder {
    state: [u32; 8],
    buf: Vec<u8>,
    len: u64,
}

const K: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

impl Sha256Builder {
    fn new() -> Self {
        Self {
            state: [
                0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
                0x5be0cd19,
            ],
            buf: Vec::with_capacity(64),
            len: 0,
        }
    }

    fn update(mut self, bytes: &[u8]) -> Self {
        self.len = self.len.wrapping_add(bytes.len() as u64);
        let mut rem = bytes;
        while !rem.is_empty() {
            let take = 64 - self.buf.len();
            let n = take.min(rem.len());
            self.buf.extend_from_slice(&rem[..n]);
            rem = &rem[n..];
            if self.buf.len() == 64 {
                let block: [u8; 64] = self.buf.as_slice().try_into().unwrap();
                self.state = compress(self.state, &block);
                self.buf.clear();
            }
        }
        self
    }

    fn finalize(mut self) -> [u8; 32] {
        let bit_len = self.len.wrapping_mul(8);
        self.buf.push(0x80);
        while self.buf.len() % 64 != 56 {
            self.buf.push(0);
        }
        self.buf.extend_from_slice(&bit_len.to_be_bytes());
        let chunks = self.buf.len() / 64;
        for c in 0..chunks {
            let block: [u8; 64] = self.buf[c * 64..(c + 1) * 64].try_into().unwrap();
            self.state = compress(self.state, &block);
        }
        let mut out = [0u8; 32];
        for (i, w) in self.state.iter().enumerate() {
            out[i * 4..(i + 1) * 4].copy_from_slice(&w.to_be_bytes());
        }
        out
    }
}

fn compress(state: [u32; 8], block: &[u8; 64]) -> [u32; 8] {
    let mut w = [0u32; 64];
    for i in 0..16 {
        w[i] = u32::from_be_bytes(block[i * 4..(i + 1) * 4].try_into().unwrap());
    }
    for i in 16..64 {
        let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
        let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16]
            .wrapping_add(s0)
            .wrapping_add(w[i - 7])
            .wrapping_add(s1);
    }
    let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = state;
    for i in 0..64 {
        let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
        let ch = (e & f) ^ (!e & g);
        let temp1 = h
            .wrapping_add(s1)
            .wrapping_add(ch)
            .wrapping_add(K[i])
            .wrapping_add(w[i]);
        let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
        let mj = (a & b) ^ (a & c) ^ (b & c);
        let temp2 = s0.wrapping_add(mj);
        h = g;
        g = f;
        f = e;
        e = d.wrapping_add(temp1);
        d = c;
        c = b;
        b = a;
        a = temp1.wrapping_add(temp2);
    }
    [
        state[0].wrapping_add(a),
        state[1].wrapping_add(b),
        state[2].wrapping_add(c),
        state[3].wrapping_add(d),
        state[4].wrapping_add(e),
        state[5].wrapping_add(f),
        state[6].wrapping_add(g),
        state[7].wrapping_add(h),
    ]
}
