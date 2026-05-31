# Atman 1.1.0 — pre-release checklist

Status of the gating items between the cut `v1.1.0` tag and a confident public
release. `[x]` = done in this repo; `[ ]` = needs your action, a decision, or
another machine. See `docs/release-checklist.md` for the generic publish flow.

## 1. Done in this pass (verify, don't redo)

- [x] Version bumped to 1.1.0 — `Cargo.toml`, `crates/atman/Cargo.toml` pin,
      `Cargo.lock`, `CITATION.cff` (+ date 2026-05-30), README docker tags.
- [x] `CHANGELOG.md` `[1.1.0]` section dated; fresh `[Unreleased]` opened.
- [x] Annotated tag `v1.1.0` created (local).
- [x] Kanban orchestration scaffolding removed from the product tree.
- [x] Docs reconciled with code (README command list, ensemble grading,
      recover-plex `--infer-pairs`, detectability flags).
- [x] Full suite green locally (490 passed / 0 failed) and the byte-identity
      determinism gates pass in release profile — on **darwin only**.

## 2. Local verification — DONE (2026-05-30)

- [x] `cargo fmt --all -- --check` — clean (the tree was not fmt-clean; ran
      `cargo fmt --all`, layout-only, no behavior change).
- [x] `cargo clippy --workspace --all-targets -- -D warnings` — **clean (exit 0)**.
      ~40 pre-existing style lints (28 in `atman-core` numeric kernels, ~12 in
      the command layer + test targets) were fixed: real iterator/slice/closure
      rewrites where behavior-identical, scoped `#[allow(...)]` with comments on
      load-bearing lockstep-index loops. Verified byte-identical via the
      determinism + golden suites; reviewed for numeric preservation.
- [x] `cargo test --workspace --release` — 490 passed, 0 failed.
- [x] `cargo package -p atman-core` — clean. `cargo package -p atman` fails
      only on the documented publish-order constraint (needs `atman-core` on
      crates.io first); moot for a local/non-crates.io release.
- [x] `docker build -t atman:1.1.0 .` + `--version`/`--help` — success. Builds
      on `rust:1.94-slim-bookworm`, so this also validated a **Linux release
      build** (see §3).

## 3. Cross-platform / CI (needs another machine or a remote)

This is the biggest open risk: everything so far was validated on one darwin
host, and the repo is local-only so CI has never actually run.

- [ ] Build + test on **Linux** (the CI target / Dockerfile base, rust 1.94).
- [ ] Build + test on **Windows** if it's a supported target (the build.rs
      git-dir resolution and BLAS/thread-pinning are the likely sore spots).
- [ ] If/when a remote is added: push and confirm `.github/workflows/ci.yml`
      goes green with the 1.94 toolchain + thread-pinning env.

## 4. Decisions only you can make — LARGELY DONE (2026-05-30)

- [x] **Docs verified against code** for the four new commands (`detectability`,
      `absence-topology`, `recover-plex`, `robust-paired`). Every documented
      claim, default, output column, and skip/fail condition was checked. The
      one overclaim found is fixed (detectability does a logistic regression, not
      a "contingency" test). The non-blocking polish is also done: recipes.md now
      documents the previously-omitted flags (`detectability`
      `--q-threshold`/`--design`/`--max-iter`/`--tol`, `recover-plex`
      `--expected-cluster-size`, `robust-paired --min-pairs`), the missing
      diagnostic values (`single_cluster`, `no_clusterable_proteins`), the
      stem-derived `*_quality.tsv` / `*_pairs.tsv` / `*_summary.tsv` filenames,
      the `--infer-pairs` numeric-stem requirement + synthetic `pair_NNNN` ids,
      and the strict-failure hard-error paths. **Your remaining call:** confirm
      the descriptions match *intent* (the precision is now verified; the science
      is yours).
- [x] **Behavior changes release-noted** — CHANGELOG `[1.1.0]` Changed/Fixed
      sections cover RNG output changes (regenerate cached outputs), strict QC,
      `recover-plex --infer-pairs` hard-fail, and the stricter ensemble grade;
      trust-boundary changes are in `SECURITY.md`.
- [x] **clippy denial policy** — resolved by making `-D warnings` clean (§2).

## 5. Security / supply chain — DONE (2026-05-30)

- [x] **Security review** of input parsing / paths / network / `run`. No
      CRITICAL/HIGH; no `unsafe`. Fixed the one MEDIUM (data-derived output
      filename path-traversal in `write_wide_panel`/`write_fold_change_panel`,
      now sanitized + regression-tested). Trust model documented in `SECURITY.md`.
- [x] **Adversarial-input robustness** — `crates/atman/tests/adversarial_input.rs`
      (ragged/missing-column/binary/truncated/empty → graceful error, no panic),
      plus a `fuzz/` cargo-fuzz scaffold for the TSV readers. **Remaining:** run
      `cargo +nightly fuzz run <target>` in CI (needs nightly; see `fuzz/README.md`).
- [x] **SBOM + license audit** — `docs/sbom-1.1.0.cdx.json` (CycloneDX,
      179 components); all transitive licenses permissive (MIT/Apache/BSD/ISC/
      Zlib/Unicode), no copyleft. `cargo audit` run.
- [x] **`README.pdf` / `mdpdf.log`** confirmed untracked + gitignored — not shipped.
- [x] **`cargo audit` gate in place.** CI runs `cargo audit` (the `audit` job).
      The two transitive-and-unreachable advisories (RUSTSEC-2026-0097 `rand`,
      RUSTSEC-2024-0436 `paste`, via `statrs → nalgebra`) are scoped in
      `.cargo/audit.toml` with justification (any *new* advisory still fails).
      The gate also caught a **reachable** one — RUSTSEC-2026-0104
      (`rustls-webpki` CRL-parsing panic on the gprofiler TLS path) — now fixed
      by bumping `rustls-webpki` to 0.103.13. `cargo audit` is clean.

## 6. Publish mechanics (when 2–4 are satisfied)

- [ ] `git push` + `git push --tags` (once a remote exists).
- [ ] Create the GitHub release for `v1.1.0` (paste the CHANGELOG section).
- [ ] If publishing to crates.io: `cargo publish -p atman-core` first, then
      `cargo publish -p atman` (the CLI depends on the published core).
