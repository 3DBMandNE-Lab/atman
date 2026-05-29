# TASK-015 — CI: build the real toolchain + enforce determinism pinning

<!-- kanban-status:start -->
**Status:** review
**Owner:** orchestrator
**Branch:** (pending)
**Started:** 2026-05-29T09:49:45+00:00
**Updated:** 2026-05-29T10:18:00+00:00
<!-- kanban-status:end -->

## Plan

1. Confirm version numbers: MSRV (`Cargo.toml` rust-version) and `Dockerfile` both pin 1.94; CI pinned 1.75.
2. Bump CI `toolchain:` 1.75 -> 1.94 in both jobs that set it (lint, test).
3. Add `rust-toolchain.toml` at repo root pinning channel 1.94 (with rustfmt, clippy).
4. Set BLAS/OMP/OPENBLAS thread-pinning env at the test job level in ci.yml.
5. Validate YAML, rebuild workspace, run determinism test.

## Notes

- Confirmed: `Cargo.toml` `rust-version = "1.94"`, `Dockerfile` `FROM rust:1.94-slim-bookworm`. CI was 1.75 in both `lint` and `test` jobs — drift fixed.
- Did NOT touch the stray top-level `github/` directory (separate cleanup task).
- `rust-toolchain.toml` channel set to "1.94" to match MSRV/Dockerfile and prevent CI/Docker/local drift.
- Test-job thread-pinning env mirrors what `determinism_new_methods.rs` sets internally, so the byte-identity guarantee is now enforced at the documented job level too.

## Verification

- [x] CI toolchain set to 1.94 in all jobs (lint + test) — confirmed via YAML parse.
- [x] rust-toolchain.toml committed at repo root pinning 1.94.
- [x] BLAS/OMP/OPENBLAS thread-pinning env on the test job in ci.yml — confirmed via YAML parse.
- [x] ci.yml is valid YAML (parsed with PyYAML: `YAML OK`; test env and both toolchain fields verified = 1.94).
- [x] `cargo build --workspace` succeeds locally (dev profile, 19.10s; local rustc is Homebrew 1.94.1 matching the pin).
- [x] `cargo test --workspace --release --test determinism_new_methods` — 3 passed, 0 failed.

## Hand-off

CI now builds the supported toolchain (1.94) on all jobs and pins BLAS/OMP/OPENBLAS=1 on the
test job so the byte-identity determinism claim is enforced under documented conditions.
`rust-toolchain.toml` prevents future drift across CI/Docker/local.
