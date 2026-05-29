# TASK-021 — split god-module de/mod.rs (2,507 LOC)

<!-- kanban-status:start -->
**Status:** review
**Owner:** orchestrator
**Branch:** worktree-agent-a5f9a83a1142451a5
**Started:** 2026-05-29T21:07:03+00:00
**Updated:** 2026-05-29T21:33:41Z
<!-- kanban-status:end -->

## Plan

Pure refactor of the ~2,507-LOC `de/mod.rs` god-module. Move cohesive,
self-contained pieces verbatim into three new submodules under `de/`,
leaving `mod.rs` as the dispatch + thin orchestration (the `run` entry,
validation, per-method routing, result-row assembly, sidecar write).

New submodules:
- `args.rs` — the clap `Args` struct + `apply_per_subject_proxy` /
  `canonical_proxy_name` arg preprocessing.
- `ols_design.rs` — design-matrix construction and OLS plumbing:
  `DesignTerm`, `OlsSetup`, `OlsDesign`, `CovKind`, `build_ols_setup`,
  `parse_design_terms`, `covariate_names_from_terms`, `resolve_comparisons`,
  `validate_random_intercept_design`, `read_covariate_columns`,
  `classify_covariates`, `build_ols_design`, `push_encoded_covariate`,
  `ensure_full_rank`, `matrix_rank`, `raw_group_means`,
  `ols_to_paired_t_result`, `build_two_group_design`,
  `augment_ols_design_with_external`, `paired_t_covariate_adjusted`.
- `output_rows.rs` — sidecar row shapes + TSV writers (`OmnibusRow`,
  `CovariateRow`, `DesignReportRow`, `ReportAccumulator`,
  `write_omnibus_rows`, `write_covariate_rows`, `write_proxy_summary`,
  `format_opt`, `fmt_opt`) and `override_subject_id`.

## Notes

- Strictly a wiring change: code moved verbatim. No logic, ordering,
  formatting, FDR/stat, or provenance content changed.
- `commands::de::run` and the module path are preserved; `Args` is now
  `pub use args::Args` so the external CLI surface is unchanged.
- Visibility: items previously bare-private in `mod.rs` but referenced by
  sibling submodules became `pub(super)` in their new home. `mod.rs`
  re-exports the shared surface via plain `use` (private items in a parent
  module are visible to descendant submodules), so existing
  `use super::{...}` in `limma.rs`/`msqrob.rs`/`ensemble.rs`/`posthoc.rs`
  keep resolving unchanged.
- `CovKind::col_labels` promoted from private to `pub(super)` (now called
  cross-module from `posthoc.rs`, previously a descendant of `mod.rs`).

## Resulting layout (LOC)

```
mod.rs        1154  (was 2507; dispatch + run orchestration)
args.rs        296  (new)
ols_design.rs  934  (new)
output_rows.rs 196  (new)
```
Other pre-existing de/ submodules unchanged (ensemble 301, limma 560,
moderated 76, msqrob 469, posthoc 1520, preflight 91, robust_stats 213).

## Verification

- `cargo build --workspace` — clean, no warnings.
- `cargo clippy -p atman` — no new warnings; none point at the new
  submodules or `de/mod.rs` (the 36 warnings are pre-existing, in
  unrelated modules: NMF/decompose loops, etc.).
- `cargo test -p atman` — 218 passed, 0 failed. Full de_* suite ran:
  de_adjust_for, de_below_lod_gate, de_deqms_cptac, de_ensemble,
  de_ensemble_dube, de_limma, de_mixed, de_msqrob, de_msqrob_cptac,
  de_ols_formula, de_omnibus_factor, de_paired_by, de_posthoc_dunnett,
  de_posthoc_sidak, de_posthoc_tukey, de_robust_stats, dube_de_sanity,
  dube_reproduction — plus determinism_new_methods (byte-identity).
- `cargo test --workspace` — 483 passed, 0 failed.
- Byte-identical output preserved: determinism_new_methods + golden de_*
  tests pass unchanged.

## Hand-off

Ready for review. Pure refactor; output byte-identical, all tests green.
Only `de/mod.rs` (modified) + `args.rs`/`ols_design.rs`/`output_rows.rs`
(new) touched. No CLI/behavior change.
