# TASK-007 — limma: BH-adjust the per-feature F-test p-values

<!-- kanban-status:start -->
**Status:** review
**Owner:** orchestrator
**Branch:** (pending)
**Started:** 2026-05-29T07:53:57+00:00
**Updated:** 2026-05-29T07:55:53Z
<!-- kanban-status:end -->

## Plan

Populate `f_bh_q` (previously hardcoded `None` at line 508) by running the
existing `bh_fdr` helper over the finite per-feature F-test p-values, grouped
per (comparison, panel), mirroring exactly how `p_value` → `bh_q` is computed
at lines 345-351. Chose the "add BH adjustment" path (not column removal) for
consistency: `f_p_value` is a reportable per-protein hypothesis and was the
only p-column lacking multiplicity control.

## Notes

- Added an `f_q_values` computation block immediately after the existing
  `q_values` block in `crates/atman/src/commands/de/limma.rs`. It builds
  `fp_by_panel: BTreeMap<String, Vec<Option<f64>>>` in the same
  `panel_by_feature` iteration order used to build `panel_positions`, so the
  splice-back via `panel_positions[panel]` aligns 1:1 with the `p_value` path.
- Finiteness guard mirrors the `p_value` path: a feature contributes its
  F p-value only when `!row.skipped && f_p_value.is_some() && finite`,
  otherwise `None` (so `bh_fdr` excludes it from the denominator, same as
  `p_value`).
- Reused `atman_core::de::bh_fdr` — no reimplementation.
- Set `f_bh_q: f_q_values[i]` at the `DeResultRow` construction (was `None`).
- Determinism: `f_bh_q` gains values only where it was empty; the BH order is
  deterministic via `BTreeMap` panel iteration. No other column changes.
- CLI currently emits one contrast per comparison, so `f_p_value` is typically
  `None` for default runs; `f_bh_q` therefore stays empty in those cases
  (correct). It populates once a k>=2 contrast path produces finite F p-values.

## Verification

- `cargo build --workspace` — clean (Finished dev profile).
- `cargo test -p atman --test de_limma` — 7 passed; 0 failed. Includes
  `multi_group_f_test_populates_f_columns` (asserts f_statistic/f_p_value/
  f_bh_q header presence) and `limma_is_deterministic_across_runs`.
- `cargo test -p atman --test determinism_new_methods` — 3 passed; 0 failed
  (de_adjust_for, ica_missingness, nmf byte-determinism all green).

Acceptance checklist:
- [x] f_bh_q populated via `bh_fdr` over finite f_p_values per (comparison, panel)
- [x] uses the same helper/path as p_value (BTreeMap per-panel + panel_positions splice)
- [x] cargo test passes

## Hand-off

Single-file change in `crates/atman/src/commands/de/limma.rs`. No schema, no
new dependencies. f_bh_q now carries BH-FDR multiplicity control consistent
with the primary `p_value`/`bh_q` column.
