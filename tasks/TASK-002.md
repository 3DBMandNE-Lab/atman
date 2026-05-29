# TASK-002 — core: guard effective_abundance() against non-finite values

<!-- kanban-status:start -->
**Status:** done
**Owner:** orchestrator
**Branch:** worktree-agent-a9a61d25457420390
**Started:** 2026-05-29T06:28:26+00:00
**Updated:** 2026-05-29T07:37:11+00:00
<!-- kanban-status:end -->

## Plan

1. Read both `effective_abundance` impls plus `MeasurementRecord` / `Abundance`
   types to match conventions.
2. Add a finiteness guard to `MeasurementRecord::effective_abundance()` so it
   returns `None` for a non-finite abundance on a non-dropped row, mirroring
   `PeptideMeasurementRecord::effective_abundance()`.
3. Add a unit test covering NaN / +inf / -inf on a non-dropped row.
4. Build the workspace and run the test suites; confirm no regressions.

## Notes

- Fix in `crates/atman-core/src/types.rs` (`MeasurementRecord::effective_abundance`):
  now computes `value = self.abundance.as_f64()` and returns `None` when
  `self.dropped_by_qc || !value.is_finite()`. This matches the peptide twin
  exactly. `Abundance::as_f64()` simply unwraps the inner f64 regardless of
  unit, so guarding the unwrapped value is correct for all variants.
- Caller audit: every consumer treats `None` as "skip this measurement"
  (`module_de`, `bootstrap`, `ratio`, `residuals`, `fold_change`, `de/*`,
  `null`, `score`, `report`, `validate`). The change only enlarges the skip
  set for malformed (NaN/inf) input; no caller relied on the old
  non-finite-passing behavior.
- `within_cohort_rank.rs:81` already double-guarded with
  `.is_some_and(f64::is_finite)`; that guard is now redundant but harmless and
  was left unchanged.
- No existing test encoded the buggy behavior, so none required fixing.

## Verification

Acceptance criteria:
- [x] `effective_abundance()` returns `None` for a non-finite abundance on a
      non-dropped row.
- [x] Behavior matches `PeptideMeasurementRecord::effective_abundance`
      (same `dropped_by_qc || !value.is_finite()` predicate).
- [x] Unit test covers NaN/inf:
      `effective_abundance_rejects_non_finite_on_non_dropped_row` in
      `crates/atman-core/tests/types_roundtrip.rs` (NaN, +inf, -inf all → None;
      finite value still passes).
- [x] `cargo test` passes.

Commands:
- `cargo build --workspace` → Finished, no errors.
- `cargo test -p atman-core` → all pass; `types_roundtrip` now 6 tests incl.
  the new non-finite test.
- `cargo test -p atman` → all suites pass, 0 failures.

Determinism: change affects only the malformed non-finite edge case (NaN/inf
on a non-dropped row). No existing valid-input test output changed.

## Hand-off

Single-line semantic fix plus a unit test. Ready for review. No follow-ups
required; the redundant `is_some_and(f64::is_finite)` guard in
`within_cohort_rank.rs` could be removed later as a cleanup but is harmless.
