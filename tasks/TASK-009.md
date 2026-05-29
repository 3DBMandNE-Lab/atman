# TASK-009 — detectability: per-condition minimum + dead-code cleanup

<!-- kanban-status:start -->
**Status:** done
**Owner:** orchestrator
**Branch:** worktree-agent-a81e1655881c39b4a
**Started:** 2026-05-29T07:53:57+00:00
**Updated:** 2026-05-29T08:08:41+00:00
<!-- kanban-status:end -->

## Plan

1. Read the full detection + abundance flow in `detectability.rs` and the integration
   test to understand how per-condition detected counts (`n_detected_a`,
   `n_detected_b`) are tracked.
2. Add a per-condition minimum-detected guard to the abundance layer so a one-sided
   split (e.g. 3-vs-0) no longer yields an abundance effect.
3. Hoist the standard `Normal` out of the per-protein loop.
4. Look for and remove dead code.
5. Add a test covering the degenerate one-sided-detection case.

## Notes

- **Per-condition guard.** The detection test already required
  `n_a >= min_samples && n_b >= min_samples`, but the abundance layer only checked the
  global `y_abundance.len() < min_detected`. A protein detected in (say) 3 cases and 0
  controls cleared that global guard and the OLS rank guard, yet the condition
  coefficient had no between-group contrast. Added a new knob
  `--min-detected-per-condition` (default 2) and a branch: if
  `n_detected_a < min_detected_per_condition || n_detected_b < min_detected_per_condition`
  the abundance layer is suppressed and `skip_reason` is tagged `one_sided_detection`.
  Default 2 is the minimum that admits any between-group variance estimate while staying
  consistent with the existing `min_detected`/`min_samples` semantics. The knob is
  recorded in the run sidecar for provenance.
- **Standard normal hoisted.** `normal_two_sided_p` previously built
  `Normal::new(0.0, 1.0)` per call (once per protein in the detection branch). It now
  takes `&Normal`; the distribution is constructed once in `compute_rows` before the
  comparison loop and reused.
- **Dead code.** `cargo build --workspace` produces no warnings and there are no
  `#[allow(dead_code)]` shims in the file; the compiler flags nothing unused. The
  `row.clone()` before `detection_design.push` is necessary because `row` is
  subsequently moved into `abundance_design`. No removable dead code was found.
- **Determinism.** MNAR / below-LOD handling is unchanged (`detected_abundance` still
  treats `below_lod`/`dropped_by_qc`/non-finite as not detected — never imputed).
  Outputs change only for proteins newly suppressed by the per-condition guard: those
  rows now emit empty abundance columns and `one_sided_detection` instead of an OLS
  contrast computed against a single populated group. This is the intended behavioral
  change; all other rows are byte-identical.

## Verification

`cargo build --workspace` — clean, no warnings.

`cargo test -p atman --test detectability`:

```
running 2 tests
test detectability_suppresses_abundance_for_one_sided_detection ... ok
test detectability_writes_detection_and_abundance_layers ... ok

test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out
```

In-module unit tests (`logistic_recovers_positive_detection_shift`,
`design_rejects_condition_term`) also pass after the `normal_two_sided_p` signature
change (`cargo test -p atman detectability` → 2 passed in the unit-test binary).

Acceptance checklist:
- [x] abundance layer enforces a per-condition minimum detected count (3-vs-0 style split
      no longer yields an abundance effect)
- [x] standard Normal constructed once, not per row
- [x] a test covers the degenerate one-sided-detection case (abundance effect suppressed)
- [x] cargo test passes

## Hand-off

- New CLI flag `--min-detected-per-condition` (default 2) on `atman detectability`,
  documented in `docs/recipes.md` and recorded in the sidecar.
- Behavioral change: one-sided-detection proteins now report `one_sided_detection` and
  no abundance effect/p; downstream consumers keying on `skip_reason` should expect this
  new value.
- Risk: low. Default of 2 is conservative; users wanting the prior behavior can set
  `--min-detected-per-condition 0`. No core-crate changes.
