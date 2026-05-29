# TASK-005 — ensemble: stop treating correlated method p-values as independent

<!-- kanban-status:start -->
**Status:** done
**Owner:** orchestrator
**Branch:** worktree-agent-a6183a81b19082def
**Started:** 2026-05-29T06:43:48+00:00
**Updated:** 2026-05-29T07:53:20+00:00
<!-- kanban-status:end -->

## Plan

1. Read `de/ensemble.rs`, `atman_core::ensemble`, and both integration tests to
   understand how VALIDATED is produced today.
2. Confirm the grade currently hard-gates on `ensemble_q` (BH of the
   independence-assuming Stouffer p) AND sign consistency, while `n_significant`
   (per-method BH-q hits) is computed but unused for grading.
3. Apply Option A: re-derive the grade from method agreement
   (`n_significant`/`n_applied`) + sign consistency; demote
   `ensemble_p`/`ensemble_q` to a documented heuristic that no longer feeds the
   grade.
4. Document `ensemble_p` semantics in code (`atman_core::ensemble`, command
   module) and in `docs/recipes.md`.
5. Dedup the duplicated `"measurements.tsv"` in the provenance hash list.
6. Update unit + integration tests to the corrected (more conservative) grades;
   document which grades changed and why.
7. Build + test; verify byte-identity of non-`ensemble_p` columns.

## Notes

**Option A chosen** (per card preference). The grade already had the
independence-free ingredients: `n_significant` (count of methods whose own
per-method BH-q clears the threshold) and sign consistency
(`n_sign_consistent`/`n_applied`). The old `assign_grade` instead hard-gated on
`ensemble_q`, the BH-FDR of the Stouffer-combined p. Stouffer assumes
independence; every ensemble method is fit on the same abundance matrix, so the
per-method p-values are positively correlated and the combined p is
anti-conservative — it manufactured significance the per-method evidence did not
support, then fed VALIDATED.

New `assign_grade(n_significant, n_applied, n_sign_consistent, thresholds)`:
- INSUFFICIENT when `n_applied == 0`.
- VALIDATED when BOTH the significant fraction `n_significant/n_applied` AND the
  sign-consistency fraction reach `validated_sign_fraction` (default 1.00 —
  every method individually significant and agreeing on direction).
- PROVISIONAL when BOTH fractions reach `provisional_sign_fraction`
  (default 0.50 — a majority significant and agreeing).
- Otherwise INSUFFICIENT.

`ensemble_p`/`ensemble_q` are retained as a ranking-convenience heuristic only.
Documented as NOT a calibrated p-value in: `combine_stouffer` doc, module-level
docs of `atman_core::ensemble`, the command module header, and `docs/recipes.md`.

**Grade outputs that changed (expected, documented):**
- Dube PT2-PR2 cohort: canonical HSPs (HSPA1A, DNAJB1, HSPB1) move VALIDATED →
  INSUFFICIENT. Inspected the rows: all three have `n_significant = 0` — NO
  single method clears its own per-method BH-q across ~2938 proteins (n≈9–20
  paired subjects). Their old VALIDATED came solely from the anti-conservative
  combined `ensemble_q` (~0.001–0.02). Sign is consistently positive (PT2 > PR2,
  biologically correct), which the test still asserts. Across the whole PT2-PR2
  comparison, the honest grading yields 2 PROVISIONAL and 0 VALIDATED (vs. the
  old anti-conservative spread). This is the intended correction, not a
  regression.
- Synthetic 3-protein fixture (`de_ensemble.rs`): strong UP/DN effects are
  significant under both methods individually, so they remain VALIDATED; null ST
  remains INSUFFICIENT. No change.

**Provenance fix:** removed the duplicate `"measurements.tsv"` from the
`hash_canonical_inputs` argument in `de/ensemble.rs` (was listed twice).

**Pre-existing non-determinism flagged (out of scope, NOT introduced here):**
`ensemble_p`/`ensemble_q` differ by floating-point ULPs run-to-run. Verified by
`git stash` that the UNMODIFIED code has the same instability. All other columns
(grade, counts, signs, ordering) are byte-identical across runs. Root cause is
upstream of this task (z-score accumulation order in the combined p). Follow-up
candidate, but since those two columns are now demoted to a non-load-bearing
heuristic, the grade/decision output is fully deterministic.

## Verification

Acceptance criteria:
- [x] VALIDATED no longer relies on an independence-assuming combined p —
  `assign_grade` now takes `n_significant` (per-method agreement) + sign
  consistency; `ensemble_q` is not read by the grader.
- [x] `ensemble_p` semantics documented in code (module docs,
  `combine_stouffer`, command header) and in `docs/recipes.md`.
- [x] Duplicate `measurements.tsv` removed from the provenance hash list.
- [x] Tests pass.

Commands:
- `cargo build --workspace` → Finished, clean.
- `cargo test -p atman-core --lib ensemble` → 11 passed, 0 failed.
- `cargo test -p atman --test de_ensemble` → 2 passed, 0 failed.
- `cargo test -p atman --test de_ensemble_dube` → 1 passed, 0 failed.
- `cargo test -p atman-core` → 238 + (4/3/2/5/5) passed, 0 failed.
- `cargo test -p atman` (full suite) → all green, 0 failures.
- `cargo clippy -p atman-core -p atman` → no ensemble findings.

Determinism: re-ran Dube ensemble twice; diff excluding `ensemble_p`/`ensemble_q`
columns is empty (grade/counts/signs/order byte-identical). The two heuristic
columns differ by ULPs — confirmed pre-existing via `git stash` on unchanged
code.

## Hand-off

- Grade semantics changed: VALIDATED/PROVISIONAL now require per-method
  significance + sign agreement, not the combined p. Downstream consumers that
  read `de_ensemble.tsv` grades will see a more conservative, honest call set.
- `ensemble_p`/`ensemble_q` remain in the schema but are heuristic ranking
  columns only — do not treat as calibrated significance.
- Follow-up (optional, separate task): make `ensemble_p` accumulation
  order-stable to restore full byte-identity on those two columns.

### Codex adversarial-review fix (orchestrator)

Codex flagged (LOW): `docs/analytical-roadmap.md` §16 still described the OLD Stouffer-`ensemble_q`-gated grade and claimed the Dube HSPs were all VALIDATED — contradicting the implemented honest grade. Updated §16: the description, the grading-logic bullets (now anchored to `n_significant`/`n_applied` + sign fraction, not the combined p), the `ensemble_p`/`ensemble_q` output description (marked non-calibrated heuristic), and the Dube acceptance bullet (HSPs show consistent sign but are NOT VALIDATED at this n because no method clears its own BH-q). No code change.
