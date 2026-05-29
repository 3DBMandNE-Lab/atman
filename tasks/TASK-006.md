# TASK-006 — meta: fix binomial sign-test trial count

<!-- kanban-status:start -->
**Status:** done
**Owner:** orchestrator
**Branch:** worktree-agent-ab33c1cabc856b588
**Started:** 2026-05-29T06:43:48+00:00
**Updated:** 2026-05-29T07:37:11+00:00
<!-- kanban-status:end -->

## Plan

1. Read `summarize_module`, `summarize`, and `binomial_upper_tail` to confirm the
   trial-count bug and how `n_positive`/`n_negative`/`effects` are computed.
2. Fix the trial count in `sign_binomial_p` to `n_positive + n_negative`.
3. Add a test with a zero-effect cohort proving the corrected p-value and that it
   differs from the old behavior.
4. Build the workspace and run `cargo test -p atman --test meta`.

## Notes

- The bug lives only in `summarize_module` (the module/`ModuleRow` path), which is
  the only summarizer that emits `sign_binomial_p`. The per-gene `summarize`
  (`MetaRow`) path computes `sign_consistency` only and has no binomial p, so it
  needed no change.
- `binomial_upper_tail(k, n)` sums `C(n,i)/2^n` for `i in k..=n`. Previously called
  with `k = max(n_positive, n_negative)` (zeros excluded) but `n = effects.len()`
  (zeros included). The fix sets `n = n_positive + n_negative`, so k and n are drawn
  from the same population (zero-effect cohorts excluded from both).
- Zero-effect cohorts (`point_delta == 0.0`) are retained by `read_module_effects`
  and still counted in `n_cohorts`; they are excluded only from the binomial sign
  test, which is the intended definition of a sign test.
- Behavior change: `sign_binomial_p` (and downstream BH q) changes for any program
  whose cohort set contains a zero effect — the p-value increases (less significant),
  correcting the previous deflation. No zero-effect cohorts => no change. No other
  meta outputs change.

## Verification

Acceptance criteria:
- [x] binomial sign test uses `n = n_positive + n_negative` as the number of trials
      (`crates/atman/src/commands/meta.rs:456`).
- [x] test with a zero-effect cohort confirms the corrected p-value and that it
      differs from the old behavior
      (`meta_module_sign_binomial_excludes_zero_effect_cohorts`): 3 positive + 1 zero
      cohort => corrected p = `binomial_upper_tail(3,3)` = 0.125, vs the old buggy
      `binomial_upper_tail(3,4)` = 0.3125; test asserts both.
- [x] cargo test passes.

Commands:
- `cargo build --workspace` => Finished (clean).
- `cargo test -p atman --test meta` => `test result: ok. 4 passed; 0 failed`.

## Hand-off

Single-line fix in `summarize_module` plus one new test. Reviewers verifying against
real meta runs: any module/program output containing a zero-effect cohort will show a
larger (correct) `sign_binomial_p` than before; this is the intended fix, not a
regression. Recompute any cached module-level meta fixtures/goldens if they pinned the
old `sign_binomial_p`/`bh_q` values.
