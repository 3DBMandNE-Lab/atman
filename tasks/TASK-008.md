# TASK-008 — post-hoc/msqrob: use correct critical values for CIs

<!-- kanban-status:start -->
**Status:** review
**Owner:** orchestrator
**Branch:** worktree-agent-a93ede88c16a6dbb4
**Started:** 2026-05-29T07:53:57+00:00
**Updated:** 2026-05-29T08:03:05+00:00
<!-- kanban-status:end -->

## Plan

1. Read the correct model (`limma.rs` CI path), the `msqrob.rs` path, and the
   posthoc tukey/dunnett paths. Inventory quantile helpers in `atman-core`.
2. msqrob: replace the hardcoded `1.96` z-multiplier with the two-sided 97.5%
   Student-t quantile at `fit.df`, mirroring limma.
3. posthoc: emit simultaneity-adjusted CIs using the SAME distribution as each
   test's adjusted p-value (studentized-range for Tukey, Dunnett for Dunnett),
   since both inverse quantile functions exist in core.
4. Remove the dead scaffolding loop in the dunnett path.
5. Build + run targeted tests; fix any tests that hard-code the old 1.96 CI.

## Notes

Quantile helpers available in `atman-core` (decided to compute exact adjusted
CIs rather than document them as unadjusted):
- `studentized_range::qtukey(p, n, df)` — inverse of `ptukey`.
- `multivariate_t::qdunnett(p, m, df, rho)` — inverse of `pdunnett` (balanced).
- `multivariate_t::pdunnett_hsu(q, df, corr, n_mc, seed)` — Hsu CDF (unbalanced);
  no closed-form inverse exists, but it is deterministic (fixed seed/draws), so
  the critical value is recovered by bisecting the SAME MC CDF used for the
  adjusted p. This is the literal inverse of the p-value's CDF — not a
  fabricated quantile.

Changes:
- `crates/atman/src/commands/de/msqrob.rs`:
  - Added `use statrs::distribution::{ContinuousCDF, StudentsT};`.
  - CI now uses `StudentsT::new(0,1,fit.df).inverse_cdf(0.975)` (fallback
    1.959963985 only if the distribution can't be constructed), not 1.96.
    Mirrors `limma.rs:357`. p-value unchanged.
- `crates/atman/src/commands/de/posthoc.rs`:
  - Tukey: added free fn `tukey_ci_halfwidth(est, se, k, df)` = `qtukey(0.95,
    k, df)·se/√2` (the studentized-range stat is `diff·√2/se`, so the
    estimate-scale half-width divides by √2). CI now uses this; the interval
    excludes 0 iff the Tukey adjusted p < 0.05.
  - Dunnett: added a `dunnett_crit(df)` closure — `qdunnett(0.95, m, df, RHO)`
    for balanced designs, bisection of `pdunnett_hsu` (same seed/draws/matrix
    as `dunnett_cdf`) for the Hsu/unbalanced branch. CI half-width = `c·se`.
  - Removed the dead `for (sid, _) in &design_rows { let _ = sid; }` scaffolding
    loop (posthoc.rs ~1055); folded its comment into the live `samples` loop.
  - Sidak path left untouched — out of TASK-008 scope (see Hand-off).
- `crates/atman/tests/de_posthoc_tukey.rs`: the self-consistency test recovered
  SE from the CI assuming `ci_low = est − 1.96·se`; updated it to the new Tukey
  convention `se = (est − ci_low)·√2 / qtukey(0.95, 3, df)`. This test asserts
  on adjusted p (which did NOT change); only its SE-recovery arithmetic needed
  to track the intended CI change.

Determinism: CI columns (`ci_low`/`ci_high`) for msqrob, tukey, and dunnett
change where the critical value changed — intended and documented. No p-value,
adjusted-p, t, df, or effect column changes. `determinism_new_methods.rs` does
not hash any posthoc/msqrob output, so no golden fixtures shifted.

## Verification

Acceptance checklist:
- [x] msqrob CI uses the t-quantile at `fit.df`, not 1.96 (msqrob.rs).
- [x] post-hoc CIs use the matching simultaneity-adjusted critical value:
      Tukey via `qtukey`, Dunnett via `qdunnett` (balanced) / bisected
      `pdunnett_hsu` (unbalanced). Implemented exact adjusted CIs because the
      inverse quantiles exist in core; nothing fabricated.
- [x] dead scaffolding loop at posthoc.rs ~1055 removed.
- [x] cargo test passes.

Commands:
- `cargo build --workspace` → Finished, no errors.
- `cargo test -p atman --test de_msqrob --test de_msqrob_cptac --test
  de_posthoc_tukey --test de_posthoc_dunnett --test de_posthoc_sidak`:
  - de_msqrob: 4 passed
  - de_msqrob_cptac: 1 passed
  - de_posthoc_tukey: 3 passed (self-consistency test now green after the
    SE-recovery fix)
  - de_posthoc_dunnett: 3 passed (incl. balanced + Hsu/unbalanced switch)
  - de_posthoc_sidak: 2 passed
- `cargo test -p atman` (full suite): 0 failed.
- `cargo clippy -p atman`: no new warnings on the changed paths (dead loop
  gone, no unused imports).

## Hand-off

Follow-up / risk: the **Sidak** post-hoc path (`run_posthoc_sidak`, posthoc.rs
~316) still emits a `1.96·se` CI while its adjusted p is `1 − (1 − p)^m` on a
t-based raw p — the same class of inconsistency, but outside TASK-008's named
scope (task specified msqrob + tukey + dunnett only). A natural fix is a
per-contrast t-interval at `fit.df` (raw, unadjusted), or a Šidák-adjusted
two-sided t critical value `qt(1 − α_sidak/2, df)` with `α_sidak = 1 − (1 −
α)^(1/m)` for a simultaneity-consistent interval. `de_posthoc_sidak.rs:130`
also recovers SE via `/1.96`, so that test must be updated in lockstep.
Recommend a separate task.
