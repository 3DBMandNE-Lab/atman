# TASK-023 — post-hoc Sidak: t-based critical value for CI

<!-- kanban-status:start -->
**Status:** done
**Owner:** orchestrator
**Branch:** worktree-agent-ae80512b0311fe4b9
**Started:** 2026-05-29T21:42:28+00:00
**Updated:** 2026-05-29T21:46:55+00:00
<!-- kanban-status:end -->

## Plan

TASK-008 fixed the hardcoded `1.96` z-multiplier in the msqrob/Tukey/Dunnett
CIs but left the Sidak post-hoc path (`run_posthoc_sidak` in
`crates/atman/src/commands/de/posthoc.rs`) emitting a `1.96*se` CI while its
adjusted p (`1-(1-p)^m`) is computed on a t-based raw p. At small df the
z-based CI disagrees with the t-based p-value. Replace the `1.96` multiplier
with the per-contrast Student-t critical value at `fit.df`, mirroring the
established TASK-008 pattern (msqrob.rs L306-309: `StudentsT::new(0,1,fit.df)
.ok().map(|d| d.inverse_cdf(0.975)).unwrap_or(1.959963985)`). Update the
parity test's SE-recovery arithmetic in lockstep.

## Notes

- `run_posthoc_sidak`: added `use statrs::distribution::{ContinuousCDF,
  StudentsT}`; compute `t_crit` once per OLS fit (constant across the
  contrast family since all contrasts share `fit.df`); used it for
  `ci_low = est - t_crit*se` / `ci_high = est + t_crit*se`.
- The Sidak adjusted p (`1-(1-p)^m`) is UNCHANGED — only the CI critical
  value moved. raw p / adjusted-p / df / effect / t / mean_diff unchanged.
- Test `de_posthoc_sidak.rs`: SE was recovered as `(est-ci_low)/1.96`; now
  recovered as `(est-ci_low)/t_crit` where `t_crit` is the 97.5% Student-t
  quantile at the `df` already emitted on the row. Added the statrs import.
  The test still asserts on the unchanged posthoc_p / posthoc_adj_p / est;
  the `se` tolerance (1e-4) vs the R reference is unaffected because the
  half-width-to-SE conversion now uses the matching t-quantile.

## Verification

- `cargo build --workspace` — clean.
- `cargo test -p atman --test de_posthoc_sidak` — 2 passed.
- `cargo test -p atman --test de_posthoc_tukey` — 3 passed.
- `cargo test -p atman --test de_posthoc_dunnett` — 3 passed.

Acceptance checklist:
- [x] Sidak CI uses Student-t quantile at fit.df (not 1.96), matching the
      TASK-008 Tukey/Dunnett/msqrob fix.
- [x] Sidak adjusted p-value `1-(1-p)^m` unchanged.
- [x] `de_posthoc_sidak.rs` SE-recovery updated to the t-quantile; test still
      asserts on the unchanged adjusted p.
- [x] Only ci_low/ci_high move; p-values/adjusted-p/df/effect unchanged.
- [x] cargo test passes.

## Determinism

Only Sidak `ci_low`/`ci_high` byte-output changes (z→t critical value at
small df). All other columns of `de_results.tsv` and the run sidecar are
byte-identical. Tukey/Dunnett/msqrob paths untouched.

## Hand-off

Self-contained CI-consistency fix mirroring TASK-008. No new follow-ups; the
trio of OLS post-hoc CIs (Tukey/Dunnett/Sidak) plus msqrob now all use a
distribution-consistent critical value.
