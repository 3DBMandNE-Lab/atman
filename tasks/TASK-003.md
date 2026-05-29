# TASK-003 — io: strict QC-flag parsing + unify float formatters

<!-- kanban-status:start -->
**Status:** review
**Owner:** orchestrator
**Branch:** worktree-agent-add24c59d666539a7
**Started:** 2026-05-29T06:28:26+00:00
**Updated:** 2026-05-29T07:40:54+00:00
<!-- kanban-status:end -->

## Plan

1. Make `parse_qc` strict: return `Result<QcFlag>`, `bail!` on any value other
   than the canonical `PASS`/`WARN`/`FAIL` spellings (the exact set emitted by
   `QcFlag::as_str`). Propagate the error from `read_measurements_long` with
   per-column context.
2. Collapse the two byte-identical formatters `format_f64`/`format_f64_fc`
   (both `format!("{}", v)`) into a single `format_f64`; update the fold-change
   panel caller.
3. Remove the local `format_float` shadow in `report.rs` (which was
   `format!("{v}")` — full precision, NOT the io 6dp `format_float`) and route
   its call sites through io's `format_f64`, which produces byte-identical
   output.
4. Route the bare-`{}` f64 writing in `enrich.rs` (GSEA + ORA writers) through
   `format_f64`. `{}`, `f64::to_string()`, and `format_f64` are all the same
   Display impl, so output bytes are unchanged.
5. Add precision-contract doc comments on both io formatters.

## Notes

- Byte-identity is preserved everywhere:
  - `format_f64` / `format_f64_fc` were literally identical bodies → safe merge.
  - report.rs local `format_float` = `format!("{v}")` = io `format_f64`. It is
    NOT io's `format_float` (which rounds to 6dp), so it was correctly swapped to
    `format_f64`, not to io's `format_float`. Swapping to io `format_float` would
    have changed bytes and was deliberately avoided.
  - enrich.rs fields (`es`, `nes`, `p_value`, `odds_ratio`, `bh_q`) all went from
    inline `{}` / `to_string()` to `format_f64` — same Display, same bytes.
    `bh_q.map(|q| q.to_string())` → `bh_q.map(format_f64)` keeps the
    empty-on-`None` behavior.
- io's `format_float` (fixed 6dp, NaN/Inf-aware) is a DIFFERENT precision
  contract and was left for its existing human-facing report/decompose columns;
  doc comment now states this divergence explicitly. The two contracts are not
  unified because doing so would change output bytes.
- Added two unit tests locking the strict `parse_qc` behavior.

## Verification

- `cargo build --workspace`: clean, no warnings.
- `cargo test -p atman`: 209 passed, 0 failed (across all integration suites)
  plus `--lib` 36 passed, 0 failed including the two new `parse_qc` tests.
- Guardrail / determinism tests pass:
  - `determinism_new_methods` — pass
  - `report_qc` (report_qc_writes_summary_sample_protein_and_condition_tables) — pass
  - `enrich_gsea` incl. `enrich_gsea_two_runs_same_seed_match_byte_for_byte` — pass
  - `enrich_ora` — pass
  - io round-trip tests (`de_results_round_trip_preserves_limma_columns`,
    `de_ensemble_tsv_roundtrips_all_columns`, etc.) — pass
- New tests:
  - `parse_qc_accepts_canonical_spellings` — pass
  - `parse_qc_rejects_unrecognized_values` (lowercase/truncated/empty/trailing-space) — pass
- `grep` confirms no remaining references to `format_f64_fc` and no leftover
  `format_float` shadow in report.rs.

Acceptance criteria:
- [x] parse_qc errors on any non-canonical flag
- [x] format_f64 / format_f64_fc collapsed into one
- [x] report.rs / enrich.rs no longer define/use a replaceable shadow formatter
- [x] no change to output bytes for valid inputs (determinism tests green)
- [x] cargo test passes

## Hand-off

Code change is isolated to `crates/atman/src/io.rs`,
`crates/atman/src/commands/report.rs`, `crates/atman/src/commands/enrich.rs`.
No public output format changes. The two distinct float-precision contracts
(`format_f64` full round-trip vs `format_float` fixed-6dp) remain intentionally
separate and are now documented; a future unification would require deciding
which columns change bytes and is out of scope here.

### Codex adversarial-review fix (orchestrator)

Codex flagged (HIGH): the strict `parse_qc` rejected empty QC cells, but empty is a previously-valid "no flag" form that `validate_measurements` (validate.rs:360) still accepts as `"PASS" | "WARN" | "FAIL" | ""` — an over-strict regression that made the typed reader and the raw validator disagree. Fixed: `parse_qc` now maps `""` to `Pass` (the documented unspecified default) and bails only on *non-empty* unrecognized values, so garbage/typos are still rejected. Updated the unit tests accordingly (empty accepted; "pass"/"fail"/"FA"/"OK"/"PASS " rejected). `report_qc`, `determinism_new_methods`, `enrich_gsea` all still green.
