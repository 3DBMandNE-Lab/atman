# TASK-001 — ingest-matrix: preserve raw abundance + fail loud on dropped rows/cols

<!-- kanban-status:start -->
**Status:** review
**Owner:** orchestrator
**Branch:** (pending)
**Started:** 2026-05-29T06:28:26+00:00
**Updated:** 2026-05-29T06:45:00Z
<!-- kanban-status:end -->

## Plan

1. Add a `abundance_raw: Option<f64>` field to `MeasurementRow`, set once at construction to the pre-normalization working value and never mutated by `apply_normalization`.
2. Emit that field in the `abundance_raw` column of `measurements.tsv` (was incorrectly emitting the normalized `abundance`).
3. Strict-failure: warn on matrix value columns matching no sample id; hard-bail on a proteins/rows length mismatch that would silently truncate.
4. Add a fidelity test proving `abundance_raw` survives normalization unchanged.

## Notes

- `measurement_row()` sets `abundance_raw: abundance` at the same point it sets `abundance` (both = the parsed, optionally-log2-transformed input). `apply_normalization` only mutates `r.abundance` (confirmed: `crates/atman/src/commands/ingest_matrix.rs` `apply_normalization` assigns `r.abundance.as_mut()` exclusively), so `abundance_raw` retains the pre-normalization value.
- The verbatim pre-transform input string remains in the `source` (`npx_source_str`) column; `abundance_raw` is the pre-normalization value at the declared `abundance_unit` scale (i.e. post `--log2-transform`). Documented this on the struct field.
- Unmatched matrix columns: previously only the all-columns-unmatched case bailed; individual stray columns were silently dropped. Now warns loudly listing them.
- proteins/rows length mismatch: `proteins.iter().zip(&matrix.rows)` would silently truncate; now bails with counts.
- Recovered from an orchestrator-side crash (subagent socket drop after implementing the fix but before committing/testing). The implementation was complete and correct; the orchestrator added the fidelity test, verified, and committed.

## Verification

| Acceptance criterion | Result |
|---|---|
| abundance_raw holds the pre-transform, pre-normalization input value | PASS — emitted from the unmutated `abundance_raw` field |
| round-trip read shows abundance_raw != normalized abundance when normalized | PASS — new test `ingest_matrix_abundance_raw_is_pre_normalization` asserts S1 raw = {14,16,18,20} while normalized differs |
| unmatched matrix sample columns produce a warning | PASS — `eprintln!` WARNING listing unmatched columns |
| length mismatch detected, not silently truncated | PASS — `bail!` with row/protein counts |
| cargo test passes; a test asserts abundance_raw fidelity | PASS — `cargo test -p atman --test ingest_matrix` → 7 passed; 0 failed |

```
running 7 tests
test ingest_matrix_abundance_raw_is_pre_normalization ... ok
...
test result: ok. 7 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out
```

## Hand-off

- `abundance_raw` is now a faithful provenance column. Any downstream QC that compares raw vs effective abundance will work as intended.
- Determinism: the `abundance_raw` column bytes change (that is the fix); no other output changes.
- Follow-up (out of scope): the `read_measurements_long` round-trip backfill of an empty `abundance` from `abundance_raw` (io.rs) is a separate concern tracked under the io review findings.
