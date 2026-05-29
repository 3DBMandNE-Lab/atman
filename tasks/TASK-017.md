# TASK-017 — unify the six hand-rolled RNGs into atman_core::rng::SplitMix64

<!-- kanban-status:start -->
**Status:** review
**Owner:** orchestrator
**Branch:** (pending)
**Started:** 2026-05-29T13:08:21+00:00
**Updated:** 2026-05-29T21:05:26+00:00
<!-- kanban-status:end -->

## Plan

1. Add one reviewed `atman_core::rng` module with `SplitMix64` and `Xoshiro256pp`,
   each exposing `next_u64()` and an UNBIASED `bounded(upper)` (Lemire). Document a
   single central seed rule (`fold_zero_seed`) and one `derive_sub_seed`.
2. Make `ica.rs` re-export `Xoshiro256pp` from `rng` so the whole core decomposition
   family keeps its `crate::ica::Xoshiro256pp` import path.
3. Migrate the three CLI generators (`null`, `bootstrap`, `ratio`) to `SplitMix64`.
4. Migrate the two remaining core generators (`multivariate_t::McRng`,
   `decompose_unmix` sub-seeds) and fix the buggy bounded draws scattered across the
   decomposition family (`gsea`, `ica_null`, `align_bootstrap`, `decompose_unmix`).
5. Preserve the raw u64 stream where call sites only consume it; change only the
   bounded-integer draws.

## Notes

### Single primitive
New `crates/atman-core/src/rng.rs` is the only PRNG home:
- `SplitMix64` — splitmix64 stream; used by CLI `null`, `bootstrap`, `ratio`.
- `Xoshiro256pp` — xoshiro256++ seeded via splitmix; used by the core decomposition
  family and Dunnett-Hsu MC. `ica.rs` now `pub use`s it (no second definition).
- Unbiased bounded draw: `bounded(upper)` via Lemire's nearly-divisionless method
  with rejection — the only sanctioned uniform-integer draw. No `next_u64() % upper`
  remains anywhere.
- Seed rule, applied once centrally in `fold_zero_seed`: a raw `0` seed folds to
  `0x9E3779B97F4A7C15`; every other seed passes through. This subsumes the divergent
  ad-hoc `if seed == 0` checks (null/bootstrap had them; ratio/ica/multivariate_t did
  not).

### Six (and the latent extras) migrated
All listed sites migrated; no hand-rolled `Rng64`/`McRng` structs remain:
- `crates/atman/src/commands/null.rs` — `Rng64` → `SplitMix64`.
- `crates/atman/src/commands/bootstrap.rs` — `Rng64` → `SplitMix64`.
- `crates/atman/src/commands/ratio.rs` — `Rng64` (LCG+splitmix) → `SplitMix64`.
- `crates/atman-core/src/ica.rs` — local `Xoshiro256pp` def removed, re-exported from `rng`.
- `crates/atman-core/src/multivariate_t.rs` — `McRng` → `rng::Xoshiro256pp`.
- `crates/atman-core/src/decompose_unmix.rs` — local `derive_sub_seed`/`bootstrap_sub_seed`
  → central `rng::derive_sub_seed`; `sample_indices_with_replacement` rewritten to use
  `bounded`.
- Also fixed the same latent bug in the rest of the family that shares `Xoshiro256pp`:
  `gsea.rs` (`rng_index`), `ica_null.rs` (`permutation`/`rng_u64`), `align_bootstrap.rs`
  (`bootstrap_indices`). These previously drew bounded integers from
  `rng.next_normal().to_bits()` (float bit-pattern, not the integer stream) with a
  modulo reduction — doubly wrong; now `next_u64`-backed unbiased `bounded`.

### Absolute-output changes (intentional correctness fixes)
Per the task, the unbiased bounded draw changes the integer sequence, so the absolute
outputs of these RNG-driven commands change vs the prior release. Two-run byte-identity
at a fixed seed still holds.
- `atman null` — permutation/sign-flip sequence now unbiased; empirical_p / null
  summaries shift slightly.
- `atman bootstrap protein|module|program` — resample index sequence now unbiased;
  CI bounds / bootstrap means shift slightly.
- `atman ratio` (bootstrap + permutation modes) — generator changed (LCG+splitmix →
  splitmix) and bounded draw unbiased; median-diff CI shifts. Point estimate, medians,
  and analytic p-values are RNG-independent and unchanged.
- `atman decompose unmix --n-boot` — bootstrap CI resampling now uses the integer
  stream with an unbiased draw (was float-bits + modulo); loading/abundance CI bounds
  shift. Point endmembers/abundances unchanged (VCA/FCLS init stream preserved).
- `atman decompose null` / GSEA permutation / `align bootstrap` — permutation/subset
  draws now unbiased; null distributions shift slightly.
- `de --posthoc dunnett-hsu` (multivariate_t MC) — `Xoshiro256pp` MC stream replaces
  the `McRng` stream; MC CDF estimates shift within Monte-Carlo noise (tolerance-bounded
  tests pass).

### Stream preservation
The xoshiro `next_u64`/`next_normal` raw stream is byte-identical to the legacy
`ica.rs` generator: `Xoshiro256pp::new` pre-advances the splitmix init by one
`GOLDEN_GAMMA` to reproduce the old `seed + GAMMA` seeding exactly. This keeps ICA,
NMF, VCA/FCLS point estimates, and `decompose ica`/`nmf` byte-identical (verified by
`nmf::select_k_cophenetic` and the `decompose_*` golden/property tests).

## Verification

- single `atman_core::rng` with unbiased Lemire bounded draw + documented seed rule — done.
- all six call sites migrated; `grep "struct Rng64\|McRng"` over non-test src is empty — done.
- `cargo build --workspace` — clean.
- `cargo test -p atman-core` — all pass (incl. `rng::tests` unbiasedness, determinism,
  zero-seed fold).
- `cargo test --workspace` — 75 test binaries, 0 failures.
- byte-identity / determinism explicitly re-run and green: `determinism_new_methods`,
  `null`, `ratio`, `bootstrap_{protein,module,program}`, `decompose_unmix`,
  `decompose_null`, `align_bootstrap`, `de_posthoc_dunnett`, `e2e_demo_paper_chain`.
- No golden RNG-derived fixture required regeneration: the affected CLI tests assert
  two-run identity + structural/inequality properties (e.g. ratio `sign_stability>0.95`,
  bootstrap counts), and the decomposition golden tests depend only on the preserved
  raw stream. Nothing encodes an externally-published RNG value.

## Hand-off

Reviewable in one module: `crates/atman-core/src/rng.rs`. Follow-ups (out of scope):
`nmf.rs` derives uniforms via `next_normal().abs()` (works, deterministic) and could
move to `next_f64`; the `Xoshiro256pp::new` GAMMA pre-advance is a stream-compat shim
that a future "clean break" release could drop alongside a one-time fixture refresh.

### Codex adversarial-review fix (orchestrator)

Codex flagged (MEDIUM): Xoshiro256pp::new applied fold_zero_seed BEFORE the GOLDEN_GAMMA pre-advance, so for seed==0 the stream became mix(3G..6G) instead of the legacy mix(2G..5G) — silently changing ICA/NMF/VCA/decompose point estimates at seed 0 (the one case the "point estimates unchanged" claim missed). Fix: the xoshiro stream-compat path now uses the RAW seed (`seed + GAMMA`, no zero-fold), reproducing legacy seeding byte-for-byte for EVERY seed including 0 (the splitmix expansion can't yield an all-zero xoshiro state, so the fold is unnecessary here). SplitMix64 (null/bootstrap bounded-draw RNGs) keeps fold_zero_seed, matching their legacy seed==0 special-casing. Replaced the test that encoded the buggy fold (asserted new(0)==new(GAMMA)) with `xoshiro_seed_zero_preserves_legacy_stream` (pins seed-0 state to the legacy mix(2G..5G) expansion and asserts new(0) != new(GAMMA)). atman-core (245) + decompose/ica/nmf/unmix/null/bootstrap determinism suites all green.
