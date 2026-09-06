# Atman feature requests

Inbox for cross-instance feature proposals based on project needs.
Shipped work is documented in `CHANGELOG.md`; this file tracks only
what has been requested but is **not yet shipped**.

Add new proposals as sections below. Remove entries once they land
(the commit message + CHANGELOG are the permanent record).

---

## Open

### Parallel `align bootstrap` iterations (deterministic rayon map)

Parked from the GBM manuscript session (2026-09-06): on the real
six-cohort run, `--n-boot 200` with `--decomposition nmf` takes 9+
hours single-threaded (200 iterations × 6 per-resample NMF fits,
plus the pooled-subject jackknife for BCa acceleration).

The iteration loop in `crates/atman-core/src/align_bootstrap.rs`
(`for iter in 0..params.n_boot`) is already independent per
iteration: each iteration builds its own `Xoshiro256pp` from
`derive_sub_seed(params.seed, iter)` and derives the per-cohort fit
seeds from that sub-seed, so no RNG state crosses iterations. The
only shared state is the per-PE-archetype `BootstrapAcc` (sums plus
the `n_cohorts` series). The jackknife loop in
`jackknife_n_cohorts` (per cohort × per dropped subject) has the
same shape.

Resolution: add `rayon` as a workspace dependency, replace both
loops with a parallel map over the index that returns each
iteration's per-archetype match result, then fold into the
accumulators sequentially in index order. Byte-identical output to
the serial path (a determinism test alongside
`crates/atman/tests/determinism_rng_commands.rs` should assert
this); speedup is roughly core-count×. Expose `--threads` (default:
all cores) so replay manifests can record it. Optionally split
`--nmf-tol` into point-estimate and per-resample tolerances — the
bootstrap aggregates over `n_boot` fits, so a looser per-resample
tol is defensible — but that changes numerics and needs its own
justification.

Not for the current paper's run tree: any binary change re-pins and
breaks the single-commit replay story (everything stays at PIN2).
Pays off for submission-time public replay and the CSF cross-cohort
work.

### Hierarchical / nested alignment

Site → cohort → cross-cohort alignment for consortium data where a
single cohort has multiple acquisition sites. Each level of the
hierarchy carries its own bootstrap-based confidence interval
(building on the existing `atman align bootstrap` surface).

### PMF (positive matrix factorization) with per-cell uncertainty weighting

Port from atmospheric source apportionment (Paatero & Tapper 1994).
Sibling to `atman decompose unmix`: where VCA+FCLS finds geometric
endmembers, PMF finds non-negative factors while natively weighting
each cell by its reciprocal analytical uncertainty. Useful for DIA-MS
inputs that ship per-protein-per-sample CV or detection-limit flags.
Implementable as `atman decompose pmf --uncertainty measurements_cv.tsv`.

### MCR-ALS (multivariate curve resolution — alternating least squares)

From chemometrics. Sibling to `atman decompose unmix`: NMF with
chemistry-aware constraints (non-negativity, unimodality, closure)
applied to the protein-by-sample matrix. Implementable as
`atman decompose mcr`.

### Pooled-QC drift correction

Port from untargeted metabolomics. Standard practice there; near-
absent in proteomics. Fit a per-protein LOESS on injection order
using interleaved pooled QC samples, correct analytical drift before
DE or decomposition. Implementable as `atman qc drift
--injection-order injection_order.tsv --pooled-qc qc_samples.tsv`.

### Longitudinal tensor decomposition

Samples × proteins × timepoints. Relevant for any consortium with
longitudinal arms. Random-subject variance becomes the signal, not
the nuisance — sits naturally on top of the existing
`decompose variance` surface.

### Compositional effect size for DE

A `--effect-size-scale clr` flag on `atman de` that reports CLR-space
coefficients alongside the standard log-FC.

### CSF plasma-endmember sanity check (deferred test)

Real-data sanity check on `sih/qc_measurements.tsv` for
`atman decompose unmix`. Requires the live `atman_inputs_albnorm/`
dataset, which lives downstream of atman's in-repo test surface.
Tracked here so the next in-repo pass picks it up when the dataset
is available.
