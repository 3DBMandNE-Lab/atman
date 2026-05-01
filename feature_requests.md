# Atman feature requests

Inbox for cross-instance feature proposals based on project needs.
Shipped work is documented in `CHANGELOG.md`; this file tracks only
what has been requested but is **not yet shipped**.

Add new proposals as sections below. Remove entries once they land
(the commit message + CHANGELOG are the permanent record).

---

## Open

### Bounded-heap truncation for `network differential --mode edge-pairwise`

`atman_core::network_differential::edge_pairwise_differential` builds
the full per-edge × per-cohort-pair `Vec<EdgePairwiseRow>` in memory
before the caller in `crates/atman/src/commands/network.rs` applies
`--top-rows` truncation. On the CPTAC pan-cancer paper workload (6701
shared features × 15 cohort pairs ≈ 336 million rows pre-truncation
at ~150 bytes per row) this fails to complete on a 128GB machine
regardless of `--top-rows` value because the algorithm allocates
before the cap. `--mode edge-summary` is unaffected (its caller's
`top_rows` truncation is reached after a per-edge bounded-state
aggregation, not after a full row-vector materialization).

Resolution: push truncation into the algorithm via a bounded
`BinaryHeap` of size `top_rows` keyed by `|z_diff|` descending
(mirror the pattern used in `figures/fig5_pancancer_diffcoex/render.py`
in the atman-paper repo for the top-K most-divergent edges over a
22.4M-row stream). Existing CLI surface and output schema unchanged.

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
