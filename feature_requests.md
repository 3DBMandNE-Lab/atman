# Atman feature requests

Inbox for cross-instance feature proposals based on project needs.
Shipped work is documented in `CHANGELOG.md`; this file tracks only
what has been requested but is **not yet shipped**.

Add new proposals as sections below. Remove entries once they land
(the commit message + CHANGELOG are the permanent record).

---

## Open

### Hierarchical / nested alignment

Site → cohort → cross-cohort alignment for consortium data where a
single cohort has multiple acquisition sites. Each level of the
hierarchy carries its own bootstrap-based confidence interval
(building on the existing `atman align bootstrap` surface).

**Designed, not scheduled.** Spec at
`docs/superpowers/specs/2026-09-06-hierarchical-alignment-design.md`.
Requested by the CSF CrossDisease session as a revision tool (a
reviewer asking whether Axis 1 replicates within Bader's four
collection sites); not needed for that paper's submission, and no
other manuscript session needs it. Build it when a reviewer asks.

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

---

## Demand survey, 2026-09-06

Every manuscript session on this machine was asked which of these it
needs. Meningioma (Text/v4) and CAR-T do not use atman at all. The GBM
methods paper needs none of them — its run tree is complete at PIN2 and
replays from a snapshotted binary. The CSF CrossDisease paper needs none
for submission and rules out three permanently:

- **Pooled-QC drift correction** has no possible input: no cohort in that
  study records injection order, acquisition date, or pooled-QC samples,
  and the SIH `plate_id` / `panel_lot` columns are empty. `ingest_order`
  is row order at ingest, not acquisition order.
- **Compositional effect size for DE** conflicts with a stated Methods
  choice (log2 + median normalization, Cohen d on that scale). CLR
  coefficients would describe a different pipeline.
- **PMF**, **MCR-ALS** and **longitudinal tensor decomposition** would each
  be a new sensitivity analysis against Methods that commit to FastICA +
  ICASSO; the only longitudinal arm (12 subjects × 3 timepoints) is
  already covered by the existing one-way ICC.

They stay open because a future project may want them, not because
anything is waiting on them. Re-ask before building.
