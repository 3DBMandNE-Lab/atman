# Atman feature requests

Inbox for cross-instance feature proposals based on project needs.
Shipped work is documented in `CHANGELOG.md`; this file tracks only
what has been requested but is **not yet shipped**.

Add new proposals as sections below. Remove entries once they land
(the commit message + CHANGELOG are the permanent record).

---

## Open

### Per-sample signature scoring (singscore)

`atman score signatures` — per-sample signature scores from a
gene-set library, sibling to `atman score modules`. Same canonical
abundance input, same `gene_sets.tsv` schema as `enrich ora` /
`enrich gsea`. Output is a samples × signatures TSV. v1 implements
singscore (Foroutan 2018): per-sample protein ranking, normalized
mean-rank score per signature; deterministic, no permutations,
handles missing values by per-sample exclusion. `--method` reserved
for future ssGSEA / GSVA implementations. Validated against
`singscore::simpleScore` on a planted fixture.

### Missingness-aware ICA

Extend the compositional-transform surface on `atman decompose ica`
with a joint abundance + detection likelihood at decomposition time
— same philosophy as `atman de --test msqrob` or proDA but at the
component level. Lets archetypes be defined partly by "which proteins
were detectable in this sample," which is real biology for MNAR-heavy
proteomics.

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

### Counterfactual archetype simulation

`atman decompose counterfactual --set-archetype A0004 --to 0 --input
activations.tsv` reconstructs a predicted abundance matrix with one
archetype zeroed out. Interpretive aid.

### Compositional effect size for DE

A `--effect-size-scale clr` flag on `atman de` that reports CLR-space
coefficients alongside the standard log-FC.

### CSF plasma-endmember sanity check (deferred test)

Real-data sanity check on `sih/qc_measurements.tsv` for
`atman decompose unmix`. Requires the live `atman_inputs_albnorm/`
dataset, which lives downstream of atman's in-repo test surface.
Tracked here so the next in-repo pass picks it up when the dataset
is available.
