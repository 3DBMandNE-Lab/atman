# Atman Showcase: 4-Week Execution Plan

## Objective

Produce a public-release-quality showcase manuscript around the Dube heat dataset that demonstrates:

1. deterministic end-to-end reproducibility from raw Olink NPX files,
2. biologically coherent acute heat-stress proteomic response,
3. disciplined, sensitivity-aware interpretation of acclimation and phenotype linkage.

## Constraints

1. No expansion of the public CLI surface for one-off analyses.
2. Keep core engine behavior stable and test-backed.
3. Separate confirmed claims from exploratory claims.

## Week 1 (Executed now): Repro + Core DE Baseline

### Deliverables

1. Fresh run of `ingest -> qc -> matrix -> fold-change -> de` on raw NPX in `example_data/dube_heat_2023`.
2. Quantitative baseline memo with comparison-level discovery counts and canonical heat-shock checks.
3. Freeze baseline artifacts in `/tmp/karna_showcase` for reproducible downstream analyses.

### Success Criteria

1. Pipeline runs cleanly from raw inputs.
2. DE output stable and interpretable.
3. Canonical acute heat signal present in both acute contrasts.

## Week 2: Robustness and Sensitivity

### Deliverables

1. Leave-one-subject-out sensitivity on DE direction/significance.
2. Threshold sensitivity on `min_pairs` and `q` cutoffs.
3. Ranking-stability tables for top proteins and pathways.

### Success Criteria

1. Core acute signal remains directionally stable under subject perturbation.
2. Report explicitly lists stable vs unstable findings.

## Week 3: Biology Layer and Narrative Assembly

### Deliverables

1. Curated biological interpretation for acute heat proteins and top restricted-universe pathways.
2. Effect-size-first result panels (not p-value-first).
3. Figure-ready tables for manuscript plotting (volcano inputs, heatmaps, comparison summaries).

### Success Criteria

1. Main claim phrasing is specific and defensible.
2. Exploratory sections are clearly labeled and bounded.

## Week 4: Manuscript Package and Release Artifacts

### Deliverables

1. Showcase manuscript draft with claim tiers: Confirmed / Suggestive / Null.
2. Reproducibility appendix with exact command transcript and environment assumptions.
3. Final QA pass: tests green, docs coherent, no accidental API/tooling changes.

### Success Criteria

1. Someone external can rerun the paper pipeline from raw files.
2. Claims map directly to generated tables/artifacts.
3. Public tool surface remains clean.

