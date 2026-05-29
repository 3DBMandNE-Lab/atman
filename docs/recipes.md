# Atman DE Recipes

Atman's `de` command ships ten test paths. The README Quick Start covers
`paired-t` end-to-end; this file is the cookbook for everything else.
All recipes assume `out/` already contains the canonical TSVs produced
by `atman ingest-matrix` (or a Python adapter from `adapters/`) plus
the appropriate QC, validate, and report steps from the Quick Start.

## Two-group paired and unpaired

```bash
# Paired (already in the Quick Start; repeated here for context)
atman de \
    --input-dir out --output-dir out \
    --test paired-t \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" \
    --min-pairs 5

# Welch two-sample (unequal variance)
atman de \
    --input-dir out --output-dir out_welch \
    --test welch-t \
    --groups "Case-Control" \
    --min-pairs 5
```

## Covariate-adjusted OLS

```bash
atman de \
    --input-dir out --output-dir out_ols \
    --test ols \
    --groups "Case-Control" \
    --design "~ condition + age + sex + batch" \
    --contrast conditionCase \
    --min-pairs 5
```

## Multi-level omnibus F

Reports per-protein F for whether a categorical factor (here `stage`
with 3 levels CN/MCI/AD) has **any** effect, alongside the standard
pairwise contrast. Emits `de_omnibus.tsv` with BH-adjusted q within each
`(comparison, panel)`.

```bash
atman de \
    --input-dir out --output-dir out_ols_omnibus \
    --test ols \
    --groups "Case-Control" \
    --design "~ condition + stage" \
    --contrast conditionCase \
    --omnibus-factor stage \
    --min-pairs 5
```

## Post-hoc pairwise contrasts (Sidak / Tukey / Dunnett)

Single OLS fit, all `choose(k, 2)` pairwise contrasts corrected for
multiple comparisons.

```bash
# Sidak step-down
atman de --post-hoc sidak \
    --input-dir out --output-dir out_sidak \
    --test ols --groups "stage" \
    --design "~ stage + age + sex" \
    --min-pairs 5

# Tukey HSD via from-scratch studentized range
atman de --post-hoc tukey  --input-dir out --output-dir out_tukey ...

# Dunnett against a control level (balanced or Dunnett–Hsu unbalanced,
# auto-detected from observed per-group n)
atman de --post-hoc dunnett --control-level CN \
    --input-dir out --output-dir out_dunnett ...
```

## Repeated-measures mixed model

```bash
atman de \
    --input-dir out --output-dir out_mixed \
    --test mixed \
    --groups "PT2-PT1" \
    --fixed "condition + age + sex" \
    --random "1|subject_id" \
    --min-pairs 5
```

## limma eBayes + TREAT

```bash
# limma-grade eBayes with parametric mean-variance trend + robust prior
atman de \
    --input-dir out --output-dir out_limma \
    --test limma \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" \
    \
    --trend true --robust true \
    --min-pairs 5

# TREAT (minimum-effect test) at log2-FC threshold 0.5
atman de \
    --input-dir out --output-dir out_limma_treat \
    --test limma \
    --groups "PT1-PR1" \
    --lfc-threshold 0.5 \
    --min-pairs 5
```

## DEqMS (peptide-count-weighted trend)

Swaps limma's mean-variance trend covariate for `log(peptide_count + 1)`
via a tricube kernel smoother (Zhu et al. 2020). Requires `peptides.tsv`
for the peptide→protein mapping; proteins absent from the catalog get
count 0.

```bash
atman de \
    --input-dir out --output-dir out_deqms \
    --test limma \
    --peptide-metadata out/peptides.tsv \
    --groups "Case-Control" \
    --trend true \
    --min-pairs 5
```

## msqrob (peptide-level ridge mixed model)

One protein → one model. Peptide-level observations with a random
intercept per peptide and a ridge penalty on non-intercept fixed
effects. Requires both `peptide_measurements.tsv` and `peptides.tsv`.

```bash
atman de \
    --input-dir out --output-dir out_msqrob \
    --test msqrob \
    --peptide-measurements out/peptide_measurements.tsv \
    --peptide-metadata out/peptides.tsv \
    --groups "Case-Control" \
    --ridge-lambda 0.5 \
    --min-peptides 2 \
    --min-pairs 5
```

## Cross-method ensemble consensus

Runs every listed method on the same canonical inputs and, within each
`(comparison, protein)`, assigns a VALIDATED / PROVISIONAL /
INSUFFICIENT grade from method agreement — how many of the fitted
methods individually clear per-method BH-q significance
(`n_significant`) and how consistent their effect sign is
(`n_sign_consistent`). Methods with missing inputs (msqrob without
`--peptide-measurements`, etc.) are auto-skipped, not an error.

The grade is:

- **VALIDATED** — every applied method is individually significant and
  agrees on the effect direction.
- **PROVISIONAL** — a majority of applied methods are significant and
  agree on direction.
- **INSUFFICIENT** — otherwise.

(The exact fractions are tunable via `--ensemble-sign-fraction` and
`--ensemble-provisional-fraction`; the same fraction gates both the
significant-method fraction and the sign-consistency fraction.)

**`ensemble_p` / `ensemble_q` are a heuristic, not a calibrated
p-value.** The `de_ensemble.tsv` columns `ensemble_p` (Stouffer-combined
per-method p) and `ensemble_q` (BH-adjusted across proteins) are emitted
as a convenience for ranking only. Stouffer's method assumes the
combined p-values are independent; `paired-t`, `welch-t`, `ols`,
`mixed`, `limma`, and `msqrob` all share most of their signal on the
same measurement matrix, so these combined p-values are
anti-conservative. They do NOT drive the grade and must not be read as
statistical significance.

**Interpretive caveat.** Ensemble is a within-dataset robustness check,
not independent-study meta-analysis. Read VALIDATED as "the finding
survives method swap on this dataset," not "independently replicated."

```bash
atman de \
    --input-dir out --output-dir out_ensemble \
    --test ensemble \
    --ensemble-methods "paired-t,welch-t,limma,msqrob" \
    --peptide-measurements out/peptide_measurements.tsv \
    --peptide-metadata out/peptides.tsv \
    --groups "PT2-PR2" \
    \
    --min-pairs 5
```

## Protein-level bootstrap uncertainty

```bash
atman bootstrap protein \
    --input-dir out \
    --groups "PT2-PT1" \
    --test paired-t \
    --n 2000 \
    --seed 1 \
    --output out/protein_bootstrap.tsv
```

## Module-level bootstrap uncertainty

```bash
atman bootstrap module \
    --input-dir out \
    --modules-tsv modules.tsv \
    --groups "PT2-PT1" \
    --test paired-t \
    --n 2000 \
    --seed 1 \
    --output out/module_bootstrap.tsv
```

## Permutation / sign-flip null calibration

```bash
atman null \
    --input-dir out \
    --output-dir out/null \
    --groups "PT2-PT1" \
    --test paired-t \
    --n 1000 \
    --seed 1
```

## Detection-aware proteomics

`atman detectability` keeps abundance and detection as separate outcomes.
For each requested comparison it fits a logistic detection model
(`detected ~ condition + covariates`) and, separately, a detected-only
abundance model (`abundance ~ condition + covariates`). Missing or
below-LOD values are not imputed.

```bash
atman detectability \
    --input-dir out \
    --groups "Case-Control" \
    --design "~ age + sex + batch" \
    --min-samples 5 \
    --min-detected 3 \
    --min-detected-per-condition 2 \
    --min-missing 3 \
    --output out/detectability.tsv
```

The abundance layer additionally requires `--min-detected-per-condition`
(default 2) detected samples in *each* condition. This suppresses one-sided
splits — for example a protein detected in 3 cases and 0 controls — where the
condition coefficient carries no real between-group contrast even though the
global detected count clears `--min-detected`. Such rows are tagged
`one_sided_detection` in `skip_reason` and report no abundance effect.

Output columns include detection counts by group, `detection_or`,
`detection_q`, `abundance_effect_detected_only`, `abundance_q`, and a
classification: `detection_shifted`, `abundance_shifted`,
`coupled_shift`, `not_significant`, or `uninformative_sparse`. A companion
`*_summary.tsv` flags designs where the detection layer is structurally
uninformative, for example plex-balanced TMT data where almost every protein
has identical detected/missing counts across conditions.

## TMT plex recovery from absence pattern

`atman recover-plex` reconstructs the TMT plex assignment when ingest dropped
`plate_id`. It exploits the fact that on TMT data a peptide that fails to
quantify is missing for *every* sample in its plex — so plex membership is
latent in the per-sample absence pattern. Sample-sample Jaccard distance on
the binary absence vector clusters samples back into their plexes via
connected-components at a strict threshold (default 0.05).

```bash
atman recover-plex \
    --input-dir out \
    --output out/recovered_plex.tsv \
    --augmented-samples-output out/samples_with_plex.tsv
```

Outputs:

- `recovered_plex.tsv` — per-sample `recovered_plex_id`, cluster size, mean
  within-cluster Jaccard distance.
- `recovered_plex_quality.tsv` — global metrics including
  `mean_within_jaccard_distance`, `mean_between_jaccard_distance`,
  `separation_ratio`, and a `diagnostic` of `plex_coherent`,
  `weak_plex_signal`, `atypical_cluster_size`, or
  `absence_not_plex_coherent`.
- Optional `--augmented-samples-output` writes a copy of `samples.tsv` with
  a new `recovered_plex_id` column. Only emitted when
  `diagnostic == plex_coherent`, to avoid downstream poisoning.

Use the recovered id as a covariate in any downstream design, e.g.
`atman detectability --design "~ recovered_plex_id"` or
`atman de --design "~ recovered_plex_id"`. On non-TMT (or any data where
absence is not plex-categorical) the `diagnostic` will fail honestly and the
augmented samples file is suppressed.

When the cohort has a paired design (e.g. tumor / paired_non_tumor) but the
canonical schema doesn't carry a real patient id, `--infer-pairs A-B` matches
each sample of condition A to its nearest-stem-distance partner of condition
B *within the same recovered plex*. Inference succeeds only when every A has
a unique B partner at distance ≤ `--pair-stem-distance` (default 1); on
success the augmented samples file gains a `patient_id` column and a
companion `*_pairs.tsv` lists the pairs and their stem distances. Use the
recovered patient id as the pairing covariate in DE, e.g.
`atman detectability --design "~ patient_id"`.

## Protein absence topology

`atman absence-topology` is the protein-side dual of `recover-plex`: it
clusters proteins by the set of samples in which they fail to quantify.
Two proteins that go missing in the same samples (e.g. share a peptide-level
fragility, or always drop out of the same plexes) end up in the same
cluster. Proteins with fewer than `--min-absent` (default 5) absences are
reported as the single `complete` group.

```bash
atman absence-topology \
    --input-dir out \
    --output out/absence_topology.tsv \
    --min-absent 5
```

Outputs:

- `absence_topology.tsv` — per-protein `cluster_id`, `cluster_size`,
  `n_absent`, `mean_within_jaccard_distance`.
- `absence_topology_quality.tsv` — global metrics including
  `separation_ratio`, `n_complete`, `n_singletons`, and a `diagnostic` of
  `protein_coherent`, `weak_signal`, or `absence_not_protein_coherent`.

On CPTAC TMT data the dominant pattern is per-protein-unique fragility
(most clusterable proteins land in singletons), but small fragility cliques
do emerge — typically 2-30 proteins that all fail in exactly the same plex.

## Leave-one-pair-out robustness for paired DE

`atman robust-paired` runs a paired t on the full set of pairs and on every
n-1 jackknife replicate (one pair dropped). For each protein it reports
`loso_sign_stability` (fraction of replicates retaining the sign of
`mean_diff`), `loso_p_lt_05_stability` (fraction with `p < 0.05`),
`loso_min_t` / `loso_max_t`, `loso_max_abs_t_delta`, and
`most_influential_pair_id` (the pair whose deletion changes |t| most).

```bash
atman robust-paired \
    --input-dir out \
    --groups "tumor-paired_non_tumor" \
    --paired-by patient_id \
    --output out/robust_paired.tsv
```

Pairing comes from a column in `samples.tsv` (default `patient_id`) — for
example, the column emitted by `atman recover-plex --infer-pairs`. The
companion `*_summary.tsv` reports how many proteins are simultaneously
significant under the full test (q < 0.05) and robust under jackknife
(`sign_stability ≥ 0.999` AND `p_lt_05_stability ≥ 0.95`); the gap exposes
calls that depend on a small number of subjects.

## Module discovery (WGCNA soft-threshold)

Data-driven module discovery from canonical measurements, written as a
`modules.tsv` that drops straight into `atman score modules`. Soft-power
β is auto-picked by scale-free `R²`; the run sidecar logs a similarity
audit (how many feature pairs hit an undefined metric and fell back to
zero).

```bash
atman modules discover \
    --input-dir out \
    --similarity pearson \
    --method wgcna-soft \
    --soft-power auto \
    --r2-target 0.8 \
    --min-module-size 5 \
    --cut-height 0.5 \
    --output-dir out/modules
```

Outputs: `modules_discovered.tsv` (gene → module assignment),
`module_discovery_report.tsv` (size, hub feature, eigengene PC1
variance per module), and `soft_power_diagnostics.tsv` for the
auto-β sweep.

## Module scoring → module-level DE

```bash
# Score modules per sample (mean, median, z-score, or PC1)
atman score modules \
    --input-dir out \
    --modules-tsv modules.tsv \
    --method mean \
    --output out/module_scores.tsv \
    --canonical-output-dir out/module_score_canonical

# Module-level DE uses the canonical TSVs written above
atman module-de \
    --input-dir out/module_score_canonical \
    --groups "PT2-PT1" \
    --test paired-t \
    --output out/module_de.tsv
```

## ORA enrichment

```bash
atman enrich ora \
    --de-results out/de_results.tsv \
    --gene-sets gene_sets.tsv \
    --comparison "PT2-PT1" \
    --output out/ora.tsv
```

## GSEA (pre-ranked gene-set enrichment)

`atman enrich gsea` runs the weighted Kolmogorov-Smirnov enrichment of
Subramanian et al. (2005) using the position-mask permutation
formulation of `fgsea::fgseaSimple`. Operates on the full ranked DE
list — no q-cutoff, no foreground/background split — so coordinated
subtle shifts are captured. Validated against `fgsea::fgseaSimple` on
a planted fixture: ES agreement within `1e-12`. See
[reference.md](reference.md#fgsea-parity-note) for the full parity
statement.

```bash
atman enrich gsea \
    --de-results out/de_results.tsv \
    --gene-sets gene_sets.tsv \
    --comparison "PT2-PT1" \
    --output out/gsea.tsv \
    --rank-by signed_log10_p \
    --n-permutations 1000 \
    --seed 42
```

Output columns: `set_name`, `set_size`, `es`, `nes`, `p_value`,
`bh_q`, `leading_edge`. Sorted ascending by `bh_q`, then descending
by `|nes|`. Alternative ranking statistics: `--rank-by t` or
`--rank-by log2_fc`. Set-size filtering: `--min-set-size` (default
5), `--max-set-size` (default 500).

## Per-sample signature scoring (singscore)

`atman score signatures` produces a per-sample score for each gene
set in a `gene_sets.tsv` library — useful for downstream phenotype
correlation, subtype stratification, and pathway-level visualization.
v1 implements singscore (Foroutan et al. 2018): per-sample
average-rank-based score, centered to `[-0.5, 0.5]`. Deterministic,
no permutations. Validated against `singscore::simpleScore` on a
planted fixture: TotalScore agreement within `1e-12`. See
[reference.md](reference.md#singscore-parity-note) for the parity
statement.

```bash
atman score signatures \
    --input-dir out \
    --gene-sets gene_sets.tsv \
    --output out/signature_scores.tsv \
    --method singscore \
    --min-set-size 3
```

Output columns: `sample_id`, `subject_id`, `condition`, `set_name`,
`method`, `score`, `n_genes_declared`, `n_genes_observed_in_sample`,
`n_proteins_in_sample`. Sorted by `(set_name, sample_id)`. A
signature whose observed-in-sample count falls below `--min-set-size`
emits `score = NaN` so per-sample coverage stays auditable.

## Cross-cohort differential coexpression

`atman network differential` computes per-cohort signed correlation
matrices on the shared feature universe and summarizes how each
protein-pair edge varies across cohorts. Three runtime-selectable
output modes:

```bash
# Mode 1: per-edge × per-cohort-pair Fisher-z test (DGCA-style).
atman network differential \
    --inputs BRCA=runs/BRCA,COAD=runs/COAD,GBM=runs/GBM,HCC=runs/HCC,HNSCC=runs/HNSCC,LUAD=runs/LUAD \
    --output runs/diff/edge_pairwise.tsv \
    --mode edge-pairwise \
    --min-overlap 5 \
    --top-rows 10000

# Mode 2: per-edge cross-cohort summary (one row per edge).
atman network differential \
    --inputs BRCA=runs/BRCA,COAD=runs/COAD,GBM=runs/GBM,HCC=runs/HCC,HNSCC=runs/HNSCC,LUAD=runs/LUAD \
    --output runs/diff/edge_summary.tsv \
    --mode edge-summary \
    --min-overlap 5

# Mode 3: per-module within-module connectivity rewiring.
atman network differential \
    --inputs BRCA=runs/BRCA,COAD=runs/COAD,GBM=runs/GBM,HCC=runs/HCC,HNSCC=runs/HNSCC,LUAD=runs/LUAD \
    --output runs/diff/module_rewiring.tsv \
    --mode module \
    --gene-sets hallmarks.tsv
```

Edge-summary output columns: `feature_a`, `feature_b`, `n_cohorts`,
`mean_corr`, `sd_corr`, `min_abs_corr`, `max_abs_corr`, `range_corr`,
`n_sign_flips`, `conservation_score` (min |r| across cohorts —
high when every cohort agrees), `divergence_score` (sd / (|mean| +
1e-3) — high when cohorts disagree), then one `<COHORT>_corr` column
per cohort. Sorted by `divergence_score` descending.

Module output columns: `set_name`, `set_size_declared`,
`set_size_observed`, `mean_connectivity`, `sd_connectivity`,
`range_connectivity`, `rewiring_score` (= sd of per-cohort within-
module mean |r|), then one `<COHORT>_connectivity` column per cohort.
Sorted by `rewiring_score` descending.

## Cross-cohort meta-analysis

```bash
atman meta \
    --inputs cohort1/de_results.tsv,cohort2/de_results.tsv \
    --output meta.tsv
```

## Plan-manifest reproducibility (`atman run`)

Declarative alternative to a shell script. Writes
`plan_manifest.tsv` with SHA-256 of every declared input and output
per stage, atman version, OS/arch, and per-stage wall-clock. Refuses
to overwrite a manifest whose `plan_commit` hasn't been bumped when
the plan content changes (unless `--allow-drift` is passed).

```yaml
# plan.yaml
plan_commit: "2026-04-20.v1"
stages:
  - name: ingest
    inputs: [example_data/dube_heat_2023/*.csv]
    outputs: [out/samples.tsv, out/proteins.tsv, out/measurements.tsv]
    cmd: >
      python3 adapters/generic/olink_explore_to_atman.py
      --output-dir out example_data/dube_heat_2023/*.csv
  - name: de
    inputs: [out/measurements.tsv, out/samples.tsv]
    outputs: [out/de_results.tsv, out/de_report.tsv]
    cmd: >
      atman de --input-dir out --output-dir out
      --test paired-t
      --groups "PT1-PR1,PR2-PR1" --min-pairs 5
```

```bash
atman run --plan plan.yaml --manifest out/plan_manifest.tsv
```

## Decomposition: discovering protein programs

Standard PCA tells you which samples are different. Atman's decomposition
commands ask a different question: **what underlying protein programs
explain the variation, and how much of each program is active in each
subject?**

This matters in proteomics because circulating protein panels (plasma,
CSF, synovial fluid) are mixtures of secreted programs — inflammatory
signaling, tissue leakage, acute-phase response, complement cascades.
Decomposition recovers those programs as interpretable components with
per-subject activation scores, so you can ask "which program drives the
case-control difference?" rather than "which individual proteins change?"

### ICA: statistically independent programs

FastICA finds maximally non-Gaussian (statistically independent)
components. Run across multiple random seeds to assess stability —
components that appear in most seeds are robust; components that
fragment across seeds are noise.

```bash
# Multi-seed ICA: 50 seeds, keep components explaining 80% of variance
atman decompose ica \
    --input-dir out \
    --k-selection cumulative-variance=0.80 \
    --n-seeds 50 \
    --seed 42 \
    --transform clr \
    --output-loadings out/loadings.tsv \
    --output-activations out/activations.tsv \
    --output-stability out/stability.tsv

# Null calibration: how stable are components under permutation?
# --k requires the concrete integer chosen by ICA (read from loadings column count)
atman decompose null \
    --input-dir out \
    --k 5 \
    --n-perm 200 \
    --seed 42 \
    --output out/null.tsv
```

Outputs: `loadings.tsv` (protein weights per component),
`activations.tsv` (subject activation per component),
`stability.tsv` (cross-seed reproducibility per component).

#### Filtering interpretable programs

`programs filter` flags ICA programs that fail any of three checks:
weak top annotation p-value, diffuse loading (no single-feature peak),
or contamination signature (the top loadings match a regex pattern of
known contaminants). The default contamination pattern targets keratin
(`(?i)^KRT|KERATIN`) — appropriate for plasma/serum/CSF cohorts where
keratin is a skin-shedding contaminant. Override the pattern for other
contexts (e.g. `(?i)^HB[AB]` for hemolysis, `(?i)^IG[HKL]` for
immunoglobulin carryover) or for tissues where keratin is real biology.

```bash
atman programs filter \
    --loadings out/loadings.tsv \
    --annotations out/program_annotations.tsv \
    --min-annotation-pvalue 0.05 \
    --min-top-loading 0.5 \
    --max-contamination-fraction 0.2 \
    --contamination-pattern '(?i)^KRT|KERATIN' \
    --top-n 20 \
    --output out/programs_interpretable.tsv
```

Output column `fail_reasons` is `none` for interpretable programs, or
a `;`-joined list of `missing_annotation`, `annotation_p_value`,
`diffuse_loading`, `contamination_signature`. The pattern used and
fraction threshold are recorded in the run sidecar.

### VCA+FCLS unmixing: compositional endmembers

When the biological question is "what pure sources contribute to this
mixture?" — e.g. plasma as a mixture of liver secretome, immune
activation, and tissue leakage — unmixing fits better than ICA.
VCA finds geometric endmembers (extreme compositions), FCLS
estimates fractional abundances that sum to 1 per subject.

```bash
atman decompose unmix \
    --input-dir out \
    --k auto \
    --seed 42 \
    --n-boot 1000 \
    --annotate-markers gene_sets.tsv \
    --output-dir out/unmix
```

Outputs: `endmembers.tsv` (rank-ordered proteins per endmember),
`abundances.tsv` (per-subject fractions, rows sum to 1),
`unmix_diagnostics.tsv`, `marker_enrichment.tsv`.

### Variance partition: attributing variation to covariates

Before running DE, check how much protein-level variance is explained
by your factor of interest vs batch, age, sex, or technical covariates.
Uses a mixed model per protein with Type III F-tests.

```bash
atman decompose variance \
    --activations out/activations.tsv \
    --samples out/samples.tsv \
    --factors "diagnosis + sex + age + (1|batch)" \
    --output out/variance.tsv
```

### NMF: biologically-motivated factorization

Non-negative matrix factorization identifies a fixed number of non-negative
components (latent protein programs) whose product reconstructs the measured
abundances. Unlike ICA (which seeks statistical independence), NMF produces
sparse, interpretable loadings and is appropriate when:

- You expect components with **clear membership** (feature either loads or doesn't)
- Biological signals are **non-Gaussian** in a specific way (sparse positive sources)
- You want a fixed, interpretable *k* rather than variance-determined dimension

Run multi-seed NMF with stability assessment:

```bash
atman decompose nmf \
    --input-dir out \
    --k-selection cophenetic-knee \
    --k-min 2 --k-max 8 \
    --beta-loss frobenius \
    --n-seeds 10 \
    --min-stable-seed-fraction 0.9 \
    --output-loadings out/nmf_loadings.tsv \
    --output-activations out/nmf_activations.tsv \
    --output-stability out/nmf_stability.tsv \
    --output-k-sweep out/nmf_k_sweep.tsv
```

Outputs:
- `nmf_loadings.tsv` (protein weights per component)
- `nmf_activations.tsv` (subject activation per component)
- `nmf_stability.tsv` (cross-seed reproducibility per component)
- `nmf_k_sweep.tsv` (cophenetic and RSS metrics across k values)

### Admixture-adjusted DE on bulk tumor proteomics

Workflow combining NMF decomposition with admixture-adjusted differential
abundance testing. Use this when bulk tumor samples contain unknown mixtures
of cell types and you want to control for compositional variation without
explicitly deconvolving. The four-step chain is:

**Step 1. Decompose with NMF to recover compositional programs:**

```bash
atman decompose nmf --input-dir runs/MY_COHORT \
    --k-selection cophenetic-knee --k-min 2 --k-max 8 \
    --beta-loss frobenius --n-seeds 10 \
    --min-stable-seed-fraction 0.9 \
    --output-loadings runs/MY_COHORT/nmf_loadings.tsv \
    --output-activations runs/MY_COHORT/nmf_activations.tsv \
    --output-stability runs/MY_COHORT/nmf_stability.tsv \
    --output-k-sweep runs/MY_COHORT/nmf_k_sweep.tsv
```

**Step 2. Compose NMF activations into a per-sample covariate matrix and run
adjusted DE:**

```bash
atman de --test limma --input-dir runs/MY_COHORT \
    --groups subtype_a-subtype_b \
    --adjust-for runs/MY_COHORT/nmf_activations.tsv \
    --output-dir runs/MY_COHORT/de_adjusted/
```

The `--adjust-for` TSV is auto-detected for format: wide (sample_id × covariate
columns) or long (sample_id, program, activation). NMF activation outputs use
long format and are consumed directly. Atman also accepts wide format from
custom external covariates.

**Step 3. Pre-ranked GSEA on the adjusted results:**

```bash
atman enrich gsea --de-results runs/MY_COHORT/de_adjusted/de_results.tsv \
    --gene-sets lib/hallmarks.tsv \
    --output runs/MY_COHORT/gsea_hallmarks.tsv
```

**Step 4. (Optional) Null calibration to sanity-check parametric p-values:**

```bash
atman null --test limma --input-dir runs/MY_COHORT \
    --groups subtype_a-subtype_b \
    --n 1000 --output-dir runs/MY_COHORT/null/
```

**Auditability:** Each step's `*.run.json` sidecar records the SHA-256 of the
previous step's output as input. A reader can reproduce any figure number by
replaying the logged commands in sequence. The sidecar `--adjust-for` input path
and its source hash are recorded for full traceability.
