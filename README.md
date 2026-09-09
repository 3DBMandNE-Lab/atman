# Atman

[![CI](https://github.com/3DBMandNE-Lab/atman/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/3DBMandNE-Lab/atman/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/3DBMandNE-Lab/atman)](https://github.com/3DBMandNE-Lab/atman/releases/latest)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](./LICENSE-MIT)

Atman is a command-line toolkit for proteomics analysis. It takes a
protein matrix from Olink Explore, SomaScan, MaxQuant, DIA-NN or
Spectronaut and runs the same commands on all of them. Every command
writes a sidecar that records the binary, the inputs and the outputs by
SHA-256. Any result can be traced to its source and rerun.

## Install

Atman is a single binary built from source with Rust 1.94 or newer.
Get Rust from [rustup.rs](https://rustup.rs) if you do not have it.

```bash
git clone https://github.com/3DBMandNE-Lab/atman.git
cd atman
cargo install --path crates/atman
atman --version
```

Or build a container:

```bash
docker build --build-arg ATMAN_GIT_SHA=$(git rev-parse HEAD) -t atman:1.2.1 .
docker run --rm atman:1.2.1 --help
```

**macOS is the supported operating system.** Atman is built, validated
and released there. Linux builds and passes the test suite in CI, and
the container is a Linux build. On Linux the decompositions use a
portable kernel instead of Accelerate. They run slower, and the
low-order digits differ from a macOS run. Neither changes a result. See
"What determinism atman guarantees" in [docs/reference.md](docs/reference.md).

## Quick start

Atman needs two files: a protein matrix and a sample metadata table.
The metadata needs a sample ID, a subject ID and a condition per sample.

```text
sample_id	subject_id	condition
DIA_C1	DIA_C1	Control
DIA_C2	DIA_C2	Control
DIA_K1	DIA_K1	Case
DIA_K2	DIA_K2	Case
```

Three commands take a DIA-NN protein-group matrix to differential
abundance results:

```bash
# 1. Ingest: matrix + metadata -> canonical TSVs (measurements, samples, proteins)
atman ingest-matrix \
    --matrix protein_groups.tsv --samples sample_metadata.tsv \
    --orientation proteins-rows --log2-transform \
    --platform diann_report --abundance-unit log2_diann_pg_quantity \
    --assay-id-col Protein.Group --gene-col Genes --uniprot-col Protein.Ids \
    --subject-id-col subject_id --condition-col condition --panel diann \
    --output-dir out

# 2. Validate: schema, keys, and sample support for the comparison
atman validate --input-dir out --groups Case-Control --min-pairs 2

# 3. Differential abundance: Welch t-test, Case versus Control
atman de --input-dir out --output-dir out --test welch-t --groups Case-Control --min-pairs 2
```

This writes `de_results.tsv`, `de_report.tsv` and the sidecar
`de_results.tsv.run.json`, shown here abridged:

```json
{
  "command": "de",
  "atman_version": "1.2.0",
  "atman_git_sha": "bee7653b22838129142de89c7cb08caca08372e9",
  "inputs_sha256": {
    "measurements.tsv": "sha256:6bb8e5b8…",
    "proteins.tsv": "sha256:c7913e94…",
    "samples.tsv": "sha256:c64f0458…"
  },
  "output_files": {
    "out/de_report.tsv": "sha256:e8fc981f…",
    "out/de_results.tsv": "sha256:5e2cac94…"
  },
  "reinvoke": "atman de --groups Case-Control --input-dir out --min-pairs 2 --output-dir out --test welch-t …",
  "build_env": {
    "rustc_version": "rustc 1.94.1 (e408947bf 2026-03-25) (Homebrew)",
    "target_triple": "aarch64-apple-darwin",
    "profile": "release"
  }
}
```

The bundled example runs the same steps with
`bash adapters/examples/diann/run.sh`. MaxQuant, Spectronaut and SomaScan
matrices use the same command with their own `--platform` tag and column
names. See [adapters/examples](adapters/examples) and "Input Support" in
the [reference](docs/reference.md#input-support).

For covariate-adjusted DE, use `--test limma --design "~ condition + age + sex"`.
Extra design terms are read from the metadata table. The
[recipes](docs/recipes.md) cover every test: paired, Welch, OLS, mixed,
limma, DEqMS, msqrob and ensemble.

### Olink Explore

Olink NPX exports are long-format. A Python adapter converts them, and
it needs Python 3 and nothing else:

```bash
python3 adapters/generic/olink_explore_to_atman.py \
    --output-dir out \
    example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv \
    example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv

atman qc --input-dir out --output-dir out
atman validate --input-dir out --groups "PT1-PR1,PR2-PR1" --min-pairs 5
atman de --input-dir out --output-dir out --test paired-t --groups "PT1-PR1,PR2-PR1" --min-pairs 5
```

The [tutorial](docs/tutorial.md) walks through this dataset end to end.

## Commands

`atman <command> --help` documents every flag. The
[reference](docs/reference.md#commands) describes every command and its
outputs. The stages, in the order a study runs them:

**Ingest and QC**

| Command | What it does |
|---|---|
| `ingest-matrix` | wide protein matrix + metadata to canonical TSVs |
| `validate` | check schema, keys, and sample support for a comparison |
| `qc` | apply QC masking rules |
| `report qc` | summarize QC, missingness, and condition support |

Also `matrix` (per-panel wide pivot) and `fold-change` (per-panel log2 fold-change tables).

**Differential abundance**

| Command | What it does |
|---|---|
| `de` | paired, Welch, OLS, mixed, limma, msqrob or ensemble DE |
| `bootstrap protein` | subject-level bootstrap intervals for protein effects |
| `null` | permutation and sign-flip null calibration of DE effects |
| `meta` | combine DE results across cohorts |
| `robust-paired` | leave-one-subject-out sign stability for paired DE |

Also `detectability`, `asymmetry`, `robustness`, `ratio`, `residuals` and `absence-topology`.

**Modules, signatures and enrichment**

| Command | What it does |
|---|---|
| `modules discover` | WGCNA-style module discovery |
| `module-de` | aggregate proteins into modules and test at module level |
| `score modules`, `score signatures`, `score weighted` | per-sample module, singscore and signed-weight scores |
| `enrich ora`, `enrich gsea` | over-representation and pre-ranked enrichment |
| `bootstrap module` | subject-level bootstrap intervals for module effects |

Also `module-trajectory` and `enrich gprofiler`.

**Decomposition and programs**

| Command | What it does |
|---|---|
| `decompose ica`, `decompose nmf` | multi-seed FastICA and NMF with stability reporting |
| `decompose null` | permutation-null calibration of archetype stability |
| `decompose variance` | mixed-model variance partition per archetype |
| `decompose unmix` | VCA + FCLS compartmental unmixing |
| `align programs`, `align bootstrap`, `align project` | cross-cohort alignment, its bootstrap, and projection of a new cohort |
| `programs filter` | flag ICA programs by annotation, loading and contamination |

Also `decompose counterfactual`, `bootstrap program`, `coupling` and `bench decompose`.

**Cross-cohort scores and harmonisation**

| Command | What it does |
|---|---|
| `axes build` | representative archetypes to axes, orthogonalised within cohort |
| `axes contrast`, `axes groups`, `axes anchor` | disease modulation, group summaries, clinical anchoring |
| `axes displacement`, `axes loco`, `axes icc`, `axes tree` | displacement vectors, cohort leave-out, trait stability, contrast tree |
| `harmonize fit`, `harmonize apply` | fit a cross-cohort model, then apply it to a held-out cohort |
| `concordance` | Spearman, sign and Jaccard agreement of effect tables |
| `network influence`, `network differential` | hub scoring and differential coexpression |

Also `within-cohort-rank`, `scale absolute` and `recover-plex`.

**Pipelines**

`atman run` executes a plan file (YAML or JSON) and emits a hash manifest
of every stage. Every command writes a `.run.json` sidecar next to its
primary output.

## Validation

Each statistical path is checked against an R reference implementation
on a published or planted fixture. These checks run in every CI build.

| Atman | R reference | Agreement |
|---|---|---|
| `de --test limma --peptide-metadata` | `DEqMS::spectraCounteBayes` | median Δ 0.000 log₂ on CPTAC UPS1 |
| `de --test msqrob` | `msqrob2::msqrob` | 100% sign agreement on CPTAC UPS1 |
| `de --test limma` | `limma::eBayes` | Δ < 1e-4 on t, p, df and s2_post |
| `de --post-hoc sidak` | `lm()` + `pairwise.t.test` | Δ < 1e-6 |
| `decompose variance --omnibus-factor` | `car::Anova(type = 3)` | Δ p < 1e-6 |
| `enrich gsea` | `fgsea::fgseaSimple` | ES Δ < 1e-12 |
| `score signatures` | `singscore::simpleScore` | Δ < 1e-12 |
| Olink NPX reproduction | Dube et al. 2023 tables | byte-exact filtered NPX |

Fixtures, reference scripts and the estimator notes are under
"Validation Details" in [docs/reference.md](docs/reference.md#validation-details).

## Documentation

- [Tutorial](docs/tutorial.md): the Dube Olink dataset from raw export to results
- [Recipes](docs/recipes.md): every DE test with worked examples, plus modules, enrichment, decomposition and `atman run`
- [Reference](docs/reference.md): commands, data model, sidecar schema, input support, validation, non-goals
- [Adapters](adapters/README.md): bringing other formats to the canonical TSVs
- [SECURITY.md](SECURITY.md): dependency policy, SBOM and the `unsafe` sites

## Citing

Cite the release you used. [CITATION.cff](CITATION.cff) carries the
metadata, and each release page records the commit and a source digest.

The bundled reference dataset:

Gagnon D, Barry H, Barhdadi A, Oussaid E, Mongrain I, Lemieux Perreault LP,
Dubé MP. *A dataset of proteomic changes during human heat stress and heat
acclimation.* Scientific Data (2023).
[10.1038/s41597-023-02809-5](https://doi.org/10.1038/s41597-023-02809-5).

## License

MIT or Apache-2.0, at your option. See [LICENSE-MIT](./LICENSE-MIT) and
[LICENSE-APACHE](./LICENSE-APACHE).

Unless you state otherwise, any contribution you intentionally submit for
inclusion in Atman is dual licensed as above. No additional terms or
conditions apply.

### Data

The code licence does not cover the bundled data. The files under
`example_data/dube_heat_2023/` are the five items of Figshare project
163291, deposited by the authors of the dataset under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Attribute them
to Gagnon D et al., *Scientific Data* (2023),
[10.1038/s41597-023-02809-5](https://doi.org/10.1038/s41597-023-02809-5).
The `cptac_*` test fixtures derive from the NCI CPTAC program's public
data releases and are used for validation only.
