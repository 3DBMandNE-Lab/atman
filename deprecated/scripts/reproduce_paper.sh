#!/usr/bin/env bash
#
# Reproduce the Atman Bioinformatics paper pipeline end-to-end on the
# Dube et al. 2023 heat-acclimation dataset. Produces every main-text result
# under ./out/ in a single invocation.
#
# Not covered by this script:
#   - LOO robustness (needs n+1 DE runs; see atman robustness --help)
#   - module-trajectory (requires a user-supplied modules.tsv)
#   - second-dataset supplementary figure
#
# Requirements: Rust toolchain (cargo), Python 3.11+, the sibling `genesets`
# CLI on PATH (or set KARNA_GENESETS), and the Dube example data present at
# example_data/dube_heat_2023/.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

OUT="${OUT:-out}"
DATA="example_data/dube_heat_2023"

command -v cargo >/dev/null || { echo "cargo not found on PATH"; exit 1; }
command -v python3 >/dev/null || { echo "python3 not found on PATH"; exit 1; }
[[ -d "$DATA" ]] || { echo "Dube example data missing at $DATA"; exit 1; }

echo "[1/8] build atman"
cargo build --workspace --release
export PATH="$REPO_ROOT/target/release:$PATH"

echo "[2/8] ingest"
atman ingest \
    --platform olink-explore-ngs --parser dube \
    --output-dir "$OUT" \
    "$DATA"/20212016_Dube_NPX_2021-11-30.csv \
    "$DATA"/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv

echo "[3/8] qc (Dube rule)"
atman qc --input-dir "$OUT" --output-dir "$OUT" --rule dube

echo "[4/8] matrix + fold-change"
atman matrix --input-dir "$OUT" --output-dir "$OUT" \
    --format dube-wide --split-by panel
atman fold-change --input-dir "$OUT" --output-dir "$OUT" \
    --groups "PT2-PT1,PR2-PR1,PT2-PR2,PT1-PR1"

echo "[5/8] differential abundance — paired-t + moderated"
atman de --input-dir "$OUT" --output-dir "$OUT" \
    --test paired-t --paired-by participant \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" --min-pairs 5
atman de --input-dir "$OUT" --output-dir "${OUT}_mod" \
    --test moderated --moderation-prior-df 4 --paired-by participant \
    --groups "PT1-PR1,PR2-PR1,PT2-PT1,PT2-PR2" --min-pairs 5

echo "[6/8] asymmetry"
atman asymmetry \
    --de-results "$OUT"/de_results.tsv \
    --pairs "PT2-PT1:PR2-PR1,PT2-PR2:PT1-PR1" \
    --output "$OUT"/asymmetry.tsv

echo "[7/8] pathway enrichment (full + universe-restricted)"
scripts/run_enrichment.sh "$OUT"/de_results.tsv docs/findings/pathway_results/
scripts/enrich_restricted_universe.py "$OUT"/de_results.tsv \
    docs/findings/pathway_results_restricted/

echo "[8/8] phenotype regression (per-protein + per-pathway)"
scripts/phenotype_regression.py docs/findings/phenotype_results/
scripts/phenotype_pathway_regression.py "$OUT"/de_results.tsv \
    docs/findings/phenotype_results/

echo ""
echo "Done. Artefacts under: $OUT/, ${OUT}_mod/, docs/findings/"
