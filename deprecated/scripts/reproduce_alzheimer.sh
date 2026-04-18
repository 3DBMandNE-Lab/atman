#!/usr/bin/env bash
# Reproduce the Wei et al. 2025 Alzheimer's plasma unpaired DE demonstration.
# Assumes `atman` on PATH (run `cargo install --path crates/atman` first) and
# Python 3.11+ with pandas and openpyxl.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"
export PATH="$REPO_ROOT/target/release:$PATH"

DATA=example_data/alzheimer_olink_2024
OUT=out_ad

# Fetch raw xlsx from figshare if missing.
if [[ ! -f "$DATA/NPX_values_Alzheimers_study_Samples.xlsx" ]]; then
    mkdir -p "$DATA"
    curl -sL -o "$DATA/NPX_values_Alzheimers_study_Samples.xlsx" \
        "https://ndownloader.figshare.com/files/53841275"
fi

echo "[1/2] ingest (Python adapter)"
python3 scripts/ingest_alzheimer.py "$DATA" "$OUT"

echo "[2/2] de (welch-t, three contrasts)"
atman de --input-dir "$OUT" --output-dir "$OUT" \
    --test welch-t --groups "AD-HC,MCI-HC,AD-MCI" --min-pairs 5

echo ""
echo "Done. Summary:"
python3 - <<'PY'
import pandas as pd
de = pd.read_csv("out_ad/de_results.tsv", sep="\t")
for c in ["AD-HC","MCI-HC","AD-MCI"]:
    sub = de[de.comparison==c]
    print(f"  {c}: tested={len(sub)}, q<0.05={(sub.bh_q<0.05).sum()}, q<0.10={(sub.bh_q<0.10).sum()}")
PY
