#!/usr/bin/env bash
# Reproduce the Gisby 2022 SomaLogic SomaScan v4.1 cross-platform demonstration.
# Assumes `atman` on PATH (run `cargo install --path crates/atman` first) and
# Python 3.11+ with pandas, numpy, scipy.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"
export PATH="$REPO_ROOT/target/release:$PATH"

DATA=example_data/gisby_somascan_2022
OUT=out_gisby_soma

# Fetch raw data from Zenodo if missing.
if [[ ! -f "$DATA/soma_abundance.csv" ]]; then
    mkdir -p "$DATA"
    base="https://zenodo.org/api/records/6497251/files"
    for f in soma_abundance.csv feature_meta.csv sample_technical_meta.csv w1_metadata.csv w2_metadata.csv; do
        curl -sL -o "$DATA/$f" "${base}/${f}/content"
    done
fi

echo "[1/4] ingest (Python adapter, log2-transform RFU)"
python3 scripts/ingest_gisby_somascan.py "$DATA" "$OUT"

echo "[2/4] de (paired-t, late - early, N=28 patients)"
atman de --input-dir "$OUT" --output-dir "$OUT" \
    --test paired-t --paired-by participant \
    --groups "late-early" --min-pairs 10

echo "[3/4] learn modules (K=25, complete linkage on per-subject delta correlations)"
python3 - <<'PY'
import pandas as pd
m = pd.read_csv("out_gisby_soma/qc_measurements.tsv", sep="\t")
s = pd.read_csv("out_gisby_soma/samples.tsv", sep="\t")
m = m.merge(s[["sample_id","subject_id","condition"]], on="sample_id")
wide = m[m.dropped_by_qc == 0].pivot_table(
    index=["subject_id","gene_symbol"], columns="condition",
    values="abundance", aggfunc="mean").reset_index()
wide["late_early"] = wide["late"] - wide["early"]
wide.dropna(subset=["late_early"])[["subject_id","gene_symbol","late_early"]] \
    .to_csv("out_gisby_soma/per_subject_gene_deltas.tsv", sep="\t", index=False)
PY
python3 scripts/learn_modules.py out_gisby_soma/per_subject_gene_deltas.tsv \
    out_gisby_soma/modules_learned.tsv --k-modules 25 \
    --summary out_gisby_soma/modules_learned_summary.tsv

echo "[4/4] module-DE (paired-t over K=25 modules)"
atman module-de --input-dir "$OUT" --output-dir "$OUT" \
    --modules-tsv "$OUT/modules_learned.tsv" \
    --test paired-t --paired-by participant \
    --groups "late-early" --min-pairs 10

echo ""
echo "Done. Summary:"
python3 - <<'PY'
import pandas as pd
de = pd.read_csv("out_gisby_soma/de_results.tsv", sep="\t")
md = pd.read_csv("out_gisby_soma/module_de_results.tsv", sep="\t")
print(f"  Per-protein DE: {len(de)} tests, q<0.05: {(de.bh_q<0.05).sum()}, q<0.10: {(de.bh_q<0.10).sum()}")
print(f"  Module-level DE: {len(md)} modules, q<0.05: {(md.bh_q<0.05).sum()}, q<0.10: {(md.bh_q<0.10).sum()}")
PY
