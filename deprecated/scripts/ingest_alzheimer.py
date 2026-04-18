#!/usr/bin/env python3
"""Adapter: Alzheimer plasma Olink Target 96 Inflammation dataset →
Atman internal TSVs.

Second generalizability demonstration and the first with an UNPAIRED design.
Three groups (A=AD n=10, B=MCI n=10, C=HC n=10) on the Inflammation panel
(92 assays), released under CC BY at figshare 28829720. Source format is
the Olink NPX Manager xlsx export ("Npx Data" sheet).

Usage
-----
    scripts/ingest_alzheimer.py example_data/alzheimer_olink_2024 out_ad
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd

SRC = Path(sys.argv[1] if len(sys.argv) > 1 else "example_data/alzheimer_olink_2024")
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "out_ad")
OUT.mkdir(parents=True, exist_ok=True)

xl = pd.ExcelFile(SRC / "NPX_values_Alzheimers_study_Samples.xlsx")
raw = pd.read_excel(xl, sheet_name="Npx Data", header=None)

# Fixed layout (verified by inspection):
#   row 3: "Assay"     — cols 1..92 are protein symbols, col 93 "Plate ID", col 94 "QC Warning"
#   row 4: "Uniprot ID" — UniProt accessions for the 92 proteins
#   row 5: "OlinkID"   — OlinkID per protein
#   rows 7..36: 30 samples, labelled group{A|B|C}{1..10}
#   rows 38+: LOD / normalization / legend (discarded)
ASSAY_ROW = 3
UNIPROT_ROW = 4
OLINKID_ROW = 5
SAMPLE_START = 7
SAMPLE_END = 37   # exclusive
N_PROTEINS = 92   # cols 1..92
PLATE_COL = 93
QC_COL = 94

assays   = raw.iloc[ASSAY_ROW,   1:1+N_PROTEINS].astype(str).tolist()
uniprots = raw.iloc[UNIPROT_ROW, 1:1+N_PROTEINS].astype(str).tolist()
olink_ids = raw.iloc[OLINKID_ROW, 1:1+N_PROTEINS].astype(str).tolist()

# ---- proteins.tsv ---------------------------------------------------------
prot = pd.DataFrame({
    "platform":     "olink_target_96",
    "assay_id":     olink_ids,
    "uniprot":      uniprots,
    "gene_symbol":  assays,
    "panel":        "Inflammation",
    "panel_lot":    "Inflammation_v3024",
})
prot.to_csv(OUT / "proteins.tsv", sep="\t", index=False)
print(f"proteins.tsv: {len(prot)} rows")

# ---- samples.tsv ----------------------------------------------------------
group_to_condition = {"A": "AD", "B": "MCI", "C": "HC"}
sample_rows = []
for r in range(SAMPLE_START, SAMPLE_END):
    sid = str(raw.iloc[r, 0])
    if not sid.startswith("group"):
        continue
    letter = sid[5]  # "A", "B", "C"
    condition = group_to_condition.get(letter)
    if condition is None:
        continue
    sample_rows.append({
        "sample_id":    sid,
        "subject_id":   sid,   # unpaired: each sample is its own subject
        "condition":    condition,
        "is_control":   0,
        "sample_type":  "",
        "ingest_order": 0,     # filled below
    })
samples = pd.DataFrame(sample_rows).sort_values("sample_id").reset_index(drop=True)
samples["ingest_order"] = range(1, len(samples) + 1)
samples.to_csv(OUT / "samples.tsv", sep="\t", index=False)
print(f"samples.tsv: {len(samples)} rows  "
      f"(conditions: {samples['condition'].value_counts().to_dict()})")

# ---- measurements.tsv -----------------------------------------------------
meas_rows = []
order = 0
for r in range(SAMPLE_START, SAMPLE_END):
    sid = str(raw.iloc[r, 0])
    if not sid.startswith("group"):
        continue
    plate  = str(raw.iloc[r, PLATE_COL]) if pd.notna(raw.iloc[r, PLATE_COL]) else ""
    qc_str = str(raw.iloc[r, QC_COL])    if pd.notna(raw.iloc[r, QC_COL])    else "Pass"
    # "Pass" in Olink Manager parlance == "PASS" in Atman's QC flag set.
    qc_sample = "PASS" if qc_str.strip().lower() == "pass" else "WARN"
    dropped = 0 if qc_sample == "PASS" else 1
    for col, (aid, gene, _up) in enumerate(zip(olink_ids, assays, uniprots), start=1):
        val = raw.iloc[r, col]
        if pd.isna(val):
            continue
        try:
            npx = float(val)
        except (TypeError, ValueError):
            continue
        order += 1
        meas_rows.append({
            "sample_id":       sid,
            "assay_id":        aid,
            "gene_symbol":     gene,
            "panel":           "Inflammation",
            "npx_source_str":  str(val),
            "abundance":       npx,
            "abundance_raw":   npx,
            "abundance_unit":  "log2_npx",
            "qc_sample":       qc_sample,
            "qc_assay":        "PASS",
            "detection_limit": "",
            "below_lod":       0,
            "dropped_by_qc":   dropped,
            "plate_id":        plate,
            "panel_lot":       "Inflammation_v3024",
            "ingest_order":    order,
        })

m = pd.DataFrame(meas_rows)
m.to_csv(OUT / "measurements.tsv", sep="\t", index=False)
# qc_measurements is the post-QC view. Rows where dropped_by_qc=1 keep the
# same abundance_raw but abundance is retained (atman downstream stages key
# off dropped_by_qc via effective_abundance()).
m.to_csv(OUT / "qc_measurements.tsv", sep="\t", index=False)
print(f"measurements.tsv: {len(m):,} rows ({m['sample_id'].nunique()} samples × "
      f"{m['assay_id'].nunique()} assays)")
n_dropped = (m['dropped_by_qc'] == 1).sum()
if n_dropped:
    print(f"  {n_dropped:,} rows flagged by QC Warning (Pass→PASS, Warning→WARN)")
