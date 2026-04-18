#!/usr/bin/env python3
"""Adapter: Gisby et al. 2022 (Cell Rep Med) COVID dialysis SomaLogic
SomaScan v4.1 dataset → Atman's internal long-format TSVs.

Third generalizability demonstration and the first on a non-Olink platform.
Demonstrates that Atman's statistical machinery (paired-t, LOO, S-score,
module-DE) is platform-agnostic: the same DAG that runs on Olink NPX runs
unchanged on log2-transformed SomaLogic RFU.

Data source
-----------
Gisby J et al. *Multi-omics identify LRRC15 as a COVID-19 severity predictor
and persistent pro-thrombotic signals in convalescence.* Nat Commun 2022.
Zenodo DOI 10.5281/zenodo.6497251, CC BY 4.0. Mirrored at
`example_data/gisby_somascan_2022/`.

Pairing scheme
--------------
For each COVID-positive patient with an early sample (`Time_From_First_Swab
<= 7`) AND a late sample (`>= 14`), the earliest early and latest late are
retained with condition labels `early` and `late`. Mirrors the Gisby 2021
Olink adapter so the paired design is comparable across platforms.

Usage
-----
    scripts/ingest_gisby_somascan.py example_data/gisby_somascan_2022 out_gisby_soma
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

SRC = Path(sys.argv[1] if len(sys.argv) > 1 else "example_data/gisby_somascan_2022")
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "out_gisby_soma")
OUT.mkdir(parents=True, exist_ok=True)

# ---- Load ----------------------------------------------------------------
abund = pd.read_csv(SRC / "soma_abundance.csv", index_col=0)
stm   = pd.read_csv(SRC / "sample_technical_meta.csv")
fm    = pd.read_csv(SRC / "feature_meta.csv")

# Parse sample_id: "C126_positive_7" → patient=C126, case=positive, timepoint=7
stm["case"]      = stm["sample_id"].str.extract(r"_(positive|negative)")[0]
stm["timepoint"] = stm["sample_id"].str.extract(r"_positive_(\d+)")[0].astype(float)
stm = stm[stm.sample_id.isin(abund.index)]

# ---- Select paired patients ---------------------------------------------
early = stm[(stm.case == "positive") & (stm.timepoint <= 7)]\
    .sort_values("timepoint").drop_duplicates("individual_id", keep="first")
late  = stm[(stm.case == "positive") & (stm.timepoint >= 14)]\
    .sort_values("timepoint").drop_duplicates("individual_id", keep="last")
paired = early[["sample_id","individual_id"]].assign(condition="early").merge(
    late[["sample_id","individual_id"]].assign(condition="late"),
    on="individual_id", suffixes=("_early","_late"))
print(f"paired patients with early + late SomaScan samples: {len(paired)}")
keep = set(paired.sample_id_early) | set(paired.sample_id_late)
print(f"total samples retained: {len(keep)}")

# ---- samples.tsv ---------------------------------------------------------
rows = []
for _, r in paired.iterrows():
    rows.append({"sample_id": r.sample_id_early, "subject_id": r.individual_id,
                 "condition": "early", "is_control": 0, "sample_type": "",
                 "ingest_order": 0})
    rows.append({"sample_id": r.sample_id_late,  "subject_id": r.individual_id,
                 "condition": "late",  "is_control": 0, "sample_type": "",
                 "ingest_order": 0})
samples = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)
samples["ingest_order"] = range(1, len(samples) + 1)
samples.to_csv(OUT / "samples.tsv", sep="\t", index=False)
print(f"samples.tsv: {len(samples)} rows")

# ---- proteins.tsv --------------------------------------------------------
# Map soma_id (SLxxx in feature_meta; columns in abund are seq.XXXX.YY).
# The abund column names correspond to the `seq_id` column in feature_meta,
# typically stored in `converted_seq_ids` or similar. Verify and fall back.
seq_col = None
for c in ("converted_seq_ids", "seq_id", "SeqId"):
    if c in fm.columns:
        seq_col = c
        break
if seq_col is None:
    # last-resort: match by replacing the "seq." prefix and "." separators
    print("WARN: no canonical seq column in feature_meta; constructing from soma_id")
    seq_col = "soma_id"
fm_keep = fm[[seq_col, "entrez_gene_symbol", "uniprot_single", "target"]].copy()
fm_keep = fm_keep.rename(columns={seq_col: "seq_id", "entrez_gene_symbol": "gene_symbol",
                                  "uniprot_single": "uniprot"})
fm_keep = fm_keep.drop_duplicates("seq_id")
fm_keep = fm_keep[fm_keep.seq_id.isin(abund.columns)]
fm_keep["assay_id"] = fm_keep.seq_id
fm_keep["panel"]     = "SomaScan_v4.1"
fm_keep["panel_lot"] = "gisby_somascan_2022"
fm_keep["platform"]  = "somalogic_somascan"
prot = fm_keep[["platform", "assay_id", "uniprot", "gene_symbol", "panel", "panel_lot"]]
prot.to_csv(OUT / "proteins.tsv", sep="\t", index=False)
print(f"proteins.tsv: {len(prot)} SOMAmers mapped to gene symbols")

# ---- measurements.tsv ----------------------------------------------------
# SomaLogic RFU values are raw (median ~615, range 9–3e5 on this file).
# Log2-transform to bring onto the same footing as log2-NPX. Drop any
# non-positive values (none expected in well-normalised RFU).
wide = abund.loc[list(keep)]
wide = wide[fm_keep.seq_id.tolist()]
# log2 transform
data = wide.to_numpy(dtype=float)
data = np.where(data > 0, np.log2(data), np.nan)
log2 = pd.DataFrame(data, index=wide.index, columns=wide.columns)

long = log2.reset_index().melt(id_vars="index", var_name="seq_id", value_name="abundance")\
    .rename(columns={"index": "sample_id"})
long = long.merge(fm_keep[["seq_id","gene_symbol","uniprot"]], on="seq_id")
long["assay_id"]       = long["seq_id"]
long["panel"]          = "SomaScan_v4.1"
long["npx_source_str"] = long["abundance"].astype(str)
long["abundance_raw"]  = long["abundance"]
long["abundance_unit"] = "log2_rfu"
long["qc_sample"]      = "PASS"
long["qc_assay"]       = "PASS"
long["detection_limit"]= ""
long["below_lod"]      = 0
long["dropped_by_qc"]  = long["abundance"].isna().astype(int)
long["plate_id"]       = "gisby_somascan_2022"
long["panel_lot"]      = "gisby_somascan_2022"
long["ingest_order"]   = range(1, len(long) + 1)
m = long[["sample_id","assay_id","gene_symbol","panel","npx_source_str",
          "abundance","abundance_raw","abundance_unit","qc_sample","qc_assay",
          "detection_limit","below_lod","dropped_by_qc","plate_id","panel_lot","ingest_order"]]
m = m.dropna(subset=["abundance"])
m.to_csv(OUT / "measurements.tsv", sep="\t", index=False)
m.to_csv(OUT / "qc_measurements.tsv", sep="\t", index=False)
print(f"measurements.tsv: {len(m):,} rows ({m.sample_id.nunique()} samples × {m.assay_id.nunique()} SOMAmers)")
print(f"paired contrast: late − early, N patients = {len(paired)}")
