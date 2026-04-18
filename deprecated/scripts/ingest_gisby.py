#!/usr/bin/env python3
"""Adapter: Gisby et al. 2021 (eLife) COVID dialysis Olink Target 96 dataset
→ Atman's internal long-format TSVs (`measurements.tsv`, `qc_measurements.tsv`,
`proteins.tsv`, `samples.tsv`) so that downstream `atman matrix /
fold-change / de / asymmetry / robustness` stages run unchanged.

This is the concrete generalizability demonstration: a second, independently
collected Olink dataset (Target 96 plasma; n=105 patients; 307 serial
samples) ingested via a ~100-line adapter and analysed with the same DAG
used for Dube.

Data source
-----------
Gisby J et al. *Longitudinal proteomic profiling of dialysis patients with
COVID-19 reveals markers of severity and predictors of death.* eLife 2021.
https://doi.org/10.7554/eLife.64827. Code and data mirrored at
https://github.com/jackgisby/longitudinal_olink_proteomics.

Pairing scheme
--------------
For each patient with at least one sample at `Time_From_First_Swab ≤ 7`
(early) AND at least one at `Time_From_First_Swab ≥ 14` (late), the
earliest early sample and the latest late sample are retained with
condition labels `early` and `late`. Mid-window samples are excluded. This
mirrors Dube's paired-contrast structure: within-subject comparison across
two temporally separated states.

Usage
-----
    scripts/ingest_gisby.py example_data/gisby_covid_2021 out_gisby
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd

SRC = Path(sys.argv[1] if len(sys.argv) > 1 else "example_data/gisby_covid_2021")
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "out_gisby")
OUT.mkdir(parents=True, exist_ok=True)

# ---- Load raw Gisby CSVs --------------------------------------------------
sm = pd.read_csv(SRC / "plasma_sample_level.csv")
npx = pd.read_csv(SRC / "plasma_npx_level.csv")
sm["Time_From_First_Swab"] = pd.to_numeric(sm["Time_From_First_Swab"], errors="coerce")

# ---- Pick one early + one late sample per patient with both --------------
early = sm[sm["Time_From_First_Swab"] <= 7].copy()
late  = sm[sm["Time_From_First_Swab"] >= 14].copy()
early_pick = early.sort_values("Time_From_First_Swab").drop_duplicates("Individual_ID", keep="first")
late_pick  = late .sort_values("Time_From_First_Swab").drop_duplicates("Individual_ID", keep="last")
paired = early_pick[["SampleID","Individual_ID","Plate_ID"]].assign(condition="early").merge(
    late_pick[["SampleID","Individual_ID","Plate_ID"]].assign(condition="late"),
    on="Individual_ID", suffixes=("_early","_late")
)
print(f"Patients with both early and late samples: {len(paired)}")
keep_samples = set(paired["SampleID_early"]).union(paired["SampleID_late"])
print(f"Total samples retained: {len(keep_samples)}")

# ---- Build samples.tsv ----------------------------------------------------
rows = []
for _, r in paired.iterrows():
    rows.append({"sample_id": r["SampleID_early"], "subject_id": r["Individual_ID"],
                 "condition": "early", "is_control": 0, "sample_type": "",
                 "ingest_order": len(rows)+1})
    rows.append({"sample_id": r["SampleID_late"],  "subject_id": r["Individual_ID"],
                 "condition": "late",  "is_control": 0, "sample_type": "",
                 "ingest_order": len(rows)+1})
samples = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)
samples["ingest_order"] = range(1, len(samples)+1)
samples.to_csv(OUT / "samples.tsv", sep="\t", index=False)
print(f"samples.tsv: {len(samples)} rows")

# ---- Build proteins.tsv (one row per OlinkID on Gisby's Target 96 panel) -
prot = (npx[["UniProt","GeneID","Assay","Panel"]]
        .drop_duplicates()
        .reset_index(drop=True))
prot["platform"]   = "olink_target_96"
prot["assay_id"]   = "OID_gisby_" + prot.index.astype(str).str.zfill(5)
prot["panel_lot"]  = "gisby_covid_2021"
prot = prot.rename(columns={"UniProt":"uniprot","Assay":"gene_symbol","Panel":"panel"})[
    ["platform","assay_id","uniprot","gene_symbol","panel","panel_lot"]
]
prot.to_csv(OUT / "proteins.tsv", sep="\t", index=False)
print(f"proteins.tsv: {len(prot)} rows")

# ---- Build measurements.tsv (long-format) --------------------------------
# Lookup gene_symbol → assay_id (to keep one row per Assay × SampleID)
gene_to_assay = dict(zip(prot["gene_symbol"], prot["assay_id"]))
gene_to_uniprot = dict(zip(prot["gene_symbol"], prot["uniprot"]))
gene_to_panel = dict(zip(prot["gene_symbol"], prot["panel"]))

keep_npx = npx[npx["SampleID"].isin(keep_samples)].copy()
plate_map = dict(zip(sm["SampleID"], sm["Plate_ID"]))

m = pd.DataFrame({
    "sample_id":         keep_npx["SampleID"].values,
    "assay_id":          keep_npx["Assay"].map(gene_to_assay).values,
    "gene_symbol":       keep_npx["Assay"].values,
    "panel":             keep_npx["Panel"].values,
    "npx_source_str":    keep_npx["NPX"].astype(str).values,
    "abundance":         keep_npx["NPX"].values,
    "abundance_raw":     keep_npx["NPX"].values,
    "abundance_unit":    "log2_npx",
    "qc_sample":         "PASS",
    "qc_assay":          "PASS",
    "detection_limit":   "",
    "below_lod":         0,
    "dropped_by_qc":     0,
    "plate_id":          keep_npx["SampleID"].map(plate_map).values,
    "panel_lot":         "gisby_covid_2021",
    "ingest_order":      range(1, len(keep_npx)+1),
})
# Drop rows with missing NPX (limma-style NA) so downstream stages don't choke.
m = m.dropna(subset=["abundance"])
m.to_csv(OUT / "measurements.tsv", sep="\t", index=False)
# Gisby NPX values are already QC-filtered in the published file; qc_measurements
# is identical to measurements here.
m.to_csv(OUT / "qc_measurements.tsv", sep="\t", index=False)
print(f"measurements.tsv: {len(m)} rows ({m['sample_id'].nunique()} samples × "
      f"{m['assay_id'].nunique()} assays)")
print(f"paired contrast: late − early, N patients = {len(paired)}")
