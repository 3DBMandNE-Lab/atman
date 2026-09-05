#!/usr/bin/env python3
"""Inputs + reference for tests/de_csf_extensions.rs: complete-case per-protein
OLS with controls flagged is_control=1, one missing cell, pooled-SD Cohen d."""
from pathlib import Path
import numpy as np, pandas as pd, statsmodels.api as sm

OUT = Path(__file__).parent / "de_csf"; OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(11)
n = 12
samples = pd.DataFrame(dict(
    sample_id=[f"S{i:02d}" for i in range(n)], subject_id=[f"S{i:02d}" for i in range(n)],
    condition=["Case"] * 6 + ["Ctrl"] * 6, is_control=[0] * 6 + [1] * 6, sample_type="bio",
    ingest_order=range(1, n + 1), age=rng.integers(40, 80, n).astype(float), sex=list("FMFMFM" * 2),
    site=["X", "X", "Y", "Y", "X", "Y"] * 2, grp=["A", "A", "A", "B", "B", "B", "C", "C", "C", "D", "D", "D"]))
samples.to_csv(OUT / "samples.tsv", sep="\t", index=False)
genes = ["G1", "G2", "G3"]
wide = pd.DataFrame({g: rng.normal(0, 1, n) + (0.9 if g == "G1" else 0.0) * (samples.condition == "Case").values for g in genes}, index=samples.sample_id)
wide.loc["S02", "G2"] = np.nan                      # one missing cell
wide.loc[["S00", "S01", "S02", "S03", "S04"], "G3"] = np.nan  # 5/12 missing => dropped at --max-missing-fraction 0.3
rows = []
for sid in wide.index:
    for g in genes:
        v = wide.loc[sid, g]
        if np.isnan(v):
            continue
        rows.append(dict(platform="maxquant_lfq", sample_id=sid, assay_id=g, gene_symbol=g, panel="P", npx_source_str=f"{v:.6f}",
                         abundance=v, abundance_raw=v, abundance_unit="log2_intensity", qc_sample="PASS", qc_assay="PASS",
                         detection_limit="", below_lod=0, dropped_by_qc=0, plate_id="", panel_lot="", ingest_order=len(rows) + 1))
pd.DataFrame(rows).to_csv(OUT / "measurements.tsv", sep="\t", index=False)
pd.DataFrame(dict(platform="maxquant_lfq", assay_id=genes, uniprot=["P1", "P2", "P3"], gene_symbol=genes, panel="P", panel_lot="")).to_csv(OUT / "proteins.tsv", sep="\t", index=False)

meta = samples.set_index("sample_id")
ref = []
for g in genes:
    y = wide[g]; ok = y.notna()
    m = meta[ok]; disease = (m.condition == "Case").astype(float)
    X = pd.DataFrame(dict(const=1.0, disease=disease, age=m.age, sexM=(m.sex == "M").astype(float)))
    fit = sm.OLS(y[ok].astype(float), X).fit()
    a, b = y[ok][disease == 1].values, y[ok][disease == 0].values
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    ref.append(dict(gene_symbol=g, beta=fit.params["disease"], p_value=fit.pvalues["disease"], cohen_d=(a.mean() - b.mean()) / sp,
                    n_a=len(a), n_b=len(b), n_obs=int(ok.sum())))
pd.DataFrame(ref).to_csv(OUT / "reference.tsv", sep="\t", index=False, float_format="%.10f")

# Continuous contrast: abundance ~ log10(age) + sex over all 12 subjects (complete case per gene).
cont = []
for g in genes:
    y = wide[g]; ok = y.notna(); m = meta[ok]
    X = pd.DataFrame(dict(const=1.0, log10_age=np.log10(m.age), sexM=(m.sex == "M").astype(float)))
    fit = sm.OLS(y[ok].astype(float), X).fit()
    cont.append(dict(gene_symbol=g, beta=fit.params["log10_age"], t=fit.tvalues["log10_age"], p_value=fit.pvalues["log10_age"], n_obs=int(ok.sum())))
pd.DataFrame(cont).to_csv(OUT / "reference_continuous.tsv", sep="\t", index=False, float_format="%.10f")
print("wrote", OUT)
