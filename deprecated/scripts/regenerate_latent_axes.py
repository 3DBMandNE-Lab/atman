#!/usr/bin/env python3
"""Regenerate subject latent axes + physiology mapping from atman outputs.

Feeds per-subject burden + module-trajectory features into z-scored SVD,
producing PC1/PC2 and subject archetype labels. Also recomputes Pearson
correlations between each metric (PC1, PC2, acute_burden, adapt_burden) and
Dube's physiological delta endpoints with BH correction within each metric
family.

Inputs:
  docs/findings/heterogeneity/subject_response_burden.tsv
  docs/findings/heterogeneity/module_trajectory_scores_v2.tsv
  example_data/dube_heat_2023/Physiological_data.xlsx

Outputs:
  docs/findings/heterogeneity/latent_axes_v2.tsv
  docs/findings/heterogeneity/physiology_mapping_v2.tsv
  stdout summary with PC1/PC2 percent variance + top PC1-phenotype corr
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parent.parent
HETDIR = ROOT / "docs" / "findings" / "heterogeneity"

burden = pd.read_csv(HETDIR / "subject_response_burden.tsv", sep="\t")
scores = pd.read_csv(HETDIR / "module_trajectory_scores_v2.tsv", sep="\t")

# Wide module-trajectory feature matrix: one row per subject, columns per
# (module, contrast) pair. Four contrast columns per module → 16 features.
contrast_cols = ["pt1_pr1", "pt2_pr2", "pt2_pt1", "pr2_pr1"]
long = scores.melt(
    id_vars=["subject_id", "module"],
    value_vars=contrast_cols,
    var_name="contrast",
    value_name="score",
)
long["feature"] = long["module"] + "__" + long["contrast"]
mod_wide = long.pivot(index="subject_id", columns="feature", values="score")

# Burden features
b = burden.set_index("subject_id")[["acute_burden", "adapt_burden"]]

feats = b.join(mod_wide, how="inner")
feats = feats.sort_index()
subjects = feats.index.tolist()
Xnames = feats.columns.tolist()
X = feats.to_numpy(dtype=float)

# Z-score columns (ddof=0, matches most plain SVD pipelines)
mu  = X.mean(axis=0)
std = X.std(axis=0, ddof=0)
std = np.where(std == 0, 1.0, std)
Z = (X - mu) / std

# SVD on centred matrix
U, S, Vt = np.linalg.svd(Z, full_matrices=False)
# Variance explained per component
var_expl = (S ** 2) / np.sum(S ** 2)
pc1_var = var_expl[0] * 100
pc2_var = var_expl[1] * 100

# Project onto first two components
scores_pc = U[:, :2] * S[:2]
PC1 = scores_pc[:, 0]
PC2 = scores_pc[:, 1]

# Archetype labels from sign splits of PC1, PC2, acute, adapt
acute = b["acute_burden"].reindex(subjects).to_numpy()
adapt = b["adapt_burden"].reindex(subjects).to_numpy()
amed = np.median(acute)
dmed = np.median(adapt)
p1med = np.median(PC1)
p2med = np.median(PC2)
archetype = [
    ("HighAxis1_HighAxis2" if p1 >= p1med else "LowAxis1_HighAxis2") + "|"
    + ("HighAxis1_HighAxis2"[5:] if False else "")
    for p1 in PC1
]
archetype = []
for i, s in enumerate(subjects):
    t1 = "HighAxis1" if PC1[i] >= p1med else "LowAxis1"
    t2 = "HighAxis2" if PC2[i] >= p2med else "LowAxis2"
    ac = "HighAcute"   if acute[i] >= amed else "LowAcute"
    ad = "HighAdapt"   if adapt[i] >= dmed else "LowAdapt"
    archetype.append(f"{t1}_{t2}|{ac}_{ad}")

la = pd.DataFrame({
    "subject_id":    subjects,
    "PC1":           PC1,
    "PC2":           PC2,
    "acute_burden":  acute,
    "adapt_burden":  adapt,
    "archetype":     archetype,
})
la.to_csv(HETDIR / "latent_axes_v2.tsv", sep="\t", index=False)

# PC loadings (v1 and v2)
load = pd.DataFrame({
    "feature": Xnames,
    "PC1_loading": Vt[0, :],
    "PC2_loading": Vt[1, :],
})
load.to_csv(HETDIR / "pc_loadings_v2.tsv", sep="\t", index=False)

print(f"latent_axes_v2.tsv written — {len(subjects)} subjects")
print(f"PC1 variance explained: {pc1_var:.2f}%")
print(f"PC2 variance explained: {pc2_var:.2f}%")
print(f"Top 5 PC1 loadings by |value|:")
top = load.reindex(load["PC1_loading"].abs().sort_values(ascending=False).index).head(5)
print(top.to_string(index=False))

# ---- Physiology mapping ---------------------------------------------------
# Reuse Dube physiology sheet parser from scripts/phenotype_regression.py
import sys as _sys
_sys.path.insert(0, str(ROOT / "scripts"))
from phenotype_regression import load_phenotypes, build_endpoints, PARTICIPANTS

pheno    = load_phenotypes()
endpoints = build_endpoints(pheno)   # {endpoint: {participant: delta}}
print(f"\nPhysiology endpoints loaded: {len(endpoints)} ({list(endpoints)[:3]}...)")

metrics_df = la.set_index("subject_id")[["PC1","PC2","acute_burden","adapt_burden"]]
rows = []
for metric_name in metrics_df.columns:
    metric_vals = metrics_df[metric_name]
    for ep_name, ep_map in endpoints.items():
        pairs = []
        for pid in PARTICIPANTS:
            if pid in metric_vals.index and pid in ep_map:
                mv = metric_vals.loc[pid]
                ev = ep_map[pid]
                if pd.notna(mv) and pd.notna(ev):
                    pairs.append((mv, ev))
        n = len(pairs)
        if n < 3:
            rows.append({"metric":metric_name, "phenotype":ep_name, "n":n,
                         "r":float("nan"), "p_value":float("nan"), "bh_q":float("nan")})
            continue
        xs = np.array([p[0] for p in pairs], dtype=float)
        ys = np.array([p[1] for p in pairs], dtype=float)
        if xs.std() == 0 or ys.std() == 0:
            rows.append({"metric":metric_name, "phenotype":ep_name, "n":n,
                         "r":float("nan"), "p_value":float("nan"), "bh_q":float("nan")})
            continue
        r, p = pearsonr(xs, ys)
        rows.append({"metric":metric_name, "phenotype":ep_name, "n":n,
                     "r":r, "p_value":p, "bh_q":float("nan")})

phys_map = pd.DataFrame(rows)

def bh_within(group):
    pv = group["p_value"].to_numpy()
    mask = ~np.isnan(pv)
    m = mask.sum()
    if m == 0:
        group["bh_q"] = float("nan")
        return group
    q = np.full_like(pv, float("nan"), dtype=float)
    order = np.argsort(pv[mask])
    sorted_p = pv[mask][order]
    ranks = np.arange(1, m+1)
    q_sorted = sorted_p * m / ranks
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)
    q_vals = np.empty(m)
    q_vals[order] = q_sorted
    q[mask] = q_vals
    group["bh_q"] = q
    return group

phys_map = phys_map.groupby("metric", group_keys=False).apply(bh_within)
phys_map = phys_map.sort_values(["metric","bh_q","p_value"])
phys_map.to_csv(HETDIR / "physiology_mapping_v2.tsv", sep="\t", index=False)
print(f"physiology_mapping_v2.tsv written — {len(phys_map)} metric-endpoint rows")

for metric_name in ("PC1","PC2","acute_burden","adapt_burden"):
    top = phys_map[phys_map["metric"]==metric_name].sort_values("p_value").head(3)
    print(f"\nTop 3 {metric_name}-phenotype raw associations:")
    for _, row in top.iterrows():
        print(f"  {row['phenotype']:14s}  n={row['n']:2d}  r={row['r']:+.4f}  p={row['p_value']:.4f}  q={row['bh_q']:.4f}")
