#!/usr/bin/env python3
"""Subject-bootstrap uncertainty at the module level.

For a paired contrast A−B and a module assignment, resample subjects
with replacement (B reps) and recompute per-module mean delta. Output
per module:
    - bootstrap mean (expected module effect)
    - bootstrap sd (uncertainty)
    - 95% CI (percentile)
    - P(m_k > 0) — fraction of bootstrap reps with positive module mean
    - n_genes in module

This replaces the S-score heuristic with a principled uncertainty
estimate at the scale where aggregation actually buys power (modules).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def load(input_dir: Path, contrast: str, modules_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (wide per-subject module deltas [subject × module], module_sizes)."""
    meas = pd.read_csv(input_dir / "qc_measurements.tsv", sep="\t")
    if "dropped_by_qc" in meas.columns:
        meas = meas[meas["dropped_by_qc"] == 0]
    samp = pd.read_csv(input_dir / "samples.tsv", sep="\t")
    mods = pd.read_csv(modules_path, sep="\t")

    cond_a, cond_b = contrast.split("-")
    sub = meas.merge(samp[["sample_id", "subject_id", "condition"]], on="sample_id")
    sub = sub[sub["condition"].isin([cond_a, cond_b])]

    # Per-(subject, gene, condition) mean abundance (collapse any technical reps).
    pivot = sub.pivot_table(
        index=["subject_id", "gene_symbol"], columns="condition",
        values="abundance", aggfunc="mean",
    ).dropna(subset=[cond_a, cond_b])
    pivot["delta"] = pivot[cond_a] - pivot[cond_b]
    pivot = pivot.reset_index()

    # Join modules, drop genes without a module.
    pivot = pivot.merge(mods, on="gene_symbol", how="inner")

    # Per-(subject, module) mean delta across member proteins.
    subj_mod = pivot.groupby(["subject_id", "module"])["delta"].mean().reset_index()
    wide = subj_mod.pivot(index="subject_id", columns="module", values="delta")
    wide = wide.dropna(axis=1, how="any")  # require every subject to have a value

    sizes = pivot.groupby("module")["gene_symbol"].nunique().rename("n_genes")
    return wide, sizes.to_frame()


def bootstrap(wide: pd.DataFrame, n_boot: int, rng: np.random.Generator) -> pd.DataFrame:
    subjects = wide.index.to_numpy()
    n = len(subjects)
    X = wide.to_numpy()  # (n_subjects, n_modules)
    draws = rng.integers(0, n, size=(n_boot, n))  # index matrix
    boot_means = X[draws].mean(axis=1)  # (n_boot, n_modules)
    rows = []
    for i, mod in enumerate(wide.columns):
        col = boot_means[:, i]
        rows.append({
            "module": mod,
            "n_subjects": n,
            "point_mean": float(X[:, i].mean()),
            "boot_mean": float(col.mean()),
            "boot_sd": float(col.std(ddof=1)),
            "ci_lo_95": float(np.quantile(col, 0.025)),
            "ci_hi_95": float(np.quantile(col, 0.975)),
            "p_pos": float((col > 0).mean()),
            "p_sign_stable": float(max((col > 0).mean(), (col < 0).mean())),
        })
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True,
                    help="directory with qc_measurements.tsv, samples.tsv")
    ap.add_argument("--modules", type=Path, required=True,
                    help="modules.tsv (module, gene_symbol)")
    ap.add_argument("--contrast", required=True, help="A-B")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    wide, sizes = load(args.input_dir, args.contrast, args.modules)
    if wide.shape[0] < 4 or wide.shape[1] == 0:
        print(f"[bootstrap] too few subjects or modules (n={wide.shape[0]}, "
              f"k={wide.shape[1]})", file=sys.stderr)
        return 1
    rng = np.random.default_rng(args.seed)
    out = bootstrap(wide, args.n_boot, rng)
    out = out.merge(sizes, left_on="module", right_index=True, how="left")
    out["contrast"] = args.contrast
    out["n_boot"] = args.n_boot
    out = out.sort_values("module").reset_index(drop=True)
    out.to_csv(args.output, sep="\t", index=False)
    n_stable = int((out["p_sign_stable"] >= 0.95).sum())
    print(f"wrote {args.output}  ({len(out)} modules, {n_stable} with P(sign)>=0.95)",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
