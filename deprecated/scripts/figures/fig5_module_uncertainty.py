#!/usr/bin/env python3
"""Fig 6 — module-level uncertainty via subject bootstrap.

Two standalone panels:
  fig6a_module_bootstrap_forest   per-module bootstrap mean + 95% CI,
                                  coloured by P(sign stable)
  fig6b_module_aggregation_benefit observed module bootstrap SD vs
                                  independent-noise prediction
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _style import (
    apply_style, save, set_title,
    PRIMARY, ACCENT, MUTED, SEQUENTIAL,
    PANEL_TALL, PANEL_SQUARE, REPO,
)

apply_style()


def per_protein_boot_sds(input_dir: Path, modules_path: Path, contrast: str,
                         n_boot: int, seed: int) -> pd.DataFrame:
    meas = pd.read_csv(input_dir / "qc_measurements.tsv", sep="\t")
    if "dropped_by_qc" in meas.columns:
        meas = meas[meas["dropped_by_qc"] == 0]
    samp = pd.read_csv(input_dir / "samples.tsv", sep="\t")
    mods = pd.read_csv(modules_path, sep="\t")
    cond_a, cond_b = contrast.split("-")
    sub = meas.merge(samp[["sample_id", "subject_id", "condition"]], on="sample_id")
    sub = sub[sub["condition"].isin([cond_a, cond_b])]
    pivot = sub.pivot_table(index=["subject_id", "gene_symbol"], columns="condition",
                            values="abundance", aggfunc="mean").dropna(subset=[cond_a, cond_b])
    pivot["delta"] = pivot[cond_a] - pivot[cond_b]
    pivot = pivot.reset_index().merge(mods, on="gene_symbol", how="inner")
    wide = pivot.pivot(index="subject_id", columns="gene_symbol",
                       values="delta").dropna(axis=1)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, wide.shape[0], size=(n_boot, wide.shape[0]))
    boot = wide.to_numpy()[idx].mean(axis=1)
    sds  = boot.std(axis=0, ddof=1)
    out = pd.DataFrame({"gene_symbol": wide.columns, "protein_boot_sd": sds})
    out = out.merge(mods, on="gene_symbol", how="inner")
    return out


def fig6a(bootstrap_tsv: Path, title: str):
    boot = pd.read_csv(bootstrap_tsv, sep="\t").sort_values("boot_mean")
    y = np.arange(len(boot))
    cmap = plt.get_cmap(SEQUENTIAL)
    norm = plt.Normalize(0.5, 1.0)
    colours = cmap(norm(boot["p_sign_stable"].to_numpy()))

    fig, ax = plt.subplots(figsize=(4.8, max(4.2, 0.32 * len(boot) + 1)))

    # Zero reference.
    ax.axvline(0, ls="--", lw=0.7, color="#4B5563", zorder=0)

    # Horizontal error bars, thin, muted.
    ax.errorbar(boot["boot_mean"], y,
                xerr=[boot["boot_mean"] - boot["ci_lo_95"],
                      boot["ci_hi_95"] - boot["boot_mean"]],
                fmt="none", ecolor="#9CA3AF",
                elinewidth=1.0, capsize=2.5, zorder=1)

    # Filled points at bootstrap mean, coloured by P(sign stable).
    ax.scatter(boot["boot_mean"], y, c=colours, s=72, zorder=3,
               edgecolor="white", linewidth=0.8)

    ax.set_yticks(y)
    ax.set_yticklabels([f"{m}  ({n})" for m, n in zip(boot["module"], boot["n_genes"])])
    ax.set_xlabel(r"bootstrap module effect  $\bar m_k$  (95% CI)")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.04, shrink=0.85)
    cbar.set_label(r"$P(\mathrm{sign\ stable})$", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    cbar.outline.set_linewidth(0)

    set_title(ax, title)
    save(fig, "fig5a_module_bootstrap_forest")
    plt.close(fig)


def fig6b(bootstrap_tsv: Path, input_dir: Path, modules_path: Path,
          contrast: str, n_boot: int, seed: int):
    boot = pd.read_csv(bootstrap_tsv, sep="\t")
    prot = per_protein_boot_sds(input_dir, modules_path, contrast, n_boot, seed)
    rms = prot.groupby("module")["protein_boot_sd"].apply(
        lambda s: float(np.sqrt(np.mean(s ** 2)))
    ).rename("rms_sd")
    boot = boot.merge(rms, on="module", how="left")
    boot["indep_sd"] = boot["rms_sd"] / np.sqrt(boot["n_genes"])

    fig, ax = plt.subplots(figsize=PANEL_SQUARE)
    lim = max(boot["boot_sd"].max(), boot["indep_sd"].max()) * 1.15

    # Independence prediction: observed == indep_sd.
    ax.plot([0, lim], [0, lim], ls="--", lw=0.7, color="#4B5563", zorder=0,
            label="independent-noise prediction")

    cmap = plt.get_cmap(SEQUENTIAL)
    norm = plt.Normalize(0.5, 1.0)
    colours = cmap(norm(boot["p_sign_stable"].to_numpy()))
    ax.scatter(boot["indep_sd"], boot["boot_sd"], c=colours, s=72,
               edgecolor="white", linewidth=0.8, zorder=3)

    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"predicted module SD if proteins independent   ($\sigma^{rms}/\sqrt{|G_k|}$)")
    ax.set_ylabel(r"observed module bootstrap SD   $\hat\sigma_k$")

    # Annotate modules with highest ratio (most correlated) and smallest (most independent).
    boot["ratio"] = boot["boot_sd"] / boot["indep_sd"].replace(0, np.nan)
    for _, r in boot.sort_values("ratio", ascending=False).head(3).iterrows():
        ax.annotate(r["module"].replace("module_", "m"),
                    (r["indep_sd"], r["boot_sd"]),
                    xytext=(5, 3), textcoords="offset points", fontsize=8,
                    color="#1F2937")

    ax.legend(loc="lower right")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.04, shrink=0.85)
    cbar.set_label(r"$P(\mathrm{sign\ stable})$", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    cbar.outline.set_linewidth(0)

    set_title(ax, "Aggregation benefit")
    save(fig, "fig5b_module_aggregation_benefit")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap-tsv", type=Path,
                    default=REPO / "out_gisby" / "robustness" / "module_bootstrap.tsv")
    ap.add_argument("--input-dir", type=Path, default=REPO / "out_gisby")
    ap.add_argument("--modules", type=Path,
                    default=REPO / "out_gisby" / "modules_learned.tsv")
    ap.add_argument("--contrast", default="late-early")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--title", default="Module-level uncertainty — Gisby late − early (n=23)")
    args = ap.parse_args()

    fig6a(args.bootstrap_tsv, args.title)
    fig6b(args.bootstrap_tsv, args.input_dir, args.modules,
          args.contrast, args.n_boot, args.seed)


if __name__ == "__main__":
    main()
