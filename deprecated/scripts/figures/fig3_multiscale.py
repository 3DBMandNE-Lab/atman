#!/usr/bin/env python3
"""Fig 4 — multiscale inference.

Three standalone panels:
  fig4a_multiscale_yield   per-protein vs module-level q<0.05 hit rate
  fig4b_module14_heatmap   Gisby module_14 (14 genes × 23 patients) deltas
  fig4c_module_volcano     Gisby module-level volcano (K=15)
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _style import (
    apply_style, save, set_title,
    PRIMARY, ACCENT, MUTED, HIGHLIGHT,
    DIVERGING, PANEL_SIZE, PANEL_SQUARE, PANEL_WIDE, PANEL_TALL, REPO,
)

apply_style()


def fig4a():
    dube     = pd.read_csv(REPO / "out" / "de_results.tsv", sep="\t")
    dube_mod = pd.read_csv(REPO / "out" / "module_de_results.tsv", sep="\t")
    gis      = pd.read_csv(REPO / "out_gisby" / "de_results.tsv", sep="\t")
    gis_mod  = pd.read_csv(REPO / "out_gisby" / "module_de_results.tsv", sep="\t")

    rows = []
    for label, p, m in [
        ("Dube\nPT1-PR1", dube[dube.comparison == "PT1-PR1"],
                          dube_mod[dube_mod.contrast == "PT1-PR1"]),
        ("Dube\nPT2-PR2", dube[dube.comparison == "PT2-PR2"],
                          dube_mod[dube_mod.contrast == "PT2-PR2"]),
        ("Gisby\nlate-early", gis[gis.comparison == "late-early"],
                              gis_mod[gis_mod.contrast == "late-early"]),
    ]:
        p_rate = (p.bh_q < 0.05).sum() / max(len(p), 1) * 100
        m_rate = (m.bh_q < 0.05).sum() / max(len(m), 1) * 100
        rows.append((label, p_rate, m_rate, len(p), len(m)))

    fig, ax = plt.subplots(figsize=PANEL_SIZE)
    labels = [r[0] for r in rows]
    x = np.arange(len(labels))
    w = 0.38
    prot_vals = [r[1] for r in rows]
    mod_vals  = [r[2] for r in rows]

    bars_p = ax.bar(x - w/2, prot_vals, w, color=MUTED,
                    edgecolor="white", linewidth=0.6, label="per-protein")
    bars_m = ax.bar(x + w/2, mod_vals,  w, color=PRIMARY,
                    edgecolor="white", linewidth=0.6, label="per-module")

    for b, (_, pv, _, np_, _) in zip(bars_p, rows):
        ax.annotate(f"{int(round(pv/100*np_))}/{np_}",
                    xy=(b.get_x()+b.get_width()/2, b.get_height()),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", fontsize=8, color="#374151")
    for b, (_, _, mv, _, nm) in zip(bars_m, rows):
        ax.annotate(f"{int(round(mv/100*nm))}/{nm}",
                    xy=(b.get_x()+b.get_width()/2, b.get_height()),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", fontsize=8, color="#374151")

    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("hit rate at q < 0.05  (%)")
    ax.set_ylim(0, max(prot_vals + mod_vals) * 1.30)
    ax.legend(loc="upper left")
    set_title(ax, "Per-protein vs module-level hit yield")
    save(fig, "fig3a_multiscale_yield")
    plt.close(fig)


def fig4b():
    mods   = pd.read_csv(REPO / "out_gisby" / "modules_learned.tsv", sep="\t")
    m14    = mods[mods.module == "module_14"].gene_symbol.tolist()
    deltas = pd.read_csv(REPO / "out_gisby" / "per_subject_gene_deltas.tsv", sep="\t")
    sub    = deltas[deltas.gene_symbol.isin(m14)]
    wide   = sub.pivot_table(index="gene_symbol", columns="subject_id",
                             values="late_early", aggfunc="mean")
    wide   = wide.reindex(wide.mean(axis=1).sort_values().index)
    wide   = wide.reindex(columns=wide.mean(axis=0).sort_values().index)
    data   = wide.to_numpy()
    vmax   = float(np.nanpercentile(np.abs(data), 95))

    # Tall panel: 14 rows, 23 cols → wide heatmap works best
    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    im = ax.imshow(data, aspect="auto", cmap=DIVERGING,
                   vmin=-vmax, vmax=vmax, interpolation="nearest")

    ax.set_xticks(range(wide.shape[1]))
    ax.set_xticklabels(wide.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(wide.shape[0]))
    ax.set_yticklabels(wide.index, fontsize=8)
    ax.set_xlabel("patient  (n = 23)")
    ax.set_ylabel("")
    ax.tick_params(axis="both", length=2, width=0.5)

    cbar = plt.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
    cbar.set_label(r"late $-$ early  ($\log_2$ NPX)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    cbar.outline.set_linewidth(0)

    set_title(ax, "Gisby module_14  (data-driven, 14 genes)")
    save(fig, "fig3b_module14_heatmap")
    plt.close(fig)


def fig4c():
    de = pd.read_csv(REPO / "out_gisby" / "module_de_results.tsv", sep="\t")
    de = de[de.contrast == "late-early"].dropna(
        subset=["bh_q", "mean_diff", "n_genes"]).copy()
    de["neglogq"] = -np.log10(de["bh_q"].clip(lower=1e-20))
    sig = de.bh_q < 0.05
    sizes = 20 + 80 * np.log10(de.n_genes.clip(lower=1))

    fig, ax = plt.subplots(figsize=PANEL_SQUARE)
    ax.axhline(-np.log10(0.05), ls="--", lw=0.7, color="#4B5563", zorder=0)
    ax.axvline(0, ls=":", lw=0.5, color="#9CA3AF", zorder=0)

    ax.scatter(de.loc[~sig, "mean_diff"], de.loc[~sig, "neglogq"],
               s=sizes[~sig], color=MUTED, edgecolor="white",
               linewidth=0.4, alpha=0.85, label="q ≥ 0.05")
    ax.scatter(de.loc[sig, "mean_diff"], de.loc[sig, "neglogq"],
               s=sizes[sig], color=PRIMARY, edgecolor="white",
               linewidth=0.4, label="q < 0.05")

    # Manual placement per label to avoid edge collisions.
    label_offsets = {
        "module_14": (-12, 10),
        "module_04": (8, 4),
        "module_03": (-60, 0),
        "module_09": (-62, -2),
    }
    for _, r in de[sig].nsmallest(4, "bh_q").iterrows():
        dx, dy = label_offsets.get(r["module"], (7, 3))
        ax.annotate(r["module"],
                    (r.mean_diff, r.neglogq),
                    xytext=(dx, dy), textcoords="offset points",
                    fontsize=8, color="#1F2937")

    # Widen x-axis so the rightmost point (module_09) has room.
    xmax = de.mean_diff.max() * 1.25
    xmin = de.mean_diff.min() * 1.10
    ax.set_xlim(xmin, xmax)

    ax.set_xlabel(r"module mean difference  (late $-$ early)")
    ax.set_ylabel(r"$-\log_{10}(q)$")
    ax.legend(loc="lower right")
    set_title(ax, "Gisby module-level volcano  (K = 15)")
    save(fig, "fig3c_module_volcano")
    plt.close(fig)


def main():
    fig4a()
    fig4b()
    fig4c()


if __name__ == "__main__":
    main()
