#!/usr/bin/env python3
"""Fig 3 — robustness geometry: Dube (n=10) vs Gisby (n=23).

Two panels:
  A. Distribution of per-protein LOO sign-match rates (SSR) across all
     proteins, split by Dube contrast and Gisby late-early. Shows that
     directional stability is high everywhere, and comparable at n=10 and
     n=23.
  B. Strict-threshold hit retention: of baseline q<0.05 hits, what fraction
     are retained at q<0.05 in ALL LOO reruns, vs ≥80% of runs, vs ≥50%.
     Shows that the margin is fragile at n=10 Dube and tight at n=23 Gisby
     — the quantitative expression of "stability scales with n".

Outputs: docs/findings/figures/fig3_robustness.{pdf,png}
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "docs" / "findings" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def panel_a(ax):
    """Violin of per-protein LOO SSR distributions (atman output schema)."""
    dube = pd.read_csv(REPO / "out/robustness/loo_sign_stability.tsv", sep="\t")
    gis  = pd.read_csv(REPO / "out_gisby/robustness/loo_sign_stability.tsv", sep="\t")
    groups, labels = [], []
    for c in ("PT1-PR1", "PT2-PR2", "PT2-PT1", "PR2-PR1"):
        sub = dube[dube.comparison == c]["sign_match_rate"].dropna()
        groups.append(sub.values)
        labels.append(f"Dube\n{c}\n(n=10)")
    sub = gis[gis.comparison == "late-early"]["sign_match_rate"].dropna()
    groups.append(sub.values)
    labels.append("Gisby\nlate-early\n(n=23)")

    parts = ax.violinplot(groups, positions=np.arange(len(groups)),
                          widths=0.75, showmeans=False, showmedians=True)
    for i, pc in enumerate(parts['bodies']):
        color = "#F59E0B" if i < 4 else "#1D4ED8"
        pc.set_facecolor(color); pc.set_edgecolor("black"); pc.set_linewidth(0.4); pc.set_alpha(0.6)
    for key in ("cmins", "cmaxes", "cbars", "cmedians"):
        if key in parts:
            parts[key].set_color("black"); parts[key].set_linewidth(0.6)
    # Overlay mean markers
    for i, g in enumerate(groups):
        ax.plot(i, np.nanmean(g), "o", color="black", markersize=3)
    ax.set_xticks(np.arange(len(groups))); ax.set_xticklabels(labels, fontsize=6.5)
    ax.set_ylabel("per-protein LOO sign-match rate")
    ax.set_ylim(0, 1.05)
    ax.axhline(1.0, ls=":", color="grey", linewidth=0.4, zorder=0)
    ax.set_title("A", loc="left", fontweight="bold", fontsize=10, pad=6)
    ax.text(0.5, 1.04, "directional stability across proteins and datasets",
            transform=ax.transAxes, ha="center", fontsize=8)


def panel_b(ax):
    """Strict-threshold retention of baseline q<0.05 hits under LOO."""
    # Dube: atman-format loo_sign_stability has baseline_q and n_retained_q_lt_05.
    dube = pd.read_csv(REPO / "out/robustness/loo_sign_stability.tsv", sep="\t")
    gis  = pd.read_csv(REPO / "out_gisby/robustness/loo_sign_stability.tsv", sep="\t")

    def compute_fracs(df, contrast):
        sub = df[(df.comparison == contrast) & (df.baseline_q.notna()) & (df.baseline_q < 0.05)]
        K = len(sub)
        if K == 0:
            return K, 0, 0, 0
        n_loo = int(sub["n_loo"].iloc[0]) if "n_loo" in sub.columns else 10
        frac_all = (sub["n_retained_q_lt_05"] == n_loo).sum() / K * 100
        frac_80  = (sub["n_retained_q_lt_05"] >= 0.8 * n_loo).sum() / K * 100
        frac_50  = (sub["n_retained_q_lt_05"] >= 0.5 * n_loo).sum() / K * 100
        return K, frac_all, frac_80, frac_50

    rows = []
    for c, label in [("PT1-PR1", "Dube\nPT1-PR1"), ("PT2-PR2", "Dube\nPT2-PR2"),
                     ("PT2-PT1", "Dube\nPT2-PT1"), ("PR2-PR1", "Dube\nPR2-PR1")]:
        K, a, b, d = compute_fracs(dube, c)
        rows.append((label, K, a, b, d))
    K, a, b, d = compute_fracs(gis, "late-early")
    rows.append(("Gisby\nlate-early", K, a, b, d))

    labels = [f"{r[0]}\n(K={r[1]})" for r in rows]
    x = np.arange(len(rows))
    w = 0.26
    colors = ["#1D4ED8", "#60A5FA", "#BFDBFE"]
    for i, (tag, key_idx) in enumerate([("all LOO", 2), ("≥80% LOO", 3), ("≥50% LOO", 4)]):
        vals = [r[key_idx] for r in rows]
        bars = ax.bar(x + (i - 1) * w, vals, w, color=colors[i], label=tag,
                      edgecolor="black", linewidth=0.4)
        for b, v in zip(bars, vals):
            if v > 5:
                ax.annotate(f"{v:.0f}%", xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                            xytext=(0, 2), textcoords="offset points", ha="center", fontsize=6)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=6.5)
    ax.set_ylabel("% of baseline q<0.05 hits\nretained at q<0.05 under LOO")
    ax.set_ylim(0, 110)
    ax.set_title("B", loc="left", fontweight="bold", fontsize=10, pad=6)
    ax.text(0.5, 1.04, "threshold stability scales with sample size",
            transform=ax.transAxes, ha="center", fontsize=8)
    ax.legend(loc="upper left", frameon=False, fontsize=7)


def main():
    fig = plt.figure(figsize=(10.0, 3.6))
    gs = GridSpec(1, 2, width_ratios=[1.0, 1.1], wspace=0.30)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    fig.savefig(OUT / "fig3_robustness.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig3_robustness.png", bbox_inches="tight", dpi=300)
    print(f"wrote {OUT/'fig3_robustness.pdf'} + .png")


if __name__ == "__main__":
    main()
