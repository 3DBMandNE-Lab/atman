#!/usr/bin/env python3
"""Fig 1 — method overview.

Three panels in the Nature Methods schematic + concept + data pattern:

  a. Iconographic workflow. Five small icons — raw NPX, protein DE,
     module assembly, subject bootstrap, P(sign stable) — with thin
     arrows. No CLI names, no band washes. The command DAG lives in
     the supplement.

  b. Scale hierarchy. A three-tier concept diagram stating the core
     claim of the paper: protein-level effects are aggregated into
     module-level effects, which are resampled over subjects to yield a
     distribution whose concordance gives P(sign stable), the primary
     inferential quantity.

  c. Headline data. Per-protein −log10(q) vs P(sign stable) on the
     Gisby late-vs-early contrast. The two are monotone-aligned, but
     P(sign stable) saturates earlier and resolves cases that a
     q-only threshold collapses.

Strict _style.py discipline — no rounded corners, no background washes,
0.5 pt strokes, Helvetica 7 pt.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

from _style import (
    apply_style, save, panel_letter, clean_axes,
    PRIMARY, ACCENT, SECOND, LEGACY, MUTED, STRUCTURE, REPO,
)

apply_style()


# ==================================================================
# panel a — iconographic workflow
# ==================================================================
def _draw_matrix(ax, cx, cy, w, h):
    rng = np.random.default_rng(2)
    M = rng.normal(size=(6, 4))
    vmin, vmax = -2.5, 2.5
    nrows, ncols = M.shape
    cw, ch = w / ncols, h / nrows
    cmap = plt.get_cmap("Greys")
    for i in range(nrows):
        for j in range(ncols):
            v = (M[i, j] - vmin) / (vmax - vmin)
            v = max(0.0, min(1.0, v))
            ax.add_patch(Rectangle(
                (cx - w/2 + j * cw, cy + h/2 - (i + 1) * ch),
                cw, ch,
                facecolor=cmap(v * 0.85 + 0.05), edgecolor="none",
                zorder=3,
            ))
    ax.add_patch(Rectangle(
        (cx - w/2, cy - h/2), w, h,
        facecolor="none", edgecolor=STRUCTURE, lw=0.5, zorder=4,
    ))


def _draw_de(ax, cx, cy, w, h):
    x = np.linspace(-4, 4, 200)
    a = np.exp(-0.5 * (x + 0.6) ** 2)
    b = np.exp(-0.5 * (x - 0.6) ** 2)
    xmap = cx - w/2 + (x + 4) / 8 * w
    amap = cy - h/2 + a * h * 0.92
    bmap = cy - h/2 + b * h * 0.92
    overlap = cy - h/2 + np.minimum(a, b) * h * 0.92
    ax.fill_between(xmap, cy - h/2, overlap, color=MUTED, alpha=0.35, lw=0,
                    zorder=2)
    ax.plot(xmap, amap, color=LEGACY, lw=0.9, zorder=3)
    ax.plot(xmap, bmap, color=PRIMARY, lw=0.9, zorder=3)


def _draw_modules(ax, cx, cy, w, h):
    rng = np.random.default_rng(7)
    centers = [(-1.0, 0.55), (0.4, -0.65), (1.15, 0.75)]
    cols = [PRIMARY, ACCENT, SECOND]
    xs, ys, cs = [], [], []
    for (ux, uy), c in zip(centers, cols):
        pts = rng.normal(size=(8, 2)) * 0.22
        for p in pts:
            xs.append(ux + p[0]); ys.append(uy + p[1]); cs.append(c)
    # map (-2..2, -1.5..1.5) into (cx±w/2, cy±h/2)
    xs = cx + np.array(xs) / 2.0 * (w / 2)
    ys = cy + np.array(ys) / 1.5 * (h / 2)
    ax.scatter(xs, ys, s=5, c=cs, lw=0, alpha=0.85, zorder=3)


def _draw_bootstrap_psign(ax, cx, cy, w, h):
    """Combined bootstrap densities + P(sign stable) readout."""
    x = np.linspace(-2.5, 2.5, 200)
    rng = np.random.default_rng(3)
    # densities occupy the top ~70% of the icon box
    top_y0 = cy - h/2 + h * 0.28
    top_h  = h * 0.66
    xmap0 = cx - w/2 + (0 + 2.5) / 5 * w
    ax.plot([xmap0, xmap0], [top_y0, top_y0 + top_h * 0.95],
            color=MUTED, lw=0.5, ls="--", zorder=2)
    for i in range(5):
        mu = rng.normal(0.55, 0.12)
        y = np.exp(-0.5 * ((x - mu) / 0.6) ** 2)
        xmap = cx - w/2 + (x + 2.5) / 5 * w
        ymap = top_y0 + y * top_h * 0.75 + i * 0.035 * top_h
        ax.plot(xmap, ymap, color=ACCENT, lw=0.55, alpha=0.8, zorder=3)

    # P(sign stable) bar in the lower band
    p = 0.93
    bar_h = h * 0.18
    bar_y = cy - h/2 + h * 0.04
    ax.add_patch(Rectangle(
        (cx - w/2, bar_y), w, bar_h,
        facecolor="#F1F5F9", edgecolor=STRUCTURE, lw=0.5, zorder=3,
    ))
    ax.add_patch(Rectangle(
        (cx - w/2, bar_y), w * p, bar_h,
        facecolor=PRIMARY, edgecolor="none", zorder=4,
    ))
    ax.text(cx - w/2 + w * p / 2, bar_y + bar_h / 2, f"P = {p:.2f}",
            ha="center", va="center", fontsize=6.2, color="white",
            fontweight="bold", zorder=5)


def panel_a(ax):
    """Single-axes iconographic workflow, 2x2 U-flow."""
    ax.set_xlim(0, 10); ax.set_ylim(0, 6)
    ax.set_aspect("auto")
    ax.axis("off")

    # 2x2 slot centres; U-flow: TL → TR → BR → BL
    positions = [
        (2.5, 4.6),  # 0: raw NPX       (top-left)
        (7.5, 4.6),  # 1: protein DE    (top-right)
        (7.5, 1.6),  # 2: modules       (bottom-right)
        (2.5, 1.6),  # 3: bootstrap     (bottom-left)
    ]
    icon_w, icon_h = 2.0, 1.9

    _draw_matrix          (ax, *positions[0], icon_w, icon_h)
    _draw_de              (ax, *positions[1], icon_w, icon_h)
    _draw_modules         (ax, *positions[2], icon_w, icon_h)
    _draw_bootstrap_psign (ax, *positions[3], icon_w, icon_h)

    # arrows forming the U-flow
    gap_x = icon_w / 2 + 0.20
    gap_y = icon_h / 2 + 0.40
    # TL → TR (horizontal, top row)
    ax.add_patch(FancyArrowPatch(
        (positions[0][0] + gap_x, positions[0][1]),
        (positions[1][0] - gap_x, positions[1][1]),
        arrowstyle="-|>,head_width=3.2,head_length=5",
        color=STRUCTURE, lw=0.7, shrinkA=0, shrinkB=0, zorder=5,
    ))
    # TR → BR (vertical, right column)
    ax.add_patch(FancyArrowPatch(
        (positions[1][0], positions[1][1] - gap_y),
        (positions[2][0], positions[2][1] + gap_y),
        arrowstyle="-|>,head_width=3.2,head_length=5",
        color=STRUCTURE, lw=0.7, shrinkA=0, shrinkB=0, zorder=5,
    ))
    # BR → BL (horizontal, bottom row, leftward)
    ax.add_patch(FancyArrowPatch(
        (positions[2][0] - gap_x, positions[2][1]),
        (positions[3][0] + gap_x, positions[3][1]),
        arrowstyle="-|>,head_width=3.2,head_length=5",
        color=STRUCTURE, lw=0.7, shrinkA=0, shrinkB=0, zorder=5,
    ))

    # labels placed on the outside of each icon so they never cross an
    # arrow: top-row labels go above, bottom-row labels go below.
    labels = {
        0: ("raw NPX",                                   "top"),
        1: ("protein DE",                                "top"),
        2: ("modules",                                   "bottom"),
        3: (r"bootstrap $\to$ $P(\mathrm{sign\ stable})$", "bottom"),
    }
    for idx, (px, py) in enumerate(positions):
        text, side = labels[idx]
        if side == "top":
            ly = py + icon_h/2 + 0.20
            va = "bottom"
        else:
            ly = py - icon_h/2 - 0.28
            va = "top"
        ax.text(px, ly, text,
                ha="center", va=va, fontsize=7, color=STRUCTURE)

    ax.text(-0.01, 1.02, "a", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=9, fontweight="bold",
            color=STRUCTURE)


# ==================================================================
# panel b — scale hierarchy concept diagram
# ==================================================================
def panel_b(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis("off")

    tiers = [
        (0.86, "protein scale",           r"$\Delta_i,\ i = 1..P$",                            PRIMARY),
        (0.53, "module scale",            r"$\bar d_m = \langle \Delta_i \rangle_{i\in m}$",   ACCENT),
        (0.20, "subject bootstrap",       r"$\{\bar d_m^{(b)}\}_{b=1}^{B}$",                   SECOND),
    ]
    x0, x1 = 0.05, 0.80
    for y, name, expr, col in tiers:
        ax.add_patch(Rectangle((x0, y - 0.10), x1 - x0, 0.20,
                               facecolor="white",
                               edgecolor=col, lw=0.6, zorder=2))
        ax.text((x0 + x1) / 2, y + 0.04, name,
                ha="center", va="center", fontsize=7.2, color=col)
        ax.text((x0 + x1) / 2, y - 0.055, expr,
                ha="center", va="center", fontsize=6.8, color=STRUCTURE)

    mid_x = (x0 + x1) / 2
    # arrow protein → module
    ax.add_patch(FancyArrowPatch(
        (mid_x, 0.76), (mid_x, 0.63),
        arrowstyle="-|>,head_width=2.8,head_length=4.4",
        color=STRUCTURE, lw=0.6, shrinkA=0, shrinkB=0))
    # arrow module → bootstrap
    ax.add_patch(FancyArrowPatch(
        (mid_x, 0.43), (mid_x, 0.30),
        arrowstyle="-|>,head_width=2.8,head_length=4.4",
        color=STRUCTURE, lw=0.6, shrinkA=0, shrinkB=0))

    # output label below the last tier
    ax.text(mid_x, 0.06, r"$P(\mathrm{sign\ stable})$",
            ha="center", va="center", fontsize=7.5, color=STRUCTURE,
            fontweight="bold")
    ax.add_patch(FancyArrowPatch(
        (mid_x, 0.10), (mid_x, 0.00) ,  # not drawn, placeholder
        arrowstyle="-", color="none"))
    # thin line segment from bottom tier to output label
    ax.add_patch(FancyArrowPatch(
        (mid_x, 0.10), (mid_x, 0.13),
        arrowstyle="-|>,head_width=2.8,head_length=4.4",
        color=STRUCTURE, lw=0.6, shrinkA=0, shrinkB=0))

    panel_letter(ax, "b", x=-0.04, y=1.02)


# ==================================================================
# panel c — headline data: -log10(q) vs P(sign stable), Gisby proteins
# ==================================================================
def panel_c(ax):
    de   = pd.read_csv(REPO / "out_gisby/de_results.tsv", sep="\t")
    boot = pd.read_csv(REPO / "out_gisby/robustness/protein_bootstrap.tsv", sep="\t")
    de = de[(de.comparison == "late-early") & de.skip_reason.isna()].copy()
    df = de.merge(boot[["gene_symbol", "p_sign_stable"]],
                  on="gene_symbol", how="inner")
    df["nlq"] = -np.log10(df["bh_q"].clip(lower=1e-30))

    rng = np.random.default_rng(0)
    x = df["p_sign_stable"].to_numpy() + rng.uniform(-0.003, 0.003, len(df))

    ax.scatter(x, df["nlq"],
               s=4.0, color=PRIMARY, lw=0, alpha=0.55, rasterized=True)

    q_cut = -np.log10(0.05)
    ax.axhline(q_cut, color=MUTED, lw=0.5, ls="--", zorder=0)
    ax.axvline(0.95,  color=MUTED, lw=0.5, ls="--", zorder=0)
    ax.text(0.03, q_cut + 0.25, r"$q = 0.05$",
            fontsize=6.3, color="#4B5563")
    ax.text(0.935, 0.3, r"$P = 0.95$",
            fontsize=6.3, color="#4B5563", rotation=90, va="bottom", ha="right")

    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, max(8, df["nlq"].max() * 1.04))
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_box_aspect(1.0)
    clean_axes(ax,
               xlabel=r"$P(\mathrm{sign\ stable})$",
               ylabel=r"$-\log_{10}\ q$")
    panel_letter(ax, "c", x=-0.22)


# ==================================================================
def main():
    fig = plt.figure(figsize=(7.09, 3.2))

    # explicit positions [left, bottom, width, height] in figure fraction
    #                                 (panel a is now 2x2, roughly square;
    #                                 it shrinks horizontally)
    ax_a = fig.add_axes([0.02, 0.06, 0.32, 0.90]); panel_a(ax_a)
    ax_b = fig.add_axes([0.40, 0.08, 0.26, 0.86]); panel_b(ax_b)
    ax_c = fig.add_axes([0.75, 0.18, 0.22, 0.72]); panel_c(ax_c)

    save(fig, "fig1_pipeline")
    plt.close(fig)


if __name__ == "__main__":
    main()
