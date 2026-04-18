#!/usr/bin/env python3
"""Subject-bootstrap uncertainty at the protein level.

Two modes:

- **Paired delta mode** (default). For each protein in a paired contrast
  A − B, resample subjects with replacement (B reps) and recompute the
  mean per-subject delta. Mirrors `bootstrap_module_de.py`. This is the
  mode used by the Dube and Gisby cohorts in the main text.

- **OLS regression mode** (when ``--covariates`` is supplied). Per-sample
  design `y ~ group + covariates`; resample samples with replacement; for
  each resample, refit the OLS model per protein and record β_group. The
  bootstrap distribution of β_group across resamples gives the non-
  parametric uncertainty on the covariate-adjusted group effect. Used for
  the SIH case-control cohort where the primary contrast is Leak vs
  NoLeak with age+sex covariate adjustment.

Output schema is the same in both modes:

  gene_symbol, n_subjects, point_mean, boot_mean, boot_sd,
  ci_lo_95, ci_hi_95, p_pos, p_sign_stable, contrast, n_boot,
  mode, covariates

``mode`` is ``delta`` or ``ols``; ``covariates`` is a comma-separated list
of the covariate columns applied (empty in delta mode). Downstream
readers should interpret ``point_mean`` as "mean subject-delta" in
delta mode and as "β_group coefficient" in OLS mode.

P(sign stable) is the primary per-feature inferential quantity in both
modes; q-values remain available as calibration and SSR as a legacy
LOO diagnostic.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ----------------------------------------------------------------------
# Delta mode (paired): mirror of module bootstrap
# ----------------------------------------------------------------------
def load_delta(input_dir: Path, contrast: str) -> pd.DataFrame:
    """Return wide per-subject delta matrix (subject × protein)."""
    meas = pd.read_csv(input_dir / "qc_measurements.tsv", sep="\t")
    if "dropped_by_qc" in meas.columns:
        meas = meas[meas["dropped_by_qc"] == 0]
    samp = pd.read_csv(input_dir / "samples.tsv", sep="\t")
    cond_a, cond_b = contrast.split("-")
    sub = meas.merge(samp[["sample_id", "subject_id", "condition"]],
                     on="sample_id")
    sub = sub[sub["condition"].isin([cond_a, cond_b])]
    pivot = sub.pivot_table(
        index=["subject_id", "gene_symbol"], columns="condition",
        values="abundance", aggfunc="mean",
    ).dropna(subset=[cond_a, cond_b])
    pivot["delta"] = pivot[cond_a] - pivot[cond_b]
    pivot = pivot.reset_index()
    wide = pivot.pivot(index="subject_id", columns="gene_symbol",
                       values="delta").dropna(axis=1, how="any")
    return wide


def bootstrap_delta(wide: pd.DataFrame, n_boot: int,
                    rng: np.random.Generator) -> pd.DataFrame:
    n = wide.shape[0]
    X = wide.to_numpy()  # (n_subj, n_proteins)
    draws = rng.integers(0, n, size=(n_boot, n))
    boot_means = X[draws].mean(axis=1)  # (n_boot, n_proteins)
    rows = {
        "gene_symbol":      wide.columns.to_numpy(),
        "n_subjects":       np.full(len(wide.columns), n),
        "point_mean":       X.mean(axis=0),
        "boot_mean":        boot_means.mean(axis=0),
        "boot_sd":          boot_means.std(axis=0, ddof=1),
        "ci_lo_95":         np.quantile(boot_means, 0.025, axis=0),
        "ci_hi_95":         np.quantile(boot_means, 0.975, axis=0),
        "p_pos":            (boot_means > 0).mean(axis=0),
    }
    df = pd.DataFrame(rows)
    df["p_sign_stable"] = df[["p_pos"]].apply(
        lambda r: max(r["p_pos"], 1 - r["p_pos"]), axis=1
    )
    return df


# ----------------------------------------------------------------------
# OLS mode (case-control, covariate-adjusted)
# ----------------------------------------------------------------------
def _encode_covariates(samp: pd.DataFrame,
                       covariates: list[str]) -> tuple[np.ndarray, list[str]]:
    cols: list[np.ndarray] = []
    labels: list[str] = []
    for cov in covariates:
        vals = samp[cov].astype(str).values
        # Numeric iff every value parses as float.
        try:
            nums = np.array([float(v) for v in vals])
            cols.append(nums.reshape(-1, 1))
            labels.append(cov)
        except ValueError:
            levels = sorted(set(vals))
            for lvl in levels[1:]:
                cols.append((vals == lvl).astype(float).reshape(-1, 1))
                labels.append(f"{cov}={lvl}")
    if not cols:
        return np.zeros((len(samp), 0)), []
    return np.concatenate(cols, axis=1), labels


def load_ols(input_dir: Path, contrast: str,
             covariates: list[str]) -> tuple:
    """Return (Y, gene_symbols, group, cov_mat, cov_labels, sample_ids).

    ``Y`` is (n_samples, n_proteins) with NaN where a protein is missing
    in a given sample; per-protein fits mask those out. ``group`` is 1
    for comp_a and 0 for comp_b, matching atman's OLS encoding convention.
    """
    meas = pd.read_csv(input_dir / "qc_measurements.tsv", sep="\t")
    if "dropped_by_qc" in meas.columns:
        meas = meas[meas["dropped_by_qc"] == 0]
    samp = pd.read_csv(input_dir / "samples.tsv", sep="\t")
    cond_a, cond_b = contrast.split("-")
    samp = samp[samp["condition"].isin([cond_a, cond_b])]
    samp = samp[samp["is_control"] == 0]
    for cov in covariates:
        if cov not in samp.columns:
            raise SystemExit(f"covariate {cov!r} not found in samples.tsv")
        # Drop rows where the covariate is missing (empty string or NaN).
        mask = samp[cov].notna() & (samp[cov].astype(str) != "")
        samp = samp[mask]
    sample_ids = samp["sample_id"].tolist()

    meas_sub = meas[meas["sample_id"].isin(sample_ids)]
    Y = meas_sub.pivot_table(index="sample_id", columns="gene_symbol",
                             values="abundance", aggfunc="mean")
    Y = Y.reindex(sample_ids)
    gene_symbols = Y.columns.tolist()

    group = (samp["condition"].values == cond_a).astype(float)
    cov_mat, cov_labels = _encode_covariates(samp, covariates)
    return (
        Y.to_numpy(dtype=float),
        gene_symbols,
        group,
        cov_mat,
        cov_labels,
        sample_ids,
    )


def _fit_beta_group(X: np.ndarray, y: np.ndarray,
                    group_col: int, min_n: int) -> float:
    """OLS β for column ``group_col``. Returns NaN if n < min_n, the
    subdesign is rank-deficient, or any numerical issue arises."""
    mask = ~np.isnan(y)
    n = int(mask.sum())
    p = X.shape[1]
    if n < max(min_n, p + 1):
        return np.nan
    Xm = X[mask]
    ym = y[mask]
    # Full-rank guard: if a categorical level vanished from the subsample,
    # skip rather than return a silently-biased estimate.
    try:
        rank = np.linalg.matrix_rank(Xm)
    except np.linalg.LinAlgError:
        return np.nan
    if rank < p:
        return np.nan
    try:
        beta, *_ = np.linalg.lstsq(Xm, ym, rcond=None)
    except np.linalg.LinAlgError:
        return np.nan
    return float(beta[group_col])


def bootstrap_ols(Y: np.ndarray, group: np.ndarray, cov_mat: np.ndarray,
                  gene_symbols: list[str], n_boot: int,
                  rng: np.random.Generator, min_n: int = 5) -> pd.DataFrame:
    n_samples, n_proteins = Y.shape
    X_full = np.concatenate([
        np.ones((n_samples, 1)),
        group.reshape(-1, 1),
        cov_mat,
    ], axis=1)
    group_col = 1

    # Point estimates.
    point = np.full(n_proteins, np.nan)
    for j in range(n_proteins):
        point[j] = _fit_beta_group(X_full, Y[:, j], group_col, min_n)

    # Bootstrap. Skip the inner fit for proteins whose point estimate
    # failed the min_n filter on the full sample set — there is no
    # "original direction" to report sign-stability against, and running
    # the bootstrap anyway would report stats computed on a biased subset
    # of resamples.
    valid = ~np.isnan(point)
    boot = np.full((n_boot, n_proteins), np.nan)
    for b in range(n_boot):
        idx = rng.integers(0, n_samples, size=n_samples)
        Xb = X_full[idx]
        Yb = Y[idx]
        for j in range(n_proteins):
            if not valid[j]:
                continue
            boot[b, j] = _fit_beta_group(Xb, Yb[:, j], group_col, min_n)

    # Per-protein effective n (samples with a non-missing abundance).
    eff_n = (~np.isnan(Y)).sum(axis=0)

    rows = {
        "gene_symbol": np.array(gene_symbols),
        "n_subjects":  eff_n,                # same semantic as delta mode
        "point_mean":  point,                # β_group point estimate
        "boot_mean":   np.nanmean(boot, axis=0),
        "boot_sd":     np.nanstd(boot, axis=0, ddof=1),
        "ci_lo_95":    np.nanquantile(boot, 0.025, axis=0),
        "ci_hi_95":    np.nanquantile(boot, 0.975, axis=0),
        "p_pos":       np.nanmean(boot > 0, axis=0),
    }
    df = pd.DataFrame(rows)
    df["p_sign_stable"] = df["p_pos"].apply(lambda p: max(p, 1 - p))
    return df


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument(
        "--covariates", type=str, default="",
        help="Comma-separated covariate column names from samples.tsv. "
             "When supplied, switches to OLS bootstrap mode (per-sample "
             "design y ~ group + covariates, resample samples, refit "
             "per resample, record β_group).",
    )
    ap.add_argument(
        "--min-n", type=int, default=5,
        help="Minimum effective sample count per protein before a fit is "
             "attempted.",
    )
    args = ap.parse_args()

    covariates = [c.strip() for c in args.covariates.split(",") if c.strip()]
    rng = np.random.default_rng(args.seed)

    if covariates:
        Y, genes, group, cov_mat, cov_labels, sids = load_ols(
            args.input_dir, args.contrast, covariates,
        )
        if len(sids) < max(args.min_n, 2 + cov_mat.shape[1] + 1):
            print(f"[bootstrap] too few samples for OLS "
                  f"(n={len(sids)}, covariate cols={len(cov_labels)})",
                  file=sys.stderr)
            return 1
        out = bootstrap_ols(Y, group, cov_mat, genes, args.n_boot, rng,
                            min_n=args.min_n)
        out["mode"] = "ols"
        out["covariates"] = ",".join(covariates)
        print(f"[bootstrap] OLS mode: {len(sids)} samples, "
              f"{cov_mat.shape[1]} covariate cols ({','.join(cov_labels)})",
              file=sys.stderr)
    else:
        wide = load_delta(args.input_dir, args.contrast)
        if wide.shape[0] < 4 or wide.shape[1] == 0:
            print(f"[bootstrap] too few subjects/proteins "
                  f"(n={wide.shape[0]}, p={wide.shape[1]})",
                  file=sys.stderr)
            return 1
        out = bootstrap_delta(wide, args.n_boot, rng)
        out["mode"] = "delta"
        out["covariates"] = ""

    out["contrast"] = args.contrast
    out["n_boot"] = args.n_boot
    out = out.sort_values("p_sign_stable", ascending=False).reset_index(drop=True)
    out.to_csv(args.output, sep="\t", index=False)
    n_1 = int((out["p_sign_stable"] >= 0.999).sum())
    n_95 = int((out["p_sign_stable"] >= 0.95).sum())
    print(f"wrote {args.output}  ({len(out)} proteins: "
          f"{n_1} with P(sign)>=0.999, {n_95} with P(sign)>=0.95)",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
