#!/usr/bin/env python3
"""Per-protein OLS regression of physiological heat-response endpoints on
delta-NPX, using Dube's deposited dNPX files and Physiological_data.xlsx.

This is what Dube's own analysis.py does (with a much larger configuration.yaml
of endpoints, many of which combine phenotypic covariates); we reproduce the
subset that's directly interpretable: for each (endpoint, protein) pair,
    endpoint ~ protein_delta
OLS across 10 subjects, df = 8. Report β, SE, t, p, BH-FDR within endpoint.

Endpoints (per-subject physiological deltas):
    d1_dTcore  = Tcore(hyperthermic, day1) - Tcore(normothermic, day1)
    d1_dTskin  = Tskin(hyperthermic, day1) - Tskin(normothermic, day1)
    d1_dTbody  = Tbody(hyperthermic, day1) - Tbody(normothermic, day1)
    d1_dLSR    = LSR(hyperthermic, day1)   - LSR(normothermic, day1)
    d1_dLDF    = LDF(hyperthermic, day1)   - LDF(normothermic, day1)
    d1_dCVC    = CVC(hyperthermic, day1)   - CVC(normothermic, day1)
    d1_dSSNA   = SSNA(hyperthermic, day1)  - SSNA(normothermic, day1)
    d7_d*      = same at day 7
    acc_d*     = day7 normothermic - day1 normothermic (acclimation at rest)

Predictors (per-subject protein deltas from Dube's dNPX files):
    hd1_<assay>      for d1_* endpoints   (dNPX_hd1)
    hd7_<assay>      for d7_* endpoints   (dNPX_hd7)
    accpre_<assay>   for acc_d*           (dNPX_accpre)

Usage:
    scripts/phenotype_regression.py <output_dir>

Reads:
    example_data/dube_heat_2023/Physiological_data.xlsx
    example_data/dube_heat_2023/delta_npx/dNPX_hd1.csv
    example_data/dube_heat_2023/delta_npx/dNPX_hd7.csv
    example_data/dube_heat_2023/delta_npx/dNPX_accpre.csv
"""

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.stats import t as student_t

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = REPO_ROOT / "example_data" / "dube_heat_2023"
PHENO_XLSX = DATA_ROOT / "Physiological_data.xlsx"
DELTA_DIR = DATA_ROOT / "delta_npx"

# Participants in canonical (alphabetical) order, matching Dube's dNPX column order.
PARTICIPANTS = ["001B", "007", "008", "009", "011", "014", "015", "016", "018", "019"]

# Variable column indices in Physiological_data.xlsx sheet "Heat exposure".
# Row layout (verified by inspection):
#   row 1:       "Pre-acclimation visit" / "Post-acclimation visit"
#   row 2:       "Normothermic" / "Hyperthermic" section headers
#   row 3:       column labels (Participant, Tcore initial, ..., Tcore final, ..., Change, Heat exposure, Sweat rate)
#   rows 4-13:   D1 (pre-acclimation) participant data
#   row 14:      "Post-acclimation visit" header
#   row 15:      (spacer)
#   row 16:      column labels (same as row 3)
#   rows 17-26:  D7 (post-acclimation) participant data
#
# Initial (normothermic) columns for 11 phenotypes, cols B:L (indices 1..11):
#   B=Tcore, C=Tskin, D=Tbody, E=LSR, F=LDF, G=CVC, H=HR, I=SBP, J=DBP, K=MAP, L=SSNA
# Final (hyperthermic) columns, cols N:X (indices 13..23):
#   N=Tcore_f, O=Tskin_f, P=Tbody_f, Q=LSR_f, R=LDF_f, S=CVC_f, T=HR_f, U=SBP_f, V=DBP_f, W=MAP_f, X=SSNA_f

PHENO_VARS = ["Tcore", "Tskin", "Tbody", "LSR", "LDF", "CVC", "HR", "SBP", "DBP", "MAP", "SSNA"]
INIT_COL_START = 1  # column B (0-indexed)
FINAL_COL_START = 13  # column N (0-indexed)

# Participant row ranges (0-indexed after header).
D1_PARTICIPANT_ROWS = slice(3, 13)   # rows 4..13 (10 participants)
D7_PARTICIPANT_ROWS = slice(16, 26)  # rows 17..26 (10 participants)


def load_phenotypes() -> Dict[str, Dict[str, float]]:
    """Return {participant_id: {var: value}} for 4 states:
    var keys like 'd1_PR_Tcore', 'd1_PT_Tcore', 'd7_PR_Tcore', 'd7_PT_Tcore'."""
    import openpyxl

    wb = openpyxl.load_workbook(PHENO_XLSX, data_only=True)
    ws = wb["Heat exposure"]
    rows = list(ws.iter_rows(values_only=True))

    def extract(row_slice: slice, day_tag: str) -> Dict[str, Dict[str, float]]:
        out: Dict[str, Dict[str, float]] = {}
        for row in rows[row_slice]:
            pid_raw = row[0]
            if pid_raw is None:
                continue
            pid = str(pid_raw).replace("-", "")  # "001-B" -> "001B"
            state: Dict[str, float] = {}
            for i, var in enumerate(PHENO_VARS):
                init_val = row[INIT_COL_START + i]
                final_val = row[FINAL_COL_START + i]
                # Some variables (HR, BP, SSNA) may be missing for some participants.
                state[f"{day_tag}_PR_{var}"] = _to_float(init_val)
                state[f"{day_tag}_PT_{var}"] = _to_float(final_val)
            out[pid] = state
        return out

    d1 = extract(D1_PARTICIPANT_ROWS, "d1")
    d7 = extract(D7_PARTICIPANT_ROWS, "d7")

    pheno: Dict[str, Dict[str, float]] = {}
    for pid in PARTICIPANTS:
        state: Dict[str, float] = {}
        if pid in d1:
            state.update(d1[pid])
        if pid in d7:
            state.update(d7[pid])
        pheno[pid] = state
    return pheno


def _to_float(v) -> float:
    if v is None or v == "":
        return float("nan")
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def build_endpoints(pheno: Dict[str, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
    """Return {endpoint_name: {participant: value}} where each endpoint is a
    per-subject physiological delta. Endpoints are scoped so they match the
    corresponding dNPX file:

        d1_dX   = X(d1_PT) - X(d1_PR)       pairs with dNPX_hd1
        d7_dX   = X(d7_PT) - X(d7_PR)       pairs with dNPX_hd7
        acc_dX  = X(d7_PR) - X(d1_PR)       pairs with dNPX_accpre
    """
    endpoints: Dict[str, Dict[str, float]] = {}
    for var in PHENO_VARS:
        d1 = {}
        d7 = {}
        acc = {}
        for pid in PARTICIPANTS:
            s = pheno[pid]
            d1[pid] = s.get(f"d1_PT_{var}", float("nan")) - s.get(f"d1_PR_{var}", float("nan"))
            d7[pid] = s.get(f"d7_PT_{var}", float("nan")) - s.get(f"d7_PR_{var}", float("nan"))
            acc[pid] = s.get(f"d7_PR_{var}", float("nan")) - s.get(f"d1_PR_{var}", float("nan"))
        endpoints[f"d1_d{var}"] = d1
        endpoints[f"d7_d{var}"] = d7
        endpoints[f"acc_d{var}"] = acc
    return endpoints


def load_dnpx(path: Path) -> Dict[Tuple[str, str, str], Dict[str, float]]:
    """Return {(olink_id, assay, panel): {participant: delta}}."""
    out: Dict[Tuple[str, str, str], Dict[str, float]] = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (row["OlinkID"], row["Assay"], row["Panel"])
            vals: Dict[str, float] = {}
            for pid in PARTICIPANTS:
                raw = row.get(pid, "")
                if raw == "":
                    continue
                try:
                    vals[pid] = float(raw)
                except ValueError:
                    continue
            out[key] = vals
    return out


def ols_slope(x: np.ndarray, y: np.ndarray) -> Tuple[int, float, float, float, float, float]:
    """Simple OLS with intercept: y = α + β·x.
    Returns (n, beta, se, t, p_two_sided, r_squared). n includes only pairs
    where both x and y are finite.
    """
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 3:
        return (n, math.nan, math.nan, math.nan, math.nan, math.nan)
    xm = x[mask]
    ym = y[mask]
    x_mean = xm.mean()
    y_mean = ym.mean()
    xd = xm - x_mean
    yd = ym - y_mean
    sxx = (xd * xd).sum()
    if sxx == 0.0:
        return (n, math.nan, math.nan, math.nan, math.nan, math.nan)
    sxy = (xd * yd).sum()
    beta = sxy / sxx
    alpha = y_mean - beta * x_mean
    y_hat = alpha + beta * xm
    resid = ym - y_hat
    rss = (resid * resid).sum()
    df = n - 2
    if df <= 0:
        return (n, beta, math.nan, math.nan, math.nan, math.nan)
    sigma2 = rss / df
    se_beta = math.sqrt(sigma2 / sxx) if sxx > 0 else math.nan
    if se_beta == 0.0 or math.isnan(se_beta):
        return (n, beta, se_beta, math.nan, math.nan, math.nan)
    t_stat = beta / se_beta
    p_two = 2.0 * (1.0 - student_t.cdf(abs(t_stat), df))
    tss = (yd * yd).sum()
    r_squared = 1.0 - rss / tss if tss > 0 else math.nan
    return (n, beta, se_beta, t_stat, float(p_two), r_squared)


def bh_fdr(p_values: List[float]) -> List[float]:
    m = sum(1 for p in p_values if not math.isnan(p))
    if m == 0:
        return [math.nan] * len(p_values)
    indexed = [(i, p) for i, p in enumerate(p_values) if not math.isnan(p)]
    indexed.sort(key=lambda x: x[1])
    q_sorted = [p * m / (rank + 1) for rank, (_, p) in enumerate(indexed)]
    q_sorted = [min(1.0, q) for q in q_sorted]
    for i in range(len(q_sorted) - 2, -1, -1):
        if q_sorted[i + 1] < q_sorted[i]:
            q_sorted[i] = q_sorted[i + 1]
    out = [math.nan] * len(p_values)
    for rank, (orig, _) in enumerate(indexed):
        out[orig] = q_sorted[rank]
    return out


# Map endpoint prefix to the matching dNPX file.
ENDPOINT_TO_DNPX = {
    "d1_": "hd1",
    "d7_": "hd7",
    "acc_": "accpre",
}


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: phenotype_regression.py <output_dir>", file=sys.stderr)
        return 2
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)

    pheno = load_phenotypes()
    endpoints = build_endpoints(pheno)
    print(f"loaded phenotypes for {len(pheno)} participants", file=sys.stderr)
    print(f"built {len(endpoints)} endpoint deltas", file=sys.stderr)

    dnpx_files = {
        "hd1": load_dnpx(DELTA_DIR / "dNPX_hd1.csv"),
        "hd7": load_dnpx(DELTA_DIR / "dNPX_hd7.csv"),
        "accpre": load_dnpx(DELTA_DIR / "dNPX_accpre.csv"),
    }
    print(
        "loaded dNPX files: "
        + ", ".join(f"{k}={len(v)} assays" for k, v in dnpx_files.items()),
        file=sys.stderr,
    )

    # Regression per (endpoint, protein).
    all_rows: List[dict] = []
    report_rows: List[dict] = []

    for ep_name, ep_vals in endpoints.items():
        dnpx_tag = next(
            (tag for prefix, tag in ENDPOINT_TO_DNPX.items() if ep_name.startswith(prefix)),
            None,
        )
        if dnpx_tag is None:
            continue
        dnpx = dnpx_files[dnpx_tag]
        y = np.array([ep_vals.get(p, math.nan) for p in PARTICIPANTS], dtype=float)
        if np.sum(np.isfinite(y)) < 5:
            continue

        family_rows: List[dict] = []
        family_p: List[float] = []
        for (olink_id, assay, panel), assay_vals in dnpx.items():
            x = np.array([assay_vals.get(p, math.nan) for p in PARTICIPANTS], dtype=float)
            n, beta, se, t_stat, p_val, r2 = ols_slope(x, y)
            row = {
                "endpoint": ep_name,
                "panel": panel,
                "assay_id": olink_id,
                "gene_symbol": assay,
                "dnpx_source": dnpx_tag,
                "n": n,
                "beta": beta,
                "se": se,
                "t": t_stat,
                "p_value": p_val,
                "r_squared": r2,
            }
            family_rows.append(row)
            family_p.append(p_val)

        qs = bh_fdr(family_p)
        for r, q in zip(family_rows, qs):
            r["bh_q"] = q
        all_rows.extend(family_rows)

        # Family summary.
        computed = [r for r in family_rows if not math.isnan(r["p_value"])]
        n_q05 = sum(1 for r in computed if not math.isnan(r["bh_q"]) and r["bh_q"] < 0.05)
        n_q10 = sum(1 for r in computed if not math.isnan(r["bh_q"]) and r["bh_q"] < 0.10)
        min_q = min((r["bh_q"] for r in computed if not math.isnan(r["bh_q"])), default=math.nan)
        max_abs_beta = max(
            (abs(r["beta"]) for r in computed if not math.isnan(r["beta"])),
            default=math.nan,
        )
        report_rows.append(
            {
                "endpoint": ep_name,
                "dnpx_source": dnpx_tag,
                "n_tests": len(family_rows),
                "n_computed": len(computed),
                "n_q_lt_05": n_q05,
                "n_q_lt_10": n_q10,
                "min_q": min_q,
                "max_abs_beta": max_abs_beta,
            }
        )

    # Sort: endpoint, then q ascending.
    def sort_key(r: dict):
        q = r.get("bh_q", math.nan)
        q_sort = q if not math.isnan(q) else math.inf
        return (r["endpoint"], q_sort, -(abs(r["beta"]) if not math.isnan(r["beta"]) else 0.0))

    all_rows.sort(key=sort_key)

    # Write outputs.
    results_path = out_dir / "phenotype_regression_results.tsv"
    with open(results_path, "w") as f:
        f.write(
            "endpoint\tpanel\tassay_id\tgene_symbol\tdnpx_source\t"
            "n\tbeta\tse\tt\tp_value\tr_squared\tbh_q\n"
        )
        for r in all_rows:
            f.write(
                f"{r['endpoint']}\t{r['panel']}\t{r['assay_id']}\t{r['gene_symbol']}\t"
                f"{r['dnpx_source']}\t{r['n']}\t"
                f"{_fmt(r['beta'])}\t{_fmt(r['se'])}\t{_fmt(r['t'])}\t"
                f"{_fmt(r['p_value'])}\t{_fmt(r['r_squared'])}\t{_fmt(r['bh_q'])}\n"
            )

    report_path = out_dir / "phenotype_regression_report.tsv"
    with open(report_path, "w") as f:
        f.write(
            "endpoint\tdnpx_source\tn_tests\tn_computed\tn_q_lt_05\tn_q_lt_10\tmin_q\tmax_abs_beta\n"
        )
        for r in report_rows:
            f.write(
                f"{r['endpoint']}\t{r['dnpx_source']}\t{r['n_tests']}\t"
                f"{r['n_computed']}\t{r['n_q_lt_05']}\t{r['n_q_lt_10']}\t"
                f"{_fmt(r['min_q'])}\t{_fmt(r['max_abs_beta'])}\n"
            )

    print(f"wrote {results_path}", file=sys.stderr)
    print(f"wrote {report_path}", file=sys.stderr)
    return 0


def _fmt(v: float) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    return f"{v}"


if __name__ == "__main__":
    sys.exit(main())
