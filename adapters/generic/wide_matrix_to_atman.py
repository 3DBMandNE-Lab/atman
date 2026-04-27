#!/usr/bin/env python3
"""Convert a wide proteomics matrix to Atman's canonical TSV schema.

This adapter is intentionally generic. It handles the common release-facing
case where protein abundances are already normalized in a CSV/TSV matrix and
sample metadata is available in a separate CSV/TSV file.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd


CANONICAL_MEASUREMENT_COLUMNS = [
    "platform",
    "sample_id",
    "assay_id",
    "gene_symbol",
    "panel",
    "npx_source_str",
    "abundance",
    "abundance_raw",
    "abundance_unit",
    "qc_sample",
    "qc_assay",
    "detection_limit",
    "below_lod",
    "dropped_by_qc",
    "plate_id",
    "panel_lot",
    "ingest_order",
]


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix", required=True, type=Path)
    p.add_argument("--samples", required=True, type=Path)
    p.add_argument(
        "--proteins",
        type=Path,
        help="Protein metadata table. Required for --orientation samples-rows.",
    )
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--platform", required=True)
    p.add_argument("--abundance-unit", required=True)
    p.add_argument(
        "--orientation",
        choices=["proteins-rows", "samples-rows"],
        default="proteins-rows",
    )
    p.add_argument("--assay-id-col", required=True)
    p.add_argument("--gene-col")
    p.add_argument("--uniprot-col")
    p.add_argument("--panel", default="")
    p.add_argument("--panel-col")
    p.add_argument("--sample-id-col", default="sample_id")
    p.add_argument("--subject-id-col")
    p.add_argument("--condition-col", required=True)
    p.add_argument("--sample-type-col")
    p.add_argument("--is-control-col")
    p.add_argument("--log2-transform", action="store_true")
    return p.parse_args()


def sample_sheet(meta: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    sid = meta[args.sample_id_col].astype(str)
    out = pd.DataFrame(
        {
            "sample_id": sid,
            "subject_id": (
                meta[args.subject_id_col].astype(str) if args.subject_id_col else sid
            ),
            "condition": meta[args.condition_col].astype(str),
            "is_control": (
                meta[args.is_control_col].astype(int) if args.is_control_col else 0
            ),
            "sample_type": (
                meta[args.sample_type_col].astype(str) if args.sample_type_col else ""
            ),
            "ingest_order": range(1, len(meta) + 1),
        }
    )

    for col in meta.columns:
        if col not in out.columns and col not in {
            args.sample_id_col,
            args.subject_id_col,
            args.condition_col,
            args.sample_type_col,
            args.is_control_col,
        }:
            out[col] = meta[col]
    return out


def protein_sheet(matrix: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    gene = matrix[args.gene_col].astype(str) if args.gene_col else ""
    uniprot = matrix[args.uniprot_col].astype(str) if args.uniprot_col else ""
    panel = matrix[args.panel_col].astype(str) if args.panel_col else args.panel
    return pd.DataFrame(
        {
            "platform": args.platform,
            "assay_id": matrix[args.assay_id_col].astype(str),
            "uniprot": uniprot,
            "gene_symbol": gene,
            "panel": panel,
            "panel_lot": "",
        }
    )


def maybe_log2(v: object, log2_transform: bool) -> float | None:
    if pd.isna(v):
        return None
    x = float(v)
    if not math.isfinite(x):
        return None
    if log2_transform:
        if x <= 0:
            return None
        x = math.log2(x)
    return x


def measurements_proteins_rows(
    matrix: pd.DataFrame,
    samples: pd.DataFrame,
    proteins: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    meta_cols = {
        args.assay_id_col,
        args.gene_col,
        args.uniprot_col,
        args.panel_col,
    }
    meta_cols.discard(None)
    sample_ids = set(samples["sample_id"].astype(str))
    value_cols = [c for c in matrix.columns if c not in meta_cols and str(c) in sample_ids]
    if not value_cols:
        raise SystemExit("no matrix columns matched sample IDs from samples table")

    id_cols = [args.assay_id_col]
    for optional in [args.gene_col, args.panel_col]:
        if optional and optional not in id_cols:
            id_cols.append(optional)
    long = matrix.melt(
        id_vars=id_cols,
        value_vars=value_cols,
        var_name="sample_id",
        value_name="source_value",
    )
    long["assay_id"] = long[args.assay_id_col].astype(str)
    long["gene_symbol"] = long[args.gene_col].astype(str) if args.gene_col else ""
    long["panel"] = long[args.panel_col].astype(str) if args.panel_col else args.panel
    return finalize_measurements(long, args)


def measurements_samples_rows(
    matrix: pd.DataFrame,
    samples: pd.DataFrame,
    proteins: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    assay_ids = set(proteins["assay_id"].astype(str))
    value_cols = [c for c in matrix.columns if str(c) in assay_ids]
    if not value_cols:
        raise SystemExit("no matrix columns matched assay IDs from protein metadata")
    long = matrix.melt(
        id_vars=[args.sample_id_col],
        value_vars=value_cols,
        var_name="assay_id",
        value_name="source_value",
    ).rename(columns={args.sample_id_col: "sample_id"})
    meta = proteins[["assay_id", "gene_symbol", "panel"]]
    long = long.merge(meta, on="assay_id", how="left")
    return finalize_measurements(long, args)


def finalize_measurements(long: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    abundance = long["source_value"].map(lambda v: maybe_log2(v, args.log2_transform))
    out = pd.DataFrame(
        {
            "platform": args.platform,
            "sample_id": long["sample_id"].astype(str),
            "assay_id": long["assay_id"].astype(str),
            "gene_symbol": long.get("gene_symbol", "").fillna("").astype(str),
            "panel": long.get("panel", args.panel).fillna("").astype(str),
            "npx_source_str": long["source_value"].astype(str),
            "abundance": abundance,
            "abundance_raw": abundance,
            "abundance_unit": args.abundance_unit,
            "qc_sample": "PASS",
            "qc_assay": "PASS",
            "detection_limit": "",
            "below_lod": 0,
            "dropped_by_qc": abundance.isna().astype(int),
            "plate_id": "",
            "panel_lot": "",
            "ingest_order": range(1, len(long) + 1),
        }
    )
    out.loc[out["abundance"].isna(), ["abundance", "abundance_raw"]] = ""
    return out[CANONICAL_MEASUREMENT_COLUMNS]


def main() -> None:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    matrix = read_table(args.matrix)
    meta = read_table(args.samples)
    samples = sample_sheet(meta, args)

    if args.orientation == "proteins-rows":
        proteins = protein_sheet(matrix, args)
        measurements = measurements_proteins_rows(matrix, samples, proteins, args)
    else:
        if args.proteins is None:
            raise SystemExit("--orientation samples-rows requires --proteins")
        protein_meta = read_table(args.proteins)
        if args.assay_id_col not in protein_meta.columns:
            raise SystemExit(f"missing assay-id column {args.assay_id_col!r} in proteins table")
        proteins = protein_sheet(protein_meta, args)
        measurements = measurements_samples_rows(matrix, samples, proteins, args)

    samples.to_csv(args.out / "samples.tsv", sep="\t", index=False)
    proteins.to_csv(args.out / "proteins.tsv", sep="\t", index=False)
    measurements.to_csv(args.out / "measurements.tsv", sep="\t", index=False)

    print(
        f"wrote {args.out}: {len(samples)} samples, "
        f"{len(proteins)} proteins, {len(measurements)} measurements"
    )


if __name__ == "__main__":
    main()
