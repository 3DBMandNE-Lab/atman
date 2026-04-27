#!/usr/bin/env python3
"""Convert Olink Explore NGS long-CSV NPX exports to Atman's canonical TSVs.

Atman has no built-in proteomics ingest. All upstream formats are converted
to the canonical schema (samples.tsv / proteins.tsv / measurements.tsv) by
adapters; the Rust binary only operates on the canonical schema. This
adapter handles the Olink Explore NGS NPX CSV format used by the Dube et al.
2023 reproduction fixture and similar datasets.

Sample IDs follow the Dube convention by default:

    SSNA-<subject>-<PR1|PR2|PT1|PT2>    biological samples
    CONTROL_SAMPLE_*                    technical controls

Datasets with a different sample ID convention need their own parser; pass
`--sample-id-regex` and `--control-prefix` to override.

Usage:

    python3 adapters/generic/olink_explore_to_atman.py \\
        --output-dir out \\
        path/to/file_one.csv path/to/file_two.csv

The output directory then contains samples.tsv, proteins.tsv, and
measurements.tsv ready for `atman qc`, `atman matrix`, etc.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path

EXPECTED_HEADERS = [
    "SampleID",
    "Index",
    "OlinkID",
    "UniProt",
    "Assay",
    "MissingFreq",
    "Panel",
    "Panel_Lot_Nr",
    "PlateID",
    "QC_Warning",
    "LOD",
    "NPX",
    "Normalization",
    "Assay_Warning",
]

DEFAULT_SAMPLE_ID_REGEX = r"^SSNA-(?P<subject>[^-]+)-(?P<condition>PR1|PR2|PT1|PT2)$"
DEFAULT_CONTROL_PREFIX = "CONTROL_SAMPLE_"

PLATFORM = "olink_explore_ngs"
ABUNDANCE_UNIT = "log2_npx"

MEASUREMENT_HEADER = [
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
SAMPLE_HEADER = [
    "sample_id",
    "subject_id",
    "condition",
    "is_control",
    "sample_type",
    "ingest_order",
]
PROTEIN_HEADER = [
    "platform",
    "assay_id",
    "uniprot",
    "gene_symbol",
    "panel",
    "panel_lot",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("inputs", nargs="+", type=Path, help="One or more NPX CSV files.")
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument(
        "--sample-id-regex",
        default=DEFAULT_SAMPLE_ID_REGEX,
        help=(
            "Regex matching biological sample IDs. Must define named groups "
            "`subject` and `condition`. Default matches the Dube SSNA-X-PR/PT "
            "convention."
        ),
    )
    p.add_argument(
        "--control-prefix",
        default=DEFAULT_CONTROL_PREFIX,
        help="Sample IDs starting with this prefix are technical controls.",
    )
    return p.parse_args()


def format_abundance(npx: float) -> str:
    """Mimic Rust's default `{}` Display for f64.

    Rust strips the trailing `.0` from integer-valued floats (so `0.0` → `"0"`,
    `-1.0` → `"-1"`) and preserves the IEEE-754 negative-zero sign (`-0.0` →
    `"-0"`). Python's `repr` keeps the `.0`. Match Rust so downstream byte-
    exact reproduction tests remain stable."""
    if math.isnan(npx):
        return "NaN"
    if npx == 0.0:
        return "-0" if math.copysign(1.0, npx) < 0 else "0"
    if not math.isinf(npx) and npx == int(npx):
        return str(int(npx))
    return repr(npx)


def parse_sample_id(
    sample_id: str, regex: re.Pattern, control_prefix: str
) -> tuple[str | None, str | None, bool]:
    """Returns (subject_id, condition, is_control)."""
    m = regex.match(sample_id)
    if m is not None:
        return m.group("subject"), m.group("condition"), False
    if sample_id.startswith(control_prefix):
        return None, None, True
    raise SystemExit(f"unparseable sample_id: {sample_id!r}")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    regex = re.compile(args.sample_id_regex)

    measurements: list[list[str]] = []
    samples: dict[str, dict[str, str]] = {}
    proteins: dict[str, dict[str, str]] = {}
    seen_pk: set[tuple[str, str]] = set()
    ingest_counter = 0

    for path in args.inputs:
        if not path.is_file():
            raise SystemExit(f"missing input file: {path}")
        with path.open("r", newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter=";")
            headers = reader.fieldnames or []
            missing = set(EXPECTED_HEADERS) - set(headers)
            extra = set(headers) - set(EXPECTED_HEADERS)
            if missing or extra:
                raise SystemExit(
                    f"header mismatch in {path}: missing={sorted(missing)} extra={sorted(extra)}"
                )

            for line_idx, row in enumerate(reader, start=2):
                if row["Normalization"] != "Plate control":
                    raise SystemExit(
                        f"{path}:{line_idx}: unexpected Normalization "
                        f"{row['Normalization']!r}; only 'Plate control' supported"
                    )

                npx_str = row["NPX"]
                lod_str = row["LOD"]
                try:
                    npx = float(npx_str)
                    lod = float(lod_str)
                except ValueError as exc:
                    raise SystemExit(f"{path}:{line_idx}: invalid abundance: {exc}")

                olink_id = row["OlinkID"]
                sample_id = row["SampleID"]

                if olink_id not in proteins:
                    proteins[olink_id] = {
                        "platform": PLATFORM,
                        "assay_id": olink_id,
                        "uniprot": row["UniProt"],
                        "gene_symbol": row["Assay"],
                        "panel": row["Panel"],
                        "panel_lot": row["Panel_Lot_Nr"],
                    }

                pk = (olink_id, sample_id)
                if pk in seen_pk:
                    raise SystemExit(
                        f"{path}:{line_idx}: duplicate (assay_id={olink_id}, "
                        f"sample_id={sample_id})"
                    )
                seen_pk.add(pk)

                if sample_id not in samples:
                    subject, condition, is_control = parse_sample_id(
                        sample_id, regex, args.control_prefix
                    )
                    samples[sample_id] = {
                        "sample_id": sample_id,
                        "subject_id": subject or "",
                        "condition": condition or "",
                        "is_control": "1" if is_control else "0",
                        "sample_type": "",
                        "ingest_order": str(ingest_counter),
                    }
                    ingest_counter += 1

                qc_sample = row["QC_Warning"] if row["QC_Warning"] in {"PASS", "WARN", "FAIL"} else "WARN"
                qc_assay = row["Assay_Warning"] if row["Assay_Warning"] in {"PASS", "WARN", "FAIL"} else "WARN"
                below_lod = "1" if npx < lod else "0"

                measurements.append(
                    [
                        PLATFORM,
                        sample_id,
                        olink_id,
                        row["Assay"],
                        row["Panel"],
                        npx_str,
                        format_abundance(npx),
                        format_abundance(npx),
                        ABUNDANCE_UNIT,
                        qc_sample,
                        qc_assay,
                        format_abundance(lod),
                        below_lod,
                        "0",
                        row["PlateID"],
                        row["Panel_Lot_Nr"],
                        str(ingest_counter),
                    ]
                )
                ingest_counter += 1

    # Write samples.tsv (insertion order = first-seen).
    samples_path = args.output_dir / "samples.tsv"
    with samples_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(SAMPLE_HEADER)
        for row in samples.values():
            w.writerow([row[c] for c in SAMPLE_HEADER])

    # Write proteins.tsv (sorted by assay_id, matching Rust ingest behavior).
    proteins_path = args.output_dir / "proteins.tsv"
    with proteins_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(PROTEIN_HEADER)
        for assay_id in sorted(proteins.keys()):
            row = proteins[assay_id]
            w.writerow([row[c] for c in PROTEIN_HEADER])

    # Write measurements.tsv (CSV-row order, matches Rust ingest).
    measurements_path = args.output_dir / "measurements.tsv"
    with measurements_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(MEASUREMENT_HEADER)
        w.writerows(measurements)

    print(
        f"olink-explore: {len(measurements)} rows, {len(proteins)} assays, "
        f"{len(samples)} samples, {len(args.inputs)} inputs",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
