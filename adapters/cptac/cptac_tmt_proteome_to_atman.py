#!/usr/bin/env python3
"""Convert a CPTAC TMT protein-level proteome TSV to Atman canonical TSVs.

CPTAC publishes TMT proteome quantifications as wide TSVs with this shape:

    Gene   <sample_a> Log Ratio   <sample_a> Unshared Log Ratio   ...

The first three rows hold per-sample summary statistics (Mean, Median,
StdDev) emitted by the upstream pipeline; they are skipped. Each
biological sample contributes two adjacent columns (`Log Ratio` and
`Unshared Log Ratio`); this adapter keeps `Log Ratio` (the standard
shared+razor quantification) and drops `Unshared Log Ratio`. Non-
biological columns (TumorOnlyIR, NormalOnlyIR, QC*, Pool*, Reference,
RefMix*) are dropped.

The condition assigned to each sample depends on the cohort. Some CPTAC
cohorts encode condition directly in the sample-ID prefix (e.g. CPTAC
HCC uses `T<N>` for tumor and `P<N>` for paired non-tumor); pass
`--sample-id-regex` with a `condition_key` named group plus
`--condition-key-map` to recover those labels from the column names. For
cohorts whose sample IDs are opaque (e.g. CPTAC GBM CPT-prefixed
aliquots, BRCA UUIDs), omit the regex; every biological sample is
assigned `--default-condition` (default `tumor`).

Output:
    <output-dir>/samples.tsv
    <output-dir>/proteins.tsv
    <output-dir>/measurements.tsv

Usage:

    # CPTAC HCC: T<N> / P<N> condition recovery from sample IDs
    python3 adapters/cptac/cptac_tmt_proteome_to_atman.py \\
        --proteome HCC_proteome.tsv \\
        --tumor-tag HCC \\
        --output-dir out/HCC \\
        --sample-id-regex '^(?P<condition_key>[TP])(?P<subject>\\d+)$' \\
        --condition-key-map 'T=tumor,P=paired_non_tumor'

    # CPTAC GBM: opaque CPT IDs, single-condition tumor cohort
    python3 adapters/cptac/cptac_tmt_proteome_to_atman.py \\
        --proteome GBM_proteome.tsv \\
        --tumor-tag GBM \\
        --output-dir out/GBM
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Iterable

DEFAULT_SUMMARY_ROW_LABELS = ("Mean", "Median", "StdDev")
DEFAULT_LOG_RATIO_SUFFIX = " Log Ratio"
DEFAULT_DROP_SUFFIX = " Unshared Log Ratio"
DEFAULT_NON_BIO_PATTERN = (
    r"^(?:TumorOnlyIR|NormalOnlyIR|QC\d*|Pool\d*|Reference|RefMix\d*)$"
)
DEFAULT_DEFAULT_CONDITION = "tumor"

PLATFORM = "cptac_tmt_proteome"
ABUNDANCE_UNIT = "log2_tmt_ratio"

SAMPLES_HEADER = (
    "sample_id",
    "subject_id",
    "condition",
    "is_control",
    "sample_type",
    "ingest_order",
)
PROTEINS_HEADER = (
    "platform",
    "assay_id",
    "uniprot",
    "gene_symbol",
    "panel",
    "panel_lot",
)
MEASUREMENTS_HEADER = (
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
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert a CPTAC TMT proteome TSV to Atman canonical TSVs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "--proteome",
        required=True,
        type=Path,
        help="Path to the CPTAC TMT proteome TSV (e.g. HCC_proteome.tsv).",
    )
    p.add_argument(
        "--tumor-tag",
        required=True,
        help="Short cohort identifier (e.g. HCC, GBM, BRCA). Becomes the `panel` value.",
    )
    p.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory to write samples.tsv, proteins.tsv, measurements.tsv.",
    )
    p.add_argument(
        "--gene-column",
        default="Gene",
        help="Column holding the gene symbol (default: Gene).",
    )
    p.add_argument(
        "--ncbi-gene-id-column",
        default="NCBIGeneID",
        help="Column holding NCBI gene ID, used as `assay_id` (default: NCBIGeneID; falls back to gene symbol if absent).",
    )
    p.add_argument(
        "--summary-row-labels",
        default=",".join(DEFAULT_SUMMARY_ROW_LABELS),
        help="Comma-separated row labels (in the gene column) to skip as summary stats.",
    )
    p.add_argument(
        "--log-ratio-suffix",
        default=DEFAULT_LOG_RATIO_SUFFIX,
        help="Column suffix marking the kept abundance column (default: ' Log Ratio').",
    )
    p.add_argument(
        "--drop-suffix",
        default=DEFAULT_DROP_SUFFIX,
        help="Column suffix to drop (default: ' Unshared Log Ratio').",
    )
    p.add_argument(
        "--non-bio-pattern",
        default=DEFAULT_NON_BIO_PATTERN,
        help="Regex matched against sample IDs (after stripping the log-ratio suffix); matching samples are dropped.",
    )
    p.add_argument(
        "--sample-id-regex",
        default=None,
        help="Optional regex with named groups `subject` (always) and optionally `condition_key` (looked up in --condition-key-map). When omitted, subject_id == sample_id and condition == --default-condition.",
    )
    p.add_argument(
        "--condition-key-map",
        default=None,
        help="Comma-separated KEY=value pairs mapping the regex-captured `condition_key` to a condition label (e.g. 'T=tumor,P=paired_non_tumor').",
    )
    p.add_argument(
        "--default-condition",
        default=DEFAULT_DEFAULT_CONDITION,
        help="Condition to assign when no `condition_key` is recovered (default: tumor).",
    )
    p.add_argument(
        "--abundance-unit",
        default=ABUNDANCE_UNIT,
        help=f"Value written to abundance_unit (default: {ABUNDANCE_UNIT}).",
    )
    p.add_argument(
        "--platform",
        default=PLATFORM,
        help=f"Value written to platform (default: {PLATFORM}).",
    )
    p.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress per-cohort progress messages on stderr.",
    )
    return p.parse_args(argv)


def parse_condition_key_map(spec: str | None) -> dict[str, str]:
    if spec is None:
        return {}
    out: dict[str, str] = {}
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if "=" not in token:
            raise ValueError(
                f"--condition-key-map entry {token!r} must be KEY=value"
            )
        key, value = token.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key or not value:
            raise ValueError(
                f"--condition-key-map entry {token!r} has empty key or value"
            )
        out[key] = value
    return out


def stream_proteome(path: Path) -> Iterable[list[str]]:
    """Yield raw rows from the proteome TSV. Header is the first row."""
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            yield row


def select_columns(
    header: list[str],
    log_ratio_suffix: str,
    drop_suffix: str,
    non_bio_pattern: str,
) -> tuple[list[tuple[int, str]], list[str]]:
    """Return [(column_index, sample_id), ...] for biological Log Ratio columns,
    plus the list of dropped sample IDs (for logging)."""
    drop_pat = re.compile(non_bio_pattern)
    kept: list[tuple[int, str]] = []
    dropped: list[str] = []
    for i, col in enumerate(header):
        if not col.endswith(log_ratio_suffix):
            continue
        # Discriminate the kept suffix from any wider drop suffix
        # (e.g. " Log Ratio" matches both " Log Ratio" and " Unshared Log Ratio"
        #  when log_ratio_suffix is " Log Ratio"; explicit drop_suffix wins).
        if drop_suffix and col.endswith(drop_suffix):
            continue
        sample_id = col[: -len(log_ratio_suffix)]
        if drop_pat.match(sample_id):
            dropped.append(sample_id)
            continue
        kept.append((i, sample_id))
    return kept, dropped


def assign_condition_subject(
    sample_id: str,
    sample_id_regex: re.Pattern[str] | None,
    condition_key_map: dict[str, str],
    default_condition: str,
) -> tuple[str, str]:
    """Return (condition, subject_id) for a sample ID."""
    if sample_id_regex is None:
        return default_condition, sample_id
    m = sample_id_regex.match(sample_id)
    if m is None:
        return default_condition, sample_id
    subject = m.groupdict().get("subject", sample_id)
    key = m.groupdict().get("condition_key")
    if key is not None and condition_key_map:
        condition = condition_key_map.get(key, default_condition)
    else:
        condition = default_condition
    return condition, subject


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.proteome.is_file():
        print(f"error: proteome file not found: {args.proteome}", file=sys.stderr)
        return 2
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary_labels = {s.strip() for s in args.summary_row_labels.split(",") if s.strip()}
    sample_id_regex = re.compile(args.sample_id_regex) if args.sample_id_regex else None
    if sample_id_regex is not None and "subject" not in sample_id_regex.groupindex:
        print(
            "error: --sample-id-regex must define a named group `subject`",
            file=sys.stderr,
        )
        return 2
    condition_key_map = parse_condition_key_map(args.condition_key_map)

    rows = stream_proteome(args.proteome)
    try:
        header = next(rows)
    except StopIteration:
        print(f"error: empty proteome file: {args.proteome}", file=sys.stderr)
        return 2
    if args.gene_column not in header:
        print(
            f"error: gene column {args.gene_column!r} not in proteome header",
            file=sys.stderr,
        )
        return 2
    gene_idx = header.index(args.gene_column)
    ncbi_idx = header.index(args.ncbi_gene_id_column) if args.ncbi_gene_id_column in header else None

    kept_cols, dropped_cols = select_columns(
        header,
        args.log_ratio_suffix,
        args.drop_suffix,
        args.non_bio_pattern,
    )
    if not kept_cols:
        print(
            f"error: no kept biological Log Ratio columns in {args.proteome}",
            file=sys.stderr,
        )
        return 2
    sample_ids = [sid for _, sid in kept_cols]
    if not args.quiet:
        print(
            f"[{args.tumor_tag}] kept {len(kept_cols)} biological samples, "
            f"dropped {len(dropped_cols)} non-bio columns: {dropped_cols!r}",
            file=sys.stderr,
        )

    samples_path = args.output_dir / "samples.tsv"
    proteins_path = args.output_dir / "proteins.tsv"
    measurements_path = args.output_dir / "measurements.tsv"

    with samples_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(SAMPLES_HEADER)
        for order, sample_id in enumerate(sample_ids, start=1):
            condition, subject_id = assign_condition_subject(
                sample_id,
                sample_id_regex,
                condition_key_map,
                args.default_condition,
            )
            w.writerow([sample_id, subject_id, condition, 0, "tissue", order])

    panel = f"CPTAC_{args.tumor_tag}"
    proteins_seen: set[str] = set()
    measurement_order = 0
    with proteins_path.open("w", encoding="utf-8", newline="") as fp, \
         measurements_path.open("w", encoding="utf-8", newline="") as fm:
        wp = csv.writer(fp, delimiter="\t", lineterminator="\n")
        wm = csv.writer(fm, delimiter="\t", lineterminator="\n")
        wp.writerow(PROTEINS_HEADER)
        wm.writerow(MEASUREMENTS_HEADER)
        for row in rows:
            if not row:
                continue
            gene_symbol = row[gene_idx].strip()
            if gene_symbol in summary_labels or not gene_symbol:
                continue
            assay_id = (
                row[ncbi_idx].strip() if ncbi_idx is not None and len(row) > ncbi_idx
                else gene_symbol
            )
            if not assay_id:
                assay_id = gene_symbol
            if assay_id not in proteins_seen:
                wp.writerow([args.platform, assay_id, "", gene_symbol, panel, ""])
                proteins_seen.add(assay_id)
            for col_idx, sample_id in kept_cols:
                if col_idx >= len(row):
                    continue
                cell = row[col_idx].strip()
                if not cell:
                    continue
                try:
                    value = float(cell)
                except ValueError:
                    continue
                measurement_order += 1
                wm.writerow(
                    [
                        args.platform,
                        sample_id,
                        assay_id,
                        gene_symbol,
                        panel,
                        cell,
                        f"{value:.10g}",
                        f"{value:.10g}",
                        args.abundance_unit,
                        "PASS",
                        "PASS",
                        "",
                        0,
                        0,
                        "",
                        "",
                        measurement_order,
                    ]
                )

    if not args.quiet:
        print(
            f"[{args.tumor_tag}] wrote {samples_path}, {proteins_path}, "
            f"{measurements_path} ({measurement_order} measurements, {len(proteins_seen)} proteins)",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
