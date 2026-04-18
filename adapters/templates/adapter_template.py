#!/usr/bin/env python3
"""Template for writing Atman canonical TSV adapters.

Copy this file when a source format needs custom parsing beyond
`adapters/generic/wide_matrix_to_atman.py`.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


OUT = Path("out_adapter")
OUT.mkdir(parents=True, exist_ok=True)

# 1. Load source files.
# matrix = ...
# metadata = ...

# 2. Write samples.tsv. Extra covariate columns are allowed.
samples = pd.DataFrame(
    columns=[
        "sample_id",
        "subject_id",
        "condition",
        "is_control",
        "sample_type",
        "ingest_order",
    ]
)

# 3. Write proteins.tsv.
proteins = pd.DataFrame(
    columns=["platform", "assay_id", "uniprot", "gene_symbol", "panel", "panel_lot"]
)

# 4. Write measurements.tsv and qc_measurements.tsv.
measurements = pd.DataFrame(
    columns=[
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
)

samples.to_csv(OUT / "samples.tsv", sep="\t", index=False)
proteins.to_csv(OUT / "proteins.tsv", sep="\t", index=False)
measurements.to_csv(OUT / "measurements.tsv", sep="\t", index=False)
measurements.to_csv(OUT / "qc_measurements.tsv", sep="\t", index=False)
