#!/usr/bin/env python3
"""Convert per-sample vertical FRiP metrics into one rectangular table."""

import csv
import json
from pathlib import Path


METRICS = [
    "total_reads",
    "mapped_reads",
    "usable_templates",
    "templates_in_peaks",
    "fragment_frip",
    "reads_in_peaks",
    "read_frip",
]


def read_metrics(path):
    values = {}
    with open(path, newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0] == "metric" or len(row) < 2:
                continue
            values[row[0]] = row[1]
    missing = [metric for metric in METRICS if metric not in values]
    if missing:
        raise ValueError(f"{path}: missing FRiP metrics: {', '.join(missing)}")
    return values


paths = [str(path) for path in list(snakemake.input.per_rep) + list(snakemake.input.pooled)]
Path(str(snakemake.output.tsv)).parent.mkdir(parents=True, exist_ok=True)

rows = []
with open(str(snakemake.output.tsv), "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["sample", "level", *METRICS])
    for path in paths:
        metrics = read_metrics(path)
        name = Path(path).name
        pooled = name.endswith(".pooled.frip.tsv")
        suffix = ".pooled.frip.tsv" if pooled else ".frip.tsv"
        sample = name[: -len(suffix)]
        writer.writerow([sample, "pooled" if pooled else "replicate", *[metrics[m] for m in METRICS]])
        rows.append((sample, "pooled" if pooled else "replicate", metrics))

data = {}
for sample, level, metrics in rows:
    row_id = f"{sample} (pooled)" if level == "pooled" else sample
    data[row_id] = {
        "level": level,
        "usable_templates": float(metrics["usable_templates"]),
        "templates_in_peaks": float(metrics["templates_in_peaks"]),
        "fragment_frip": float(metrics["fragment_frip"]),
        "read_frip": float(metrics["read_frip"]),
    }

payload = {
    "id": "cutrun_frip",
    "section_name": "CUT&RUN FRiP",
    "description": (
        "Fraction of usable fragments overlapping the corresponding MACS3 peak set. "
        "Pooled rows are shown only for samples with at least two biological replicates."
    ),
    "plot_type": "table",
    "pconfig": {"id": "frip_table", "title": "CUT&RUN FRiP"},
    "headers": {
        "level": {"title": "Level"},
        "usable_templates": {
            "title": "Usable fragments",
            "format": "{:,.0f}",
        },
        "templates_in_peaks": {
            "title": "Fragments in peaks",
            "format": "{:,.0f}",
        },
        "fragment_frip": {
            "title": "Fragment FRiP",
            "description": "Usable fragments overlapping peaks / usable fragments",
            "format": "{:.3f}",
        },
        "read_frip": {
            "title": "Read FRiP",
            "description": "Reads overlapping peaks / mapped reads",
            "format": "{:.3f}",
        },
    },
    "data": data,
}

with open(str(snakemake.output.mqc), "w") as handle:
    json.dump(payload, handle, indent=2)
    handle.write("\n")
