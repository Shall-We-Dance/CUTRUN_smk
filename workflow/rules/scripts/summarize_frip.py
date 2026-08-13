#!/usr/bin/env python3
"""Convert per-sample vertical FRiP metrics into one rectangular table."""

import csv
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
