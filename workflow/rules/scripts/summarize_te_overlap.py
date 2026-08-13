#!/usr/bin/env python3
"""Summarize per-sample TE overlap metrics into one TSV + MultiQC JSON.

MultiQC 1.27 custom-content tables require ``data`` as a dict keyed by
sample name (each value a {header: value} dict); a list of rows makes the
module crash with AssertionError (and MultiQC still exits 0, silently
dropping the section).
"""

import json
from pathlib import Path

rows = []
for path in snakemake.input:
    metrics = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            metric, value = line.rstrip("\n").split("\t", 1)
            metrics[metric] = value
    rows.append(metrics)

columns = [
    "sample",
    "total_alignments",
    "secondary_alignments",
    "total_templates",
    "multimapped_templates",
    "fraction_templates_multimapped",
    "templates_in_te_unweighted",
    "fraction_templates_in_te_unweighted",
    "weighted_in_te",
    "fraction_weighted_in_te",
    "multi_overlap_alignments",
    "fragments_after_processing",
]

out_tsv = snakemake.output.tsv
Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
with open(out_tsv, "w") as out:
    out.write("\t".join(columns) + "\n")
    for row in sorted(rows, key=lambda r: r.get("sample", "")):
        out.write("\t".join(str(row.get(c, "")) for c in columns) + "\n")

# MultiQC table data: {sample: {header: value}} with numeric columns
# converted to float so the header format specifiers apply.
data = {}
for row in sorted(rows, key=lambda r: r.get("sample", "")):
    sample = row.get("sample", "")
    entry = {}
    for c in columns:
        if c == "sample":
            continue
        try:
            entry[c] = float(row.get(c, ""))
        except ValueError:
            entry[c] = row.get(c, "")
    data[sample] = entry

mqc_payload = {
    "id": "te_overlap",
    "section_name": "TE-mode overlap QC (fractional 1/N)",
    "description": (
        "Fraction of templates whose alignments overlap annotated TE loci, "
        "with 1/N fractional weighting for multimappers; unit = template "
        "(query name). Secondary-alignment-based overcounting is avoided."
    ),
    "plot_type": "table",
    "pconfig": {"id": "te_overlap_table", "title": "TE-mode overlap QC"},
    "headers": {
        "sample": {"title": "Sample", "description": "Sample name"},
        "fraction_weighted_in_te": {
            "title": "frac weighted in TE",
            "description": "Weighted (1/N) fraction of templates in TE loci",
            "format": "{:0.4f}",
        },
        "fraction_templates_in_te_unweighted": {
            "title": "frac templates in TE (unweighted)",
            "description": "Unweighted fraction of templates with >=1 TE-overlapping alignment",
            "format": "{:0.4f}",
        },
        "fraction_templates_multimapped": {
            "title": "frac templates multimapped",
            "description": "Templates with >=1 secondary alignment / total templates",
            "format": "{:0.4f}",
        },
        "weighted_in_te": {
            "title": "weighted templates in TE",
            "description": "Sum of 1/N weights over TE-overlapping alignments",
            "format": "{:,.1f}",
        },
        "total_templates": {
            "title": "total templates",
            "description": "Unique query names after TE-mode filtering",
            "format": "{:,.0f}",
        },
    },
    "data": data,
}

out_mqc = snakemake.output.mqc
Path(out_mqc).parent.mkdir(parents=True, exist_ok=True)
with open(out_mqc, "w") as out:
    json.dump(mqc_payload, out, indent=2)
    out.write("\n")
