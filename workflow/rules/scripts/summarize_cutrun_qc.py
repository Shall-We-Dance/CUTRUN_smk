from pathlib import Path
import json


def read_metric_table(path):
    metrics = {}
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header != ["metric", "value"]:
            raise ValueError(f"{path} is not a metric/value TSV.")
        for line in handle:
            metric, value = line.rstrip("\n").split("\t", 1)
            metrics[metric] = value
    return metrics


def read_spikein_table(path):
    if not path:
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = {}
        for line in handle:
            values = line.rstrip("\n").split("\t")
            row = dict(zip(header, values))
            sample = row.get("sample")
            if sample:
                rows[sample] = row
    return rows


def read_fragment_lengths(path):
    counts = {}
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header != ["fragment_length", "count"]:
            raise ValueError(f"{path} is not a fragment_length/count TSV.")
        for line in handle:
            if not line.strip():
                continue
            fragment_length, count = line.rstrip("\n").split("\t", 1)
            length = str(int(fragment_length))
            counts[length] = counts.get(length, 0) + int(count)
    # JSON object keys are strings. Some MultiQC versions sort these before
    # numeric conversion (100, 101, ..., 50), creating backtracking lines.
    # Equal-width numeric keys sort correctly in both string and numeric paths;
    # the plot still uses a continuous numeric axis, not categorical spacing.
    width = max((len(length) for length in counts), default=1)
    return {length.zfill(width): counts[length] for length in sorted(counts, key=int)}


usable_paths = list(snakemake.input.usable)
fragment_paths = {Path(path).name.split(".fragment_lengths.tsv")[0]: path for path in snakemake.input.fragments}
spikein_paths = list(snakemake.input.spikein)
spikein_rows = read_spikein_table(spikein_paths[0]) if spikein_paths else {}

rows = []
columns = [
    "sample",
    "read_type",
    "total_reads_after_processing",
    "mapped_reads_after_processing",
    "mapped_fraction_after_processing",
    "usable_reads_after_processing",
    "usable_fraction_after_processing",
    "proper_pair_reads_after_processing",
    "usable_fragments_after_processing",
    "fragment_lengths",
    "spikein_reads",
    "spikein_fragments",
    "spikein_scale_factor",
    "spikein_status",
]

for path in usable_paths:
    metrics = read_metric_table(path)
    sample = metrics["sample"]
    spikein = spikein_rows.get(sample, {})
    total_reads = int(metrics.get("total_reads_after_processing", 0) or 0)
    mapped_reads = int(metrics.get("mapped_reads_after_processing", 0) or 0)
    usable_reads = int(metrics.get("usable_reads_after_processing", 0) or 0)
    rows.append(
        {
            "sample": sample,
            "read_type": metrics.get("read_type", ""),
            "total_reads_after_processing": metrics.get("total_reads_after_processing", ""),
            "mapped_reads_after_processing": metrics.get("mapped_reads_after_processing", ""),
            "mapped_fraction_after_processing": mapped_reads / total_reads if total_reads else 0,
            "usable_reads_after_processing": metrics.get("usable_reads_after_processing", ""),
            "usable_fraction_after_processing": usable_reads / total_reads if total_reads else 0,
            "proper_pair_reads_after_processing": metrics.get("proper_pair_reads_after_processing", ""),
            "usable_fragments_after_processing": metrics.get("usable_fragments_after_processing", ""),
            "fragment_lengths": fragment_paths.get(sample, ""),
            "spikein_reads": spikein.get("spikein_reads", ""),
            "spikein_fragments": spikein.get("spikein_fragments", ""),
            "spikein_scale_factor": spikein.get("scale_factor", ""),
            "spikein_status": spikein.get("status", ""),
        }
    )

Path(snakemake.output.tsv).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.output.tsv, "w", encoding="utf-8") as out:
    out.write("\t".join(columns) + "\n")
    for row in sorted(rows, key=lambda item: item["sample"]):
        out.write("\t".join(str(row.get(column, "")) for column in columns) + "\n")

fragment_length_data = {
    sample: read_fragment_lengths(path)
    for sample, path in sorted(fragment_paths.items())
}
multiqc_payload = {
    "id": "cutrun_fragment_lengths",
    "section_name": "CUT&RUN Fragment Lengths",
    "description": "Accepted fragment-span lengths after shared pair, flag, chromosome and blacklist validation; SE lengths reflect the configured extension.",
    "plot_type": "linegraph",
    "pconfig": {
        "id": "cutrun_fragment_lengths_plot",
        "title": "CUT&RUN fragment length distributions",
        "xlab": "Fragment / read length (bp)",
        "ylab": "Count",
        "xDecimals": False,
        "yDecimals": False,
        "categories": False,
        "smooth_points": None,
    },
    "data": fragment_length_data,
}

Path(snakemake.output.fragments_mqc).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.output.fragments_mqc, "w", encoding="utf-8") as out:
    json.dump(multiqc_payload, out, indent=2)
    out.write("\n")

summary_data = {}
for row in sorted(rows, key=lambda item: item["sample"]):
    entry = {
        "read_type": row["read_type"],
        "total_reads": int(row["total_reads_after_processing"]),
        "mapped_fraction": float(row["mapped_fraction_after_processing"]),
        "usable_reads": int(row["usable_reads_after_processing"]),
        "usable_fraction": float(row["usable_fraction_after_processing"]),
        "usable_fragments": int(row["usable_fragments_after_processing"]),
    }
    if row["spikein_fragments"] != "":
        entry["spikein_fragments"] = int(row["spikein_fragments"])
        entry["spikein_scale_factor"] = float(row["spikein_scale_factor"])
        entry["spikein_status"] = row["spikein_status"]
    summary_data[row["sample"]] = entry

summary_payload = {
    "id": "cutrun_qc_summary",
    "section_name": "CUT&RUN processing QC",
    "description": (
        "Alignment counts from the final BAM after configured mapping, "
        "blacklist, and duplicate filtering; usable counts additionally require shared fragment validation. Fractions use total reads in that final BAM."
    ),
    "plot_type": "table",
    "pconfig": {"id": "cutrun_qc_table", "title": "CUT&RUN processing QC"},
    "headers": {
        "read_type": {"title": "Type", "description": "PE or SE"},
        "total_reads": {
            "title": "Final reads",
            "description": "Reads in the final processed BAM",
            "format": "{:,.0f}",
        },
        "mapped_fraction": {
            "title": "Mapped %",
            "description": "Mapped reads / final BAM reads",
            "format": "{:.1%}",
        },
        "usable_reads": {
            "title": "Usable reads",
            "description": "Reads passing configured SAM-flag exclusions",
            "format": "{:,.0f}",
        },
        "usable_fraction": {
            "title": "Usable %",
            "description": "Usable reads / final BAM reads",
            "format": "{:.1%}",
        },
        "usable_fragments": {
            "title": "Usable fragments",
            "description": "Validated two-mate templates for PE; accepted extended alignments for SE",
            "format": "{:,.0f}",
        },
        "spikein_fragments": {
            "title": "Spike-in fragments",
            "description": "Mapped spike-in fragments used for scaling",
            "format": "{:,.0f}",
        },
        "spikein_scale_factor": {
            "title": "Spike-in scale",
            "description": "Configured spike-in target / spike-in fragments",
            "format": "{:.4g}",
        },
        "spikein_status": {"title": "Spike-in status"},
    },
    "data": summary_data,
}

with open(snakemake.output.summary_mqc, "w", encoding="utf-8") as out:
    json.dump(summary_payload, out, indent=2)
    out.write("\n")
