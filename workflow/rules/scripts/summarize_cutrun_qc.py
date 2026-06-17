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
    points = []
    with open(path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header != ["fragment_length", "count"]:
            raise ValueError(f"{path} is not a fragment_length/count TSV.")
        for line in handle:
            if not line.strip():
                continue
            fragment_length, count = line.rstrip("\n").split("\t", 1)
            points.append([int(fragment_length), int(count)])
    return points


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
    "usable_reads_after_processing",
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
    rows.append(
        {
            "sample": sample,
            "read_type": metrics.get("read_type", ""),
            "total_reads_after_processing": metrics.get("total_reads_after_processing", ""),
            "mapped_reads_after_processing": metrics.get("mapped_reads_after_processing", ""),
            "usable_reads_after_processing": metrics.get("usable_reads_after_processing", ""),
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
    "description": "Fragment or read length distributions calculated from the final BAM after configured filtering steps.",
    "plot_type": "linegraph",
    "pconfig": {
        "id": "cutrun_fragment_lengths_plot",
        "title": "CUT&RUN fragment length distributions",
        "xlab": "Fragment / read length (bp)",
        "ylab": "Count",
        "xDecimals": False,
        "yDecimals": False,
    },
    "data": fragment_length_data,
}

Path(snakemake.output.mqc).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.output.mqc, "w", encoding="utf-8") as out:
    json.dump(multiqc_payload, out, indent=2)
    out.write("\n")
