from pathlib import Path


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


def read_preseq_status_tables(paths):
    rows = {}
    for path in paths:
        with open(path, "r", encoding="utf-8") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            for line in handle:
                values = line.rstrip("\n").split("\t")
                row = dict(zip(header, values))
                sample = row.get("sample")
                if sample:
                    rows[sample] = row
    return rows


usable_paths = list(snakemake.input.usable)
fragment_paths = {Path(path).name.split(".fragment_lengths.tsv")[0]: path for path in snakemake.input.fragments}
preseq_paths = {Path(path).name.split(".preseq.lc_extrap.txt")[0]: path for path in snakemake.input.preseq}
preseq_status_rows = read_preseq_status_tables(snakemake.input.preseq_status)
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
    "preseq_lc_extrap",
    "preseq_status",
    "preseq_message",
    "spikein_reads",
    "spikein_fragments",
    "spikein_scale_factor",
    "spikein_status",
]

for path in usable_paths:
    metrics = read_metric_table(path)
    sample = metrics["sample"]
    preseq_status = preseq_status_rows.get(sample, {})
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
            "preseq_lc_extrap": preseq_paths.get(sample, ""),
            "preseq_status": preseq_status.get("preseq_status", ""),
            "preseq_message": preseq_status.get("preseq_message", ""),
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
