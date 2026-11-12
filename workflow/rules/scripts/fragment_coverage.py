"""Exact one-fragment-per-overlapping-bin coverage from shared BED3 fragments.

CPM/RPKM denominators count fragments, not alignment ends. BPM uses total
fragment-bin incidences; RPGC uses sum(fragment lengths)/retained genome size.
Empty inputs produce explicit zero-valued tracks, not NaNs or missing outputs.
"""
import argparse
import csv
import heapq
import itertools
import math
from pathlib import Path

from fragment_records import bed_fragments, chrom_sizes


def fragment_statistics(path, sizes, bin_size):
    count = bases = incidences = 0
    previous = None
    for chrom, start, end in bed_fragments(path):
        key = (chrom, start)
        if previous is not None and key < previous:
            raise ValueError("Fragment BED must be chromosome/start sorted")
        if chrom not in sizes or end > sizes[chrom]:
            raise ValueError(f"Fragment outside retained genome: {chrom}:{start}-{end}")
        previous = key
        count += 1
        bases += end - start
        incidences += (end - 1) // bin_size - start // bin_size + 1
    return count, bases, incidences


def normalization_factor(mode, count, bases, incidences, genome_size, bin_size):
    if mode == "None":
        return 1.0
    if count == 0:
        return 0.0
    if mode == "CPM":
        return 1e6 / count
    if mode == "RPKM":
        return 1e9 / (count * bin_size)
    if mode == "BPM":
        return 1e6 / incidences
    if mode == "RPGC":
        return genome_size / bases
    raise ValueError(f"Unsupported normalization: {mode}")


def chromosome_runs(records, length, bin_size):
    """Streaming sweep: memory scales with local fragment depth, not genome size."""
    ends = []
    cursor = 0
    for start_bin, group in itertools.groupby(records, key=lambda r: r[1] // bin_size):
        while ends and ends[0] <= start_bin:
            boundary = ends[0]
            if boundary > cursor:
                yield cursor * bin_size, min(boundary * bin_size, length), len(ends)
            while ends and ends[0] == boundary:
                heapq.heappop(ends)
            cursor = boundary
        if start_bin > cursor:
            yield cursor * bin_size, min(start_bin * bin_size, length), len(ends)
        for _, start, end in group:
            heapq.heappush(ends, (end - 1) // bin_size + 1)
        cursor = start_bin
    while ends:
        boundary = ends[0]
        if boundary > cursor:
            yield cursor * bin_size, min(boundary * bin_size, length), len(ends)
        while ends and ends[0] == boundary:
            heapq.heappop(ends)
        cursor = boundary
    if cursor * bin_size < length:
        yield cursor * bin_size, length, 0


def write_bigwig(bed, sizes, output, bin_size=10, normalization="CPM", scale_factor=1.0):
    import pyBigWig
    if bin_size < 1 or not math.isfinite(scale_factor) or scale_factor < 0:
        raise ValueError("bin_size must be positive and scale_factor finite/nonnegative")
    count, bases, incidences = fragment_statistics(bed, sizes, bin_size)
    scale = normalization_factor(normalization, count, bases, incidences, sum(sizes.values()), bin_size) * scale_factor
    groups = iter(itertools.groupby(bed_fragments(bed), key=lambda r: r[0]))
    current = next(groups, None)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with pyBigWig.open(str(output), "w") as out:
        out.addHeader(sorted(sizes.items()))
        for chrom, length in sorted(sizes.items()):
            records = current[1] if current is not None and current[0] == chrom else ()
            starts, ends, values = [], [], []
            for start, end, value in chromosome_runs(records, length, bin_size):
                starts.append(start)
                ends.append(end)
                values.append(float(value * scale))
                if len(starts) == 10000:
                    out.addEntries([chrom] * len(starts), starts, ends=ends, values=values)
                    starts, ends, values = [], [], []
            if starts:
                out.addEntries([chrom] * len(starts), starts, ends=ends, values=values)
            if current is not None and current[0] == chrom:
                current = next(groups, None)
    return count, scale


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fragments", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--remove-chrm", action="store_true")
    parser.add_argument("--out", required=True)
    parser.add_argument("--bin-size", type=int, default=10)
    parser.add_argument("--normalization", choices=["None", "CPM", "RPKM", "BPM", "RPGC"], default="CPM")
    parser.add_argument("--scale-factor", type=float, default=1.0)
    parser.add_argument("--scale-table", help="Spike-in TSV; requires exactly one data row")
    args = parser.parse_args()
    if args.scale_table:
        with open(args.scale_table) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if len(rows) != 1:
            parser.error("--scale-table must contain exactly one data row")
        args.scale_factor *= float(rows[0]["scale_factor"])
        print(f"spikein_status={rows[0].get('status', 'unknown')}")
    count, scale = write_bigwig(args.fragments, chrom_sizes(args.chrom_sizes, args.remove_chrm),
                               args.out, args.bin_size, args.normalization, args.scale_factor)
    print(f"unit=fragment fragments={count} normalization={args.normalization} scale_factor={scale}")


if __name__ == "__main__":
    main()
