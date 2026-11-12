"""Optional pairwise IDR on relaxed narrow peaks; not the full ENCODE QC workflow."""
import argparse
import csv
import itertools
import math
from pathlib import Path
import subprocess

from consensus import PeakIndex, read_peaks


def filter_idr(raw_path, output_path, threshold):
    """IDR column 12 is -log10(global IDR), NOT a p/q-value or local IDR."""
    cutoff = -math.log10(threshold)
    count = 0
    with open(raw_path) as source, open(output_path, "w") as target:
        for line in source:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError("IDR output has fewer than 12 columns")
            score = float(fields[11])
            if not math.isfinite(score):
                raise ValueError("Non-finite global IDR score")
            if score >= cutoff:
                target.write("\t".join(fields[:10]) + "\n")
                count += 1
    return count


def run_pairwise(replicate_paths, pooled_path, outdir, summary, threshold=0.05,
                 rank="p.value", min_peaks=20, executable="idr"):
    if not 0 < threshold < 1 or min_peaks < 20:
        raise ValueError("IDR threshold must be in (0,1), min_peaks >=20")
    if rank not in {"p.value", "signal.value"}:
        raise ValueError("Relaxed p-value peak calls support p.value or signal.value ranking")
    if len(replicate_paths) < 2:
        raise ValueError("IDR requires at least two independent replicates")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    peaks = {name: read_peaks(path) for name, path in replicate_paths.items()}
    pooled = read_peaks(pooled_path)
    for entries in [pooled, *peaks.values()]:
        for peak in entries:
            peak.summit  # validate a genuine narrowPeak/summit, never a BED3 consensus
    indexes = {name: PeakIndex(entries) for name, entries in peaks.items()}
    rows = []
    for i, (a, b) in enumerate(itertools.combinations(sorted(peaks), 2), 1):
        prefix = outdir / f"pair_{i:03d}"
        raw, passed, log = str(prefix) + ".idr.tsv", str(prefix) + ".narrowPeak", str(prefix) + ".log"
        matched = sum(bool(indexes[a].query(p.chrom, p.start, p.end))
                      and bool(indexes[b].query(p.chrom, p.start, p.end)) for p in pooled)
        status, count = "ok", 0
        if min(len(peaks[a]), len(peaks[b]), matched) < min_peaks:
            status = "insufficient_peaks"
            Path(raw).write_text("")
            Path(passed).write_text("")
            Path(log).write_text(f"Not fitted: fewer than {min_peaks} peaks or shared pooled candidates.\n")
        else:
            command = [executable, "--samples", str(replicate_paths[a]), str(replicate_paths[b]),
                       "--peak-list", str(pooled_path), "--input-file-type", "narrowPeak",
                       "--rank", rank, "--soft-idr-threshold", str(threshold),
                       "--use-best-multisummit-IDR", "--random-seed", "0", "--output-file", raw]
            with open(log, "w") as handle:
                subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
            count = filter_idr(raw, passed, threshold)
        rows.append([a, b, len(peaks[a]), len(peaks[b]), matched, threshold, count, status, Path(passed).name])
    Path(summary).parent.mkdir(parents=True, exist_ok=True)
    with open(summary, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["replicate_a", "replicate_b", "peaks_a", "peaks_b", "shared_candidates",
                         "idr_threshold", "passing_peaks", "status", "peak_file"])
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", nargs="+", required=True)
    parser.add_argument("--peaks", nargs="+", required=True)
    parser.add_argument("--pooled", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--threshold", type=float, default=0.05)
    parser.add_argument("--rank", choices=["p.value", "signal.value"], default="p.value")
    parser.add_argument("--min-peaks", type=int, default=20)
    args = parser.parse_args()
    if len(args.replicates) != len(args.peaks) or len(set(args.replicates)) != len(args.replicates):
        parser.error("Replicate identities must map one-to-one to input peaks")
    run_pairwise(dict(zip(args.replicates, args.peaks)), args.pooled, args.outdir, args.summary,
                 args.threshold, args.rank, args.min_peaks)


if __name__ == "__main__":
    main()
