"""Shared BED3 fragment population for MACS3, FRiP, tracks and QC."""
import argparse
from collections import Counter
import os
from pathlib import Path
import subprocess
import tempfile

from consensus import PeakIndex, read_peaks
from fragment_records import chrom_sizes, template_groups, usable_fragment


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--blacklist")
    parser.add_argument("--remove-chrm", action="store_true")
    parser.add_argument("--read-type", choices=["PE", "SE"], required=True)
    parser.add_argument("--exclude-flags", type=int, default=1804)
    parser.add_argument("--se-extension", type=int, default=200)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--lengths", required=True)
    args = parser.parse_args()
    if args.se_extension < 1:
        parser.error("--se-extension must be positive")
    sizes = chrom_sizes(args.chrom_sizes, args.remove_chrm)
    blacklist = PeakIndex(read_peaks(args.blacklist)) if args.blacklist else None
    counts, lengths = Counter(), Counter()
    for path in [args.out, args.metrics, args.lengths]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="fragments_") as tempdir:
        unsorted = Path(tempdir) / "fragments.bed"
        with unsorted.open("w") as out:
            for records in template_groups(args.bam, args.threads):
                counts["total_reads_after_processing"] += len(records)
                counts["mapped_reads_after_processing"] += sum(not int(r[1]) & (4 | 256 | 2048) for r in records)
                fragment, reason = usable_fragment(records, args.read_type, args.exclude_flags, sizes, args.se_extension)
                if fragment is None:
                    counts[f"rejected_{reason}"] += 1
                    continue
                chrom, start, end, mates = fragment
                if blacklist and blacklist.query(chrom, start, end):
                    counts["rejected_blacklist_span"] += 1
                    continue
                out.write(f"{chrom}\t{start}\t{end}\n")
                counts["usable_fragments_after_processing"] += 1
                counts["usable_reads_after_processing"] += mates
                counts["total_fragment_bases"] += end - start
                lengths[end - start] += 1
        # Identical coordinates from independent molecules are intentionally retained.
        with open(args.out, "w") as out:
            subprocess.run(["sort", "-k1,1", "-k2,2n", "-k3,3n", "-T", tempdir, str(unsorted)],
                           stdout=out, check=True, env={**os.environ, "LC_ALL": "C"})
    with open(args.metrics, "w") as out:
        out.write(f"metric\tvalue\nsample\t{args.sample}\nread_type\t{args.read_type}\n")
        for key in ["total_reads_after_processing", "mapped_reads_after_processing", "usable_reads_after_processing",
                    "usable_fragments_after_processing", "total_fragment_bases"]:
            out.write(f"{key}\t{counts[key]}\n")
        proper = counts["usable_reads_after_processing"] if args.read_type == "PE" else "NA"
        out.write(f"proper_pair_reads_after_processing\t{proper}\n")
        out.write(f"max_fragment_length\t{max(lengths, default=0)}\n")
        for key in sorted(counts):
            if key.startswith("rejected_"):
                out.write(f"{key}\t{counts[key]}\n")
    with open(args.lengths, "w") as out:
        out.write("fragment_length\tcount\n")
        for length, count in sorted(lengths.items()):
            out.write(f"{length}\t{count}\n")


if __name__ == "__main__":
    main()
