"""Count validated spike-in templates, using the same PE eligibility as targets."""
import argparse
import os
from pathlib import Path
import subprocess

from fragment_records import template_groups, usable_fragment


def count_templates(bam, read_type, exclude_flags, threads=1):
    header = subprocess.check_output([os.environ.get("SAMTOOLS", "samtools"), "view", "-H", str(bam)], text=True)
    sizes = {}
    for line in header.splitlines():
        if line.startswith("@SQ\t"):
            tags = dict(field.split(":", 1) for field in line.split("\t")[1:])
            sizes[tags["SN"]] = int(tags["LN"])
    count = reads = 0
    for records in template_groups(bam, threads):
        fragment, _ = usable_fragment(records, read_type, exclude_flags, sizes)
        if fragment is not None:
            count += 1
            reads += fragment[3]
    return reads, count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--read-type", choices=["PE", "SE"], required=True)
    parser.add_argument("--exclude-flags", type=int, default=1804)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--scale-to", type=float, default=10000)
    parser.add_argument("--min-fragments", type=int, default=1000)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    reads, count = count_templates(args.bam, args.read_type, args.exclude_flags, args.threads)
    status = "NO_SPIKEIN" if count == 0 else "LOW_SPIKEIN" if count < args.min_fragments else "OK"
    scale = args.scale_to / count if count else 1.0
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as out:
        out.write("sample\tspikein_reads\tspikein_fragments\tscale_to\tscale_factor\tstatus\n")
        out.write(f"{args.sample}\t{reads}\t{count}\t{args.scale_to}\t{scale:.10f}\t{status}\n")
    print(f"Validated templates={count}, status={status}; zero-spike-in retains scale 1, not calibrated abundance.")


if __name__ == "__main__":
    main()
