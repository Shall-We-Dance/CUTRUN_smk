#!/usr/bin/env python3
"""Compute effective genome size (bp) over contigs kept in bigWig tracks.

Contigs excluded when bigwig.remove_chrM_and_scaffolds is true:
  - names containing '_' (scaffold naming), 'scaffold', 'random', 'un',
    'alt', 'fix', 'hap' (case-insensitive)
  - chrM / MT / M
"""
import argparse
import re
import sys

EXCLUDE_RE = re.compile(r"(?i)(^|_)scaffold|random|_un$|_un_|alt|fix|hap")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--remove-chrm", action="store_true", default=False)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    total = 0
    excluded = 0
    with open(args.chrom_sizes) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            chrom, size = fields[0], int(fields[1])
            if (
                args.remove_chrm
                and (
                    "_" in chrom
                    or EXCLUDE_RE.search(chrom)
                    or chrom.upper() in {"CHRM", "MT", "M"}
                )
            ):
                excluded += size
                continue
            total += size
    with open(args.out, "w") as out:
        out.write(f"{total}\n")
    print(f"[genome_size] total={total} excluded={excluded}", file=sys.stderr)


if __name__ == "__main__":
    main()
