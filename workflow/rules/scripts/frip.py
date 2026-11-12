#!/usr/bin/env python3
"""Fragment-level FRiP (fraction of reads in peaks).

Unit of counting: TEMPLATE (query name) — NOT alignment ends.

For PE data, a fragment is the span of a properly-paired template (both
mates mapped, flag 0x2). Each template is counted exactly once:
  - both mates inside peaks  -> fragment counted once
  - exactly one mate inside  -> fragment counted once (span-based presence)
  - secondary/supplementary alignments are excluded from fragment counting
  - a fragment overlapping several peaks is still counted once
  - orphan/singleton templates are excluded from the denominator

For SE data, each mapped read (template) is one fragment.

Implementation (memory-bounded, name-sorted stream):
  1. `samtools sort -n` groups a template's records together.
  2. For each template: decide inclusion (PE: both mates present and flag 2;
     SE: mapped) and fragment span; count 1 if any part of the span overlaps
     a peak (1 bp overlap, tested via a sorted peak interval index).
  3. denominators: total fragments (usable templates) and total reads
     (alignment ends, for read-level reporting).

Outputs (metric/value):
  total_reads, mapped_reads, usable_templates,
  fragments_in_peaks, fragment_frip, reads_in_peaks, read_frip
"""

import argparse
import bisect
import os
import subprocess
import sys
import tempfile
from pathlib import Path

SAMTOOLS = os.environ.get("SAMTOOLS", "samtools")


def parse_cigar_end(cigar, pos):
    ref_len = 0
    i = 0
    n = 0
    while i < len(cigar):
        c = cigar[i]
        if c.isdigit():
            n = n * 10 + int(c)
        else:
            if c in "MDN=X":
                ref_len += n
            n = 0
        i += 1
    return pos + ref_len - 1


class PeakIndex:
    """Sorted per-chromosome peak intervals for 1 bp overlap queries."""

    def __init__(self, peaks_path):
        self.starts = {}
        self.ends = {}
        self.max_end = {}
        with open(peaks_path) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                chrom = f[0]
                s, e = int(f[1]), int(f[2])
                self.starts.setdefault(chrom, []).append(s)
                self.ends.setdefault(chrom, []).append(e)
        for chrom in self.starts:
            pairs = sorted(zip(self.starts[chrom], self.ends[chrom]))
            self.starts[chrom] = [p[0] for p in pairs]
            self.ends[chrom] = [p[1] for p in pairs]
            prefix = []
            best = 0
            for e in self.ends[chrom]:
                best = max(best, e)
                prefix.append(best)
            self.max_end[chrom] = prefix

    def overlaps(self, chrom, start, end):
        if chrom not in self.starts:
            return False
        starts = self.starts[chrom]
        max_end = self.max_end[chrom]
        i_lo = bisect.bisect_left(starts, start)
        i_hi = bisect.bisect_left(starts, end)
        for i in range(i_lo, min(i_hi, len(starts))):
            if start < self.ends[chrom][i] and end > starts[i]:
                return True
        i = i_lo - 1
        while i >= 0 and max_end[i] > start:
            if start < self.ends[chrom][i] and end > starts[i]:
                return True
            i -= 1
        return False


def iterate_name_sorted(bam, threads=1):
    cmd = [
        SAMTOOLS, "sort", "-n",
        # Reserve one CPU for this Python process and one for samtools' main
        # thread. ``-@`` counts only additional samtools threads.
        "-@", str(max(0, threads - 2)),
        "-O", "SAM", "-o", "-", bam,
    ]
    err = tempfile.TemporaryFile(mode="w+t")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, text=True, bufsize=1)
    for line in proc.stdout:
        if line.startswith("@"):
            continue
        yield line.rstrip("\n")
    rc = proc.wait()
    err.seek(0)
    message = err.read().strip()
    err.close()
    if rc != 0:
        raise RuntimeError(f"samtools sort -n failed: {' '.join(cmd)}\n{message}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--bam", help="Legacy standalone BAM mode (unextended SE reads)")
    source.add_argument("--fragments", help="Shared pipeline BED3 fragment population")
    ap.add_argument("--peaks", required=True)
    ap.add_argument("--read-type", choices=["PE", "SE"], default="PE")
    ap.add_argument("--sample", default="sample")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    peaks = PeakIndex(args.peaks)
    if args.fragments:
        from fragment_records import bed_fragments
        total = inside = 0
        for chrom, start, end in bed_fragments(args.fragments):
            total += 1
            inside += peaks.overlaps(chrom, start, end)
        with open(args.out, "w") as out:
            out.write("metric\tvalue\n")
            out.write(f"sample\t{args.sample}\ncounting_unit\tfragment\n")
            out.write(f"usable_templates\t{total}\ntemplates_in_peaks\t{inside}\n")
            out.write(f"fragment_frip\t{inside / total if total else 0:.6f}\n")
            # Read-end FRiP cannot be reconstructed from a fragment span.
            for key in ["total_reads", "mapped_reads", "reads_in_peaks", "read_frip"]:
                out.write(f"{key}\tNA\n")
        print(f"[frip] {args.sample}: shared fragments {inside}/{total}")
        return
    total_reads = 0
    mapped_reads = 0
    usable_templates = 0
    templates_in_peaks = 0
    reads_in_peaks = 0
    template_reads = 0

    group = []
    cur_qname = None

    def flush():
        nonlocal usable_templates, templates_in_peaks, reads_in_peaks
        if not group:
            return
        if args.read_type == "SE":
            # one template per alignment (already mapped by caller filter? no:
            # we count all records of the name group)
            for line in group:
                f = line.split("\t")
                flag = int(f[1])
                if flag & (0x4 | 0x100 | 0x800):
                    continue
                usable_templates += 1
                chrom, pos = f[2], f[3]
                if chrom == "*" or pos == "0":
                    continue
                pos1 = int(pos)
                end = parse_cigar_end(f[5], pos1)
                if peaks.overlaps(chrom, pos1 - 1, end):
                    templates_in_peaks += 1
            return

        # PE: proper pair present?
        recs = [g.split("\t") for g in group if len(g.split("\t")) >= 6]
        primary = [r for r in recs if not (int(r[1]) & (0x4 | 0x100 | 0x800))]
        proper = [r for r in primary if (int(r[1]) & 0x2)]
        read1 = [r for r in proper if int(r[1]) & 0x40]
        read2 = [r for r in proper if int(r[1]) & 0x80]
        if read1 and read2:
            usable_templates += 1
            r1, r2 = read1[0], read2[0]
            c1, p1, c2, p2 = r1[2], int(r1[3]), r2[2], int(r2[3])
            e1 = parse_cigar_end(r1[5], p1)
            e2 = parse_cigar_end(r2[5], p2)
            if c1 == c2:
                lo, hi = min(p1 - 1, p2 - 1), max(e1, e2)
                if peaks.overlaps(c1, lo, hi):
                    templates_in_peaks += 1
            elif peaks.overlaps(c1, p1 - 1, e1) or peaks.overlaps(c2, p2 - 1, e2):
                templates_in_peaks += 1
        # read-level: mapped primary alignment ends
        for r in primary:
            flag = int(r[1])
            if r[2] != "*" and r[3] != "0":
                pos1 = int(r[3])
                end = parse_cigar_end(r[5], pos1)
                if peaks.overlaps(r[2], pos1 - 1, end):
                    reads_in_peaks += 1

    for line in iterate_name_sorted(args.bam, args.threads):
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        qname = fields[0]
        flag = int(fields[1])
        total_reads += 1
        if not (flag & (0x4 | 0x100 | 0x800)):
            mapped_reads += 1
        if cur_qname is None:
            cur_qname = qname
        if qname != cur_qname:
            flush()
            group = []
            cur_qname = qname
        group.append(line)
    flush()

    fragment_frip = templates_in_peaks / usable_templates if usable_templates else 0.0
    read_frip = reads_in_peaks / mapped_reads if mapped_reads else 0.0
    with open(args.out, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"sample\t{args.sample}\n")
        out.write(f"total_reads\t{total_reads}\n")
        out.write(f"mapped_reads\t{mapped_reads}\n")
        out.write(f"usable_templates\t{usable_templates}\n")
        out.write(f"templates_in_peaks\t{templates_in_peaks}\n")
        out.write(f"fragment_frip\t{fragment_frip:.6f}\n")
        out.write(f"reads_in_peaks\t{reads_in_peaks}\n")
        out.write(f"read_frip\t{read_frip:.6f}\n")
    print(f"[frip] {args.sample}: fragments {templates_in_peaks}/{usable_templates} "
          f"fragment_frip={fragment_frip:.4f} read_frip={read_frip:.4f}")


if __name__ == "__main__":
    main()
