#!/usr/bin/env python3
"""TE quantification with fractional (1/N) weighting.

Unit of counting: TEMPLATE (query name). For every template, N = number of
alignments reported in the input BAM after TE-mode filtering (primary +
secondary, both mates for PE).

Memory model (bounded, no per-qname dictionary):
  Pass 1: `samtools view | cut -f1 | sort | uniq -c` -> qname -> N counts
          (disk-based; memory bounded by sort's --buffer-size).
  Pass 2: `samtools sort -n` (name-sorted) stream merged with the sorted
          counts file -> N per template without loading qnames into memory.

Assignment of an alignment to a TE record:
  * the alignment is converted to (chrom, start, end) using its CIGAR;
  * all TE records overlapping the alignment are found;
  * the record with the LONGEST overlap (bp) wins ("best-overlap" primary
    assignment); nested/overlapping TE records are therefore not counted
    multiple times per alignment (multi-overlap alignments are reported as a
    diagnostic metric).

Fractional weighting:
  * each TE-overlapping alignment contributes 1/N to the family of its
    best-overlapping TE record;
  * a template with 2 alignments hitting 2 different families contributes
    0.5 + 0.5 (i.e. the template weight is spread over its alignments);
  * no alignment is counted twice.

Metrics:
  * total_alignments              - alignments passing the exclude filter
  * secondary_alignments          - alignments flagged 0x100
  * total_templates               - unique query names
  * multimapped_templates         - templates carrying >=1 secondary alignment
  * fraction_templates_multimapped - multimapped_templates / total_templates
  * templates_in_te_unweighted    - templates with >=1 TE-overlapping alignment
  * fraction_templates_in_te_unweighted
  * weighted_in_te                - sum of 1/N over TE-overlapping alignments
  * fraction_weighted_in_te       - weighted_in_te / total_templates
  * multi_overlap_alignments      - alignments overlapping >1 TE record
  * fragments_after_processing    - proper-pair fragments (PE) or reads (SE)

RPM per family = weighted_count / (total_templates / 1e6).
"""

import argparse
import bisect
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

CIGAR_OPS = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_cigar_end(cigar, pos):
    """Return 1-based inclusive end of an alignment from CIGAR string."""
    ref_len = 0
    for n, op in CIGAR_OPS.findall(cigar):
        if op in "MDN=X":
            ref_len += int(n)
    return pos + ref_len - 1


class TEAnnotation:
    """Sorted per-chromosome TE records for fast overlap lookup.

    For each chromosome, records are sorted by start; a ``max_end`` prefix
    array allows correct backward scanning (a record with earlier start can
    still extend past the query start, e.g. nested TE copies).
    """

    def __init__(self, bed_path, family_col, subfamily_col, class_col,
                 min_overlap_frac=0.0, min_overlap_bp=1):
        self.family_col = family_col - 1
        self.subfamily_col = subfamily_col - 1
        self.class_col = class_col - 1
        self.min_overlap_frac = min_overlap_frac
        self.min_overlap_bp = max(1, int(min_overlap_bp))
        self.records = {}  # chrom -> list of (start, end, family, subfamily, te_class)
        self.starts = {}   # chrom -> sorted list of start positions
        self.max_end = {}  # chrom -> prefix max of ends
        self._load(bed_path)

    def _load(self, bed_path):
        ncol = None
        with open(bed_path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue
                if ncol is None:
                    ncol = len(fields)
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                family = fields[self.family_col] if self.family_col < len(fields) else fields[3]
                subfamily = fields[self.subfamily_col] if self.subfamily_col < len(fields) else ""
                te_class = fields[self.class_col] if self.class_col < len(fields) else ""
                family = str(family).split("|")[0].strip()
                subfamily = str(subfamily).split("|")[0].strip()
                self.records.setdefault(chrom, []).append(
                    (start, end, family, subfamily, te_class)
                )
        for chrom in self.records:
            self.records[chrom].sort(key=lambda r: (r[0], r[1]))
            self.starts[chrom] = [r[0] for r in self.records[chrom]]
            prefix = []
            best = 0
            for r in self.records[chrom]:
                best = max(best, r[1])
                prefix.append(best)
            self.max_end[chrom] = prefix

    def best_overlap(self, chrom, start, end):
        """Return best-overlap record for an alignment, or None.

        Returns (family, subfamily, te_class, overlap_bp, n_overlaps).
        """
        recs = self.records.get(chrom)
        if not recs:
            return None
        starts = self.starts[chrom]
        max_end = self.max_end[chrom]
        best = None
        best_bp = -1
        n_overlaps = 0
        align_len = max(1, end - start + 1)

        def consider(i):
            nonlocal best, best_bp, n_overlaps
            rstart, rend = recs[i][0], recs[i][1]
            if rend <= start or rstart >= end:
                return
            ov = min(end, rend) - max(start, rstart)
            if ov < self.min_overlap_bp:
                return
            if self.min_overlap_frac > 0 and ov / align_len < self.min_overlap_frac:
                return
            n_overlaps += 1
            if ov > best_bp:
                best_bp = ov
                best = recs[i]

        # forward: records starting within [start, end)
        i_lo = bisect.bisect_left(starts, start)
        i_hi = bisect.bisect_left(starts, end)
        for i in range(i_lo, i_hi):
            consider(i)
        # backward: records starting before `start` whose span may still
        # extend past it (checked via prefix max_end, correct for nesting)
        i = i_lo - 1
        while i >= 0 and max_end[i] > start:
            consider(i)
            i -= 1
        if best is None:
            return None
        return best[2], best[3], best[4], best_bp, n_overlaps


def iterate_name_sorted(bam, exclude_flags, threads=1):
    """Yield alignment lines from a name-sorted, flag-filtered stream.

    ``samtools sort -n`` sorts by query name (numeric-aware), guaranteeing
    that all alignments of one template are consecutive. The exclude-flags
    filter is applied as a separate ``samtools view -F`` stage feeding the
    name sort, so the stream only contains alignments that pass the filter
    (same semantics as the docs: unmapped/QC-fail/etc. excluded from both N
    and the denominators).
    """
    view_cmd = [
        "samtools", "view", "-b", "-F", str(exclude_flags), bam,
    ]
    view_err = tempfile.TemporaryFile(mode="w+t")
    sort_err = tempfile.TemporaryFile(mode="w+t")
    view = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=view_err)
    sort_cmd = [
        "samtools", "sort", "-n",
        # Reserve one CPU each for Python, samtools view, and the sort main
        # thread. ``-@`` counts only additional sort threads.
        "-@", str(max(0, threads - 3)),
        "-O", "SAM", "-o", "-", "-",
    ]
    proc = subprocess.Popen(
        sort_cmd, stdin=view.stdout, stdout=subprocess.PIPE,
        stderr=sort_err, text=True, bufsize=1,
    )
    view.stdout.close()
    for line in proc.stdout:
        if line.startswith("@"):
            continue
        yield line.rstrip("\n")
    sort_rc = proc.wait()
    view_rc = view.wait()
    view_err.seek(0)
    sort_err.seek(0)
    view_message = view_err.read().strip()
    sort_message = sort_err.read().strip()
    view_err.close()
    sort_err.close()
    if view_rc != 0 or sort_rc != 0:
        raise RuntimeError(
            "name-sorted BAM stream failed\n"
            f"view rc={view_rc}: {view_message}\n"
            f"sort rc={sort_rc}: {sort_message}"
        )


def quantify_streaming(bam, te_annotation, exclude_flags, read_type, threads=1):
    """Single pass over the name-sorted alignment stream.

    Templates are grouped (consecutive records share a query name), so N =
    number of alignments of the template is known without any per-qname
    dictionary or cross-stream merge. Memory is bounded: no per-qname sets
    are kept; per-template flags are group-local booleans, and totals are
    scalars.
    """
    total_alignments = 0
    secondary_alignments = 0
    total_templates = 0
    templates_in_te = 0
    multimapped_templates = 0
    weighted_in_te = 0.0
    multi_overlap_alignments = 0
    proper_pair_alignments = 0
    family_weighted = defaultdict(float)
    family_alignment = defaultdict(int)
    family_meta = {}  # family -> (subfamily, te_class)

    group = []  # current template's alignment lines
    cur_qname = None

    def flush():
        nonlocal total_templates, templates_in_te, multimapped_templates
        nonlocal total_alignments, secondary_alignments, proper_pair_alignments
        nonlocal weighted_in_te, multi_overlap_alignments
        n = len(group)
        if n == 0:
            return
        total_templates += 1
        has_secondary = False
        te_hit = False
        for line in group:
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            flag = int(fields[1])
            chrom, pos = fields[2], fields[3]
            total_alignments += 1
            if flag & 0x100:
                secondary_alignments += 1
                has_secondary = True
            if read_type == "PE" and (flag & 0x2):
                proper_pair_alignments += 1
            if chrom == "*" or chrom.startswith("_") or pos == "0":
                continue
            # SAM POS is 1-based inclusive; BED is 0-based half-open.
            pos1 = int(pos)
            start = pos1 - 1
            end = parse_cigar_end(fields[5], pos1)
            hit = te_annotation.best_overlap(chrom, start, end)
            if hit is None:
                continue
            family, subfamily, te_class, overlap_bp, n_ov = hit
            te_hit = True
            if n_ov > 1:
                multi_overlap_alignments += 1
            weight = 1.0 / n
            weighted_in_te += weight
            family_weighted[family] += weight
            family_alignment[family] += 1
            family_meta.setdefault(family, (subfamily, te_class))
        if has_secondary:
            multimapped_templates += 1
        if te_hit:
            templates_in_te += 1
        group.clear()

    for line in iterate_name_sorted(bam, exclude_flags, threads):
        qname = line.split("\t", 1)[0]
        if cur_qname is None:
            cur_qname = qname
        if qname != cur_qname:
            flush()
            cur_qname = qname
        group.append(line)
    if group:
        flush()

    metrics = {
        "total_alignments": total_alignments,
        "secondary_alignments": secondary_alignments,
        "total_templates": total_templates,
        "multimapped_templates": multimapped_templates,
        "fraction_templates_multimapped": (
            f"{multimapped_templates / total_templates:.6f}" if total_templates else "0"
        ),
        "templates_in_te_unweighted": templates_in_te,
        "fraction_templates_in_te_unweighted": (
            f"{templates_in_te / total_templates:.6f}" if total_templates else "0"
        ),
        "weighted_in_te": f"{weighted_in_te:.6f}",
        "fraction_weighted_in_te": (
            f"{weighted_in_te / total_templates:.6f}" if total_templates else "0"
        ),
        "multi_overlap_alignments": multi_overlap_alignments,
        "fragments_after_processing": (
            proper_pair_alignments // 2 if read_type == "PE" else total_alignments
        ),
    }
    return metrics, family_weighted, family_alignment, family_meta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--te-bed", required=True)
    parser.add_argument("--exclude-flags", type=int, default=516)
    parser.add_argument("--read-type", choices=["PE", "SE"], default="PE")
    parser.add_argument("--family-column", type=int, default=9)
    parser.add_argument("--subfamily-column", type=int, default=8)
    parser.add_argument("--class-column", type=int, default=7)
    parser.add_argument("--min-overlap-frac", type=float, default=0.0)
    parser.add_argument("--min-overlap-bp", type=int, default=1)
    parser.add_argument("--sample", default="sample")
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-family", required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    te = TEAnnotation(
        args.te_bed, args.family_column, args.subfamily_column, args.class_column,
        args.min_overlap_frac, args.min_overlap_bp,
    )

    metrics, family_weighted, family_alignment, family_meta = quantify_streaming(
        args.bam, te, args.exclude_flags, args.read_type, args.threads
    )

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"sample\t{args.sample}\n")
        for key, value in metrics.items():
            out.write(f"{key}\t{value}\n")

    rpm_denom = max(1.0, metrics["total_templates"] / 1e6)
    with open(args.out_family, "w") as out:
        out.write("family\tsubfamily\tclass\tweighted_count\talignment_count\trpm\n")
        for family in sorted(family_weighted, key=lambda f: -family_weighted[f]):
            subfamily, te_class = family_meta.get(family, ("", ""))
            out.write(
                f"{family}\t{subfamily}\t{te_class}\t{family_weighted[family]:.4f}\t"
                f"{family_alignment[family]}\t{family_weighted[family] / rpm_denom:.4f}\n"
            )
    print(f"[te_quant] {args.sample}: {metrics['total_templates']} templates, "
          f"weighted_in_te={metrics['weighted_in_te']}")


if __name__ == "__main__":
    main()
