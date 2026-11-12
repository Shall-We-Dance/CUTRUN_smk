"""Shared fragment eligibility and coordinates. No third-party dependencies."""
import itertools
import os
import re
import subprocess
import tempfile
from collections import defaultdict


def excluded_contig(chrom):
    return ("_" in chrom or chrom.upper() in {"CHRM", "MT", "M"}
            or re.search(r"scaffold|random|alt|fix|hap", chrom, re.I) is not None)


def chrom_sizes(path, remove_chrm=False):
    sizes = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, size = line.split()[:2]
            if chrom in sizes or int(size) <= 0:
                raise ValueError(f"Invalid or duplicate chromosome size: {line.strip()}")
            if not remove_chrm or not excluded_contig(chrom):
                sizes[chrom] = int(size)
    if not sizes:
        raise ValueError("No chromosomes remain after filtering")
    return sizes


def reference_end(record):
    start = int(record[3]) - 1
    cigar = record[5]
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if not ops or "".join(n + op for n, op in ops) != cigar:
        raise ValueError(f"Invalid CIGAR: {cigar}")
    return start + sum(int(n) for n, op in ops if op in "MDN=X")


def template_groups(bam, threads=1):
    """Disk-backed name sort; identical query names in different RGs stay separate."""
    cmd = [os.environ.get("SAMTOOLS", "samtools"), "sort", "-n", "-@",
           str(max(0, threads - 2)), "-O", "SAM", "-o", "-", str(bam)]
    with tempfile.TemporaryFile(mode="w+t") as errors:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errors, text=True)
        try:
            rows = (line.rstrip("\n").split("\t") for line in proc.stdout if not line.startswith("@"))
            for _, group in itertools.groupby(rows, key=lambda r: r[0]):
                by_rg = defaultdict(list)
                for record in group:
                    if len(record) < 11:
                        raise ValueError("Truncated SAM record")
                    rg = next((tag[5:] for tag in record[11:] if tag.startswith("RG:Z:")), "")
                    by_rg[rg].append(record)
                yield from by_rg.values()
            if proc.wait() != 0:
                errors.seek(0)
                raise RuntimeError(f"samtools name-sort failed: {errors.read()}")
        finally:
            proc.stdout.close()
            if proc.poll() is None:
                proc.terminate()
                proc.wait()


def usable_fragment(records, read_type, exclude_flags, sizes, se_extension=200):
    """PE: exactly one primary R1+R2, both proper, both passing, same chromosome.

    Span comes from both CIGARs, not stale TLEN. SE uses an explicit 5'-anchored
    extension clipped to chromosome boundaries. Returns (fragment, reason).
    """
    primary = [r for r in records if not int(r[1]) & (0x100 | 0x800)]
    if read_type == "PE":
        r1 = [r for r in primary if int(r[1]) & 0x40 and not int(r[1]) & 0x80]
        r2 = [r for r in primary if int(r[1]) & 0x80 and not int(r[1]) & 0x40]
        if len(primary) != 2 or len(r1) != 1 or len(r2) != 1:
            return None, "incomplete_or_ambiguous_pair"
        mates = [r1[0], r2[0]]
        if any(not int(r[1]) & 2 or int(r[1]) & (exclude_flags | 0xC) for r in mates):
            return None, "flags_or_improper_pair"
        if mates[0][2] != mates[1][2]:
            return None, "different_chromosomes"
        start = min(int(r[3]) - 1 for r in mates)
        end = max(reference_end(r) for r in mates)
    else:
        if len(primary) != 1:
            return None, "ambiguous_primary"
        mates = primary
        record = primary[0]
        if int(record[1]) & (exclude_flags | 4):
            return None, "flags_or_unmapped"
        if int(record[1]) & 16:
            end = reference_end(record)
            start = max(0, end - se_extension)
        else:
            start = int(record[3]) - 1
            end = start + se_extension
    chrom = mates[0][2]
    if chrom not in sizes:
        return None, "excluded_chromosome"
    if read_type == "SE":
        end = min(end, sizes[chrom])
    if not 0 <= start < end <= sizes[chrom]:
        return None, "invalid_span"
    return (chrom, start, end, len(mates)), None


def bed_fragments(path):
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            chrom, start, end = line.split()[:3]
            start, end = int(start), int(end)
            if not 0 <= start < end:
                raise ValueError(f"{path}:{number}: invalid fragment")
            yield chrom, start, end
