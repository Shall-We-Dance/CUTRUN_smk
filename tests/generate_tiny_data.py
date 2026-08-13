#!/usr/bin/env python3
"""Generate tiny test FASTQs / BAMs with known properties.

Reads are drawn from a 10 kb tiny reference (tests/data/tiny_ref.fa), which
is a slice of dm6 chr2L. We plant:
  - unique reads (map uniquely to the tiny ref)
  - a duplicate pair (identical read pair appears twice)
  - a multimapper (read drawn from a duplicated 200 bp block)
  - a read overlapping the blacklist interval
  - a read overlapping the TE interval

Outputs under tests/data/: tiny_R1.fastq.gz, tiny_R2.fastq.gz,
tiny_blacklist.bed, tiny_TE.bed, tiny.chrom.sizes.
"""
import gzip
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
random.seed(42)

REF = "".join(
    line.strip()
    for line in open(os.path.join(DATA, "tiny_ref.fa"))
    if not line.startswith(">")
)
assert len(REF) == 10000, len(REF)


def revcomp(s):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(comp)[::-1]


def extract(seq, start, length):
    """Return (sequence, end) with proper strand."""
    return seq[start : start + length], start + length


READS_R1, READS_R2 = [], []
records = []  # (qname, r1_seq, r2_seq, ref_start, description)


def add_pair(qname, start, read_len=90, insert=160, n_dup=1, desc=""):
    """Add a proper read pair starting at ref position `start`."""
    r1 = REF[start : start + read_len]
    r2_fwd_start = start + insert - read_len
    r2_fwd = REF[r2_fwd_start : r2_fwd_start + read_len]
    if len(r1) < read_len or len(r2_fwd) < read_len:
        return
    r2 = revcomp(r2_fwd)
    for i in range(n_dup):
        name = f"{qname}"
        records.append((name, r1, r2, start, desc + (f".dup{i}" if n_dup > 1 else "")))


# 1) unique reads across the ref
for i in range(200):
    start = 1000 + i * 20
    add_pair(f"unique_{i}", start)

# 2) a duplicated pair (exact duplicate fragments) - appears 4x
add_pair("dup_tmpl", 3500, n_dup=4, desc="duplicate")

# 3) a multimapper: a 200 bp sequence repeated twice in the ref -> reads
#    drawn from copy A map also to copy B
copy_a = 5000
copy_b = 5200
for i in range(20):
    start = copy_a + i * 10
    add_pair(f"multi_{i}", start, desc="multimapper")

# 4) a read overlapping the blacklist interval (planted at 6000-6100)
add_pair("bl_1", 5990, desc="blacklist")

# 5) a read overlapping the TE interval (planted at 7000-7500)
add_pair("te_1", 7000, desc="te")
add_pair("te_2", 7200, desc="te")

with open(os.path.join(DATA, "tiny_R1.fastq"), "w") as f1, open(
    os.path.join(DATA, "tiny_R2.fastq"), "w"
) as f2:
    for idx, (qname, r1, r2, start, desc) in enumerate(records):
        qual = "I" * len(r1)
        f1.write(f"@{qname}\n{r1}\n+\n{qual}\n")
        f2.write(f"@{qname}\n{r2}\n+\n{qual}\n")

with open(os.path.join(DATA, "tiny_R1.fastq.gz"), "wb") as f:
    with gzip.open(f, "wt") as gz:
        for qname, r1, r2, start, desc in records:
            gz.write(f"@{qname}\n{r1}\n+\n{'I' * len(r1)}\n")
with open(os.path.join(DATA, "tiny_R2.fastq.gz"), "wb") as f:
    with gzip.open(f, "wt") as gz:
        for qname, r1, r2, start, desc in records:
            gz.write(f"@{qname}\n{r2}\n+\n{'I' * len(r1)}\n")

# blacklist + TE annotations (1-based BED)
with open(os.path.join(DATA, "tiny_blacklist.bed"), "w") as out:
    out.write("2L\t5995\t6100\tBlacklistTest\n")
with open(os.path.join(DATA, "tiny_TE.bed"), "w") as out:
    out.write(
        "2L\t6980\t7520\tTE_copy1\t0\t+\tClass I (retrotransposon)\tLTR/Gypsy\tGypsy-1_DM\n"
    )
    out.write(
        "2L\t7300\t7600\tTE_copy2\t0\t-\tClass I (retrotransposon)\tLTR/Gypsy\tGypsy-2_DM\n"
    )
with open(os.path.join(DATA, "tiny.chrom.sizes"), "w") as out:
    out.write(f"2L\t{len(REF)}\n")

print(f"generated {len(records)} read pairs")
for d in set(r[4].split(".")[0] for r in records):
    n = sum(1 for r in records if r[4].split(".")[0] == d)
    print(f"  {d}: {n} pairs")
