#!/usr/bin/env python3
"""Unit test for fragment-level FRiP.

Tiny reference: 2L of length 10000. Peaks: 2L:1000-2000, 2L:5000-6000.

Templates (PE):
  tA  proper pair both mates inside peak1 (1100/1190, 1500/1590)
      -> 1 usable template, 1 fragment in peaks
  tB  proper pair, one mate in peak1 (1100), one mate outside (3000)
      -> fragment span 1100-3180 overlaps peak1 -> counted 1
  tC  proper pair, both mates outside (3000/3090, 3200/3290)
      -> 1 usable, 0 in peaks
  tD  proper pair with ONE secondary alignment (both primary mates at
      1200/1290 -> INSIDE peak1; secondary mates at 1250/1340)
      -> secondary excluded from fragment span; primary span 1200-1380
         overlaps peak1 -> counted once
  tE  orphan (mate mapped at 1500, mate unmapped flag 8)
      -> NOT a proper pair -> 0 usable templates
  tF  proper pair, mates in two DIFFERENT peaks (1100 in peak1, 5290 in
      peak2) -> span 1100-5380 overlaps -> counted 1
  tG  1 bp first mate at SAM POS 2000, the last 1-based base represented by
      BED peak1 [1000,2000), second mate outside -> catches BED/SAM boundary

Expected:
  total_reads          = 16
  mapped_reads         = 13 (primary mapped ends only)
  usable_templates     = 6 (tA,tB,tC,tD,tF,tG)  [tE excluded]
  templates_in_peaks   = 5 (tA,tB,tD,tF,tG)
  fragment_frip        = 5/6
  read_frip            = reads in peaks (primary mapped ends overlapping) /
                         mapped_reads
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCRIPT = os.path.join(ROOT, "workflow/rules/scripts/frip.py")
SAMTOOLS = os.environ.get("SAMTOOLS") or shutil.which("samtools")
if not SAMTOOLS:
    raise SystemExit("samtools not found; activate an environment containing samtools or set SAMTOOLS")

OUT = tempfile.mkdtemp(prefix="frip_test_")
BAM = os.path.join(OUT, "test.bam")
PEAKS = os.path.join(OUT, "peaks.bed")

with open(PEAKS, "w") as fh:
    fh.write("2L\t1000\t2000\tpeak1\n")
    fh.write("2L\t5000\t6000\tpeak2\n")


def pe(qname, flag, pos, mate_pos, tlen, seq="A"):
    return "\t".join(
        [qname, str(flag), "2L", str(pos), "42", "90M", "=", str(mate_pos), str(tlen),
         seq * 90, "I" * 90, "RG:Z:rg1"]
    )


def pe_unmapped_mate(qname, flag, pos):
    return "\t".join(
        [qname, str(flag), "2L", str(pos), "42", "90M", "*", "0", "0",
         "A" * 90, "I" * 90, "RG:Z:rg1"]
    )


def pe_one_base(qname, flag, pos, mate_pos, tlen):
    return "\t".join(
        [qname, str(flag), "2L", str(pos), "42", "1M", "=", str(mate_pos), str(tlen),
         "A", "I", "RG:Z:rg1"]
    )


SAM = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:2L\tLN:10000",
    "@RG\tID:rg1\tSM:test",
    # tA both mates in peak1
    pe("tA", 99, 1100, 1190, 180),
    pe("tA", 147, 1190, 1100, 180),
    # tB one mate in peak1, one outside
    pe("tB", 99, 1100, 3090, 180),
    pe("tB", 147, 3090, 1100, 180),
    # tC both outside
    pe("tC", 99, 3000, 3090, 180),
    pe("tC", 147, 3090, 3000, 180),
    # tD proper pair + secondary mate
    pe("tD", 99, 1200, 1290, 180),
    pe("tD", 147, 1290, 1200, 180),
    pe("tD", 355, 1250, 1340, 180),
    pe("tD", 403, 1340, 1250, 180),
    # tE orphan: mate mapped, other unmapped
    pe_unmapped_mate("tE", 73, 1500),
    pe_unmapped_mate("tE", 133, 0),
    # tF mates in different peaks
    pe("tF", 99, 1100, 5290, 180),
    pe("tF", 147, 5290, 1100, -180),
    # tG: exact 0-based BED / 1-based SAM boundary
    pe_one_base("tG", 99, 2000, 2200, 201),
    pe_one_base("tG", 147, 2200, 2000, -201),
]
sam_path = os.path.join(OUT, "test.sam")
with open(sam_path, "w") as fh:
    fh.write("\n".join(SAM) + "\n")
subprocess.run([SAMTOOLS, "view", "-b", "-o", BAM, sam_path], check=True)
sb = os.path.join(OUT, "sorted.bam")
subprocess.run([SAMTOOLS, "sort", "-o", sb, BAM], check=True)
os.replace(sb, BAM)
subprocess.run([SAMTOOLS, "index", BAM], check=True)

out_tsv = os.path.join(OUT, "frip.tsv")
env = dict(os.environ)
env["SAMTOOLS"] = SAMTOOLS
subprocess.run(
    [sys.executable, SCRIPT, "--bam", BAM, "--peaks", PEAKS,
     "--read-type", "PE", "--sample", "test", "--out", out_tsv, "--threads", "1"],
    check=True, env=env,
)
metrics = {}
with open(out_tsv) as fh:
    for line in fh:
        k, _, v = line.rstrip("\n").partition("\t")
        metrics[k] = v

expected = {
    "total_reads": "16",
    "mapped_reads": "13",
    "usable_templates": "6",
    "templates_in_peaks": "5",
    "fragment_frip": "0.833333",
}
errors = []
for k, exp in expected.items():
    got = metrics.get(k, "")
    if got != exp:
        errors.append(f"{k}: expected {exp}, got {got}")
if errors:
    print("FAIL:")
    for e in errors:
        print(" ", e)
    print("metrics:", metrics)
    sys.exit(1)
print("PASS: fragment-level FRiP matches hand-computed expectations")
print("  fragment_frip:", metrics["fragment_frip"], "read_frip:", metrics["read_frip"])
