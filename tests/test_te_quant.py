#!/usr/bin/env python3
"""Unit test for TE quantification (fractional 1/N weighting).

Builds a tiny BAM by hand with known templates (single-end style records;
the weighting logic is strand/end independent) and checks te_quant.py output
against hand-computed expected values.

TE annotation (tests/data independent, created here):
  TE_copy1: 2L:6980-7520  family Gypsy-1_DM (class Class I)
  TE_copy2: 2L:7700-8000  family Gypsy-2_DM
  TE_copy3: 2L:7200-7400  family Jockey        (nested inside TE_copy1)

Records:
  tA  flag 0  7000-7090  -> TE_copy1 only              N=1 -> 1.0 Gypsy-1_DM
  tB  flag 0  7100-7190  -> TE_copy1 only              N=2 -> 0.5 Gypsy-1_DM
      flag 256 7450-7540 -> TE_copy1 only (70 bp)            0.5 Gypsy-1_DM
  tC  flag 0  7000-7090  -> TE_copy1                    N=2 -> 0.5 Gypsy-1_DM
      flag 256 7750-7840 -> TE_copy2                          0.5 Gypsy-2_DM
  tD  flag 0  5000-5090  -> no TE                        N=1 -> 0
  tE  flag 0  7300-7390  -> overlaps TE_copy1 AND TE_copy3 (tie at 90 bp)
                            best-overlap tie -> the record considered first
                            (latest start in backward scan): TE_copy3
                            N=1 -> 1.0 Jockey, multi_overlap=1
  tF  flag 0  7250-7340  -> overlaps TE_copy1 AND TE_copy3 (tie at 90 bp)
                            -> TE_copy3                    N=1 -> 1.0 Jockey,
                            multi_overlap=1
  tG  flag 0  7520-7520  -> last 1-based base of BED [6980,7520)
                            N=1 -> 1.0 Gypsy-1_DM

Expected:
  total_alignments          9
  secondary_alignments      2            (tB 256, tC 256)
  total_templates           7
  multimapped_templates     2            (tB, tC)
  templates_in_te_unweighted 6           (all but tD)
  weighted_in_te            6.0
  multi_overlap_alignments  2            (tE, tF)
  Gypsy-1_DM  weighted 3.5   (tA 1.0 + tB 1.0 + tC 0.5 + tG 1.0)
  Gypsy-2_DM  weighted 0.5   (tC 0.5)
  Jockey      weighted 2.0   (tE 1.0 + tF 1.0)
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCRIPT = os.path.join(ROOT, "workflow/rules/scripts/te_quant.py")
SAMTOOLS = os.environ.get("SAMTOOLS") or shutil.which("samtools")
if not SAMTOOLS:
    raise SystemExit("samtools not found; activate an environment containing samtools or set SAMTOOLS")

OUT = tempfile.mkdtemp(prefix="te_quant_test_")
BAM = os.path.join(OUT, "test.bam")
TE_BED = os.path.join(OUT, "te.bed")

with open(TE_BED, "w") as fh:
    fh.write(
        "2L\t6980\t7520\tTE_copy1\t0\t+\tClass I (retrotransposon)\tLTR/Gypsy\tGypsy-1_DM\n"
    )
    fh.write(
        "2L\t7700\t8000\tTE_copy2\t0\t-\tClass I (retrotransposon)\tLTR/Gypsy\tGypsy-2_DM\n"
    )
    fh.write(
        "2L\t7200\t7400\tTE_copy3\t0\t+\tClass II (DNA transposon)\tTIR\tJockey\n"
    )


def rec(qname, flag, pos, seq="A"):
    return "\t".join(
        [qname, str(flag), "2L", str(pos), "42", "90M", "*", "0", "0",
         seq * 90, "I" * 90, "RG:Z:rg1"]
    )


SAM_LINES = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:2L\tLN:10000",
    "@RG\tID:rg1\tSM:test",
    rec("tA", 0, 7000),
    rec("tB", 0, 7100, "C"),
    rec("tB", 256, 7450, "C"),
    rec("tC", 0, 7000, "G"),
    rec("tC", 256, 7750, "G"),
    rec("tD", 0, 5000, "T"),
    rec("tE", 0, 7300),
    rec("tF", 0, 7250, "T"),
    "\t".join(["tG", "0", "2L", "7520", "42", "1M", "*", "0", "0", "A", "I", "RG:Z:rg1"]),
]

sam_path = os.path.join(OUT, "test.sam")
with open(sam_path, "w") as fh:
    fh.write("\n".join(SAM_LINES) + "\n")

subprocess.run([SAMTOOLS, "view", "-b", "-o", BAM, sam_path], check=True)
sorted_bam = os.path.join(OUT, "test.sorted.bam")
subprocess.run([SAMTOOLS, "sort", "-o", sorted_bam, BAM], check=True)
os.replace(sorted_bam, BAM)
subprocess.run([SAMTOOLS, "index", BAM], check=True)

out_tsv = os.path.join(OUT, "overlap.tsv")
out_family = os.path.join(OUT, "family.tsv")
env = dict(os.environ)
env["PATH"] = os.path.dirname(SAMTOOLS) + os.pathsep + env.get("PATH", "")
subprocess.run(
    [
        sys.executable, SCRIPT,
        "--bam", BAM,
        "--te-bed", TE_BED,
        "--exclude-flags", "516",
        "--read-type", "SE",
        "--sample", "test",
        "--out-tsv", out_tsv,
        "--out-family", out_family,
    ],
    check=True,
    env=env,
)

metrics = {}
with open(out_tsv) as fh:
    for line in fh:
        k, _, v = line.rstrip("\n").partition("\t")
        metrics[k] = v

families = {}
with open(out_family) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        vals = line.rstrip("\n").split("\t")
        row = dict(zip(header, vals))
        families[row["family"]] = row

expected_metrics = {
    "total_alignments": "9",
    "secondary_alignments": "2",
    "total_templates": "7",
    "multimapped_templates": "2",
    "templates_in_te_unweighted": "6",
    "weighted_in_te": "6.000000",
    "multi_overlap_alignments": "2",
}
expected_families = {"Gypsy-1_DM": 3.5, "Gypsy-2_DM": 0.5, "Jockey": 2.0}

errors = []
for key, exp in expected_metrics.items():
    got = metrics.get(key, "")
    if str(got) != exp:
        errors.append(f"metric {key}: expected {exp}, got {got}")

for fam, exp in expected_families.items():
    got = float(families.get(fam, {}).get("weighted_count", 0))
    if abs(got - exp) > 1e-4:
        errors.append(f"family {fam}: expected {exp}, got {got}")

if errors:
    print("FAIL:")
    for e in errors:
        print(" ", e)
    print("metrics:", metrics)
    print("families:", families)
    sys.exit(1)
print("PASS: TE quantification matches hand-computed expectations")
print("  weighted_in_te:", metrics["weighted_in_te"])
print("  families:", {k: v["weighted_count"] for k, v in families.items()})
