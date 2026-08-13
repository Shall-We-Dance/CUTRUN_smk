#!/usr/bin/env python3
"""Real paired-end unit test for TE quantification (fractional 1/N).

Constructs a tiny BAM with true PE pairs (flags 99/147 proper pairs, 355/403
secondary mates) and verifies te_quant.py output against hand-computed
expectations, including fragment-vs-mate semantics:

TE annotation:
  TE_copy1: 2L:6980-7520  family Gypsy-1_DM
  TE_copy2: 2L:7700-8000  family Gypsy-2_DM

Templates:
  tA  proper pair, both mates inside TE_copy1 (7000/7090), N=2
      -> both mates overlap TE_copy1: 0.5 + 0.5 = 1.0 Gypsy-1_DM
  tB  proper pair, one mate in TE_copy1 (7000), other mate outside (5000)
      -> 0.5 Gypsy-1_DM (fragment straddles TE boundary: partial weight)
  tC  2-way multimapper pair: mates at 7100 (primary) and 7200 (secondary),
      both in TE_copy1, N=4 (2 mates x 2 positions)
      -> 4 x 0.25 = 1.0 Gypsy-1_DM
  tD  proper pair, both mates in TE_copy2 (7750/7840), N=2 -> 1.0 Gypsy-2_DM
  tE  proper pair, both mates outside TE (5000/5090) -> 0

Expected:
  total_alignments          12   (tA2+tB2+tC4+tD2+tE2)
  secondary_alignments       2   (tC secondary mates)
  total_templates            5
  multimapped_templates      1   (tC)
  templates_in_te_unweighted 4   (tA,tB,tC,tD)
  weighted_in_te             3.5  (1.0+0.5+1.0+1.0)
  multi_overlap_alignments   0
  Gypsy-1_DM 2.5 (tA 1.0 + tB 0.5 + tC 1.0), Gypsy-2_DM 1.0 (tD)
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

OUT = tempfile.mkdtemp(prefix="te_quant_pe_test_")
BAM = os.path.join(OUT, "test.bam")
TE_BED = os.path.join(OUT, "te.bed")

with open(TE_BED, "w") as fh:
    fh.write("2L\t6980\t7520\tTE_copy1\t0\t+\tClass I (retrotransposon)\tLTR/Gypsy\tGypsy-1_DM\n")
    fh.write("2L\t7700\t8000\tTE_copy2\t0\t-\tClass I (retrotransposon)\tLTR/Gypsy\tGypsy-2_DM\n")


def pe_rec(qname, flag, pos, mate_pos, insert, tlen, seq="A"):
    """Proper-pair record. flag: 99 (R1 fwd), 147 (R2 rev), 355 (R1 sec), 403 (R2 sec)."""
    return "\t".join(
        [qname, str(flag), "2L", str(pos), "42", "90M", "=", str(mate_pos), str(tlen),
         seq * 90, "I" * 90, "RG:Z:rg1"]
    )


SAM_LINES = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:2L\tLN:10000",
    "@RG\tID:rg1\tSM:test",
    # tA: both mates in TE_copy1
    pe_rec("tA", 99, 7000, 7090, 180, 180),
    pe_rec("tA", 147, 7090, 7000, 180, -180),
    # tB: one mate in TE_copy1, other outside
    pe_rec("tB", 99, 7000, 5090, 180, 180),
    pe_rec("tB", 147, 5090, 7000, 180, -180),
    # tC: 2-way multimapper proper pair (primary + secondary mates)
    pe_rec("tC", 99, 7100, 7190, 180, 180, "C"),
    pe_rec("tC", 147, 7190, 7100, 180, -180, "C"),
    pe_rec("tC", 355, 7200, 7290, 180, 180, "C"),
    pe_rec("tC", 403, 7290, 7200, 180, -180, "C"),
    # tD: both mates in TE_copy2
    pe_rec("tD", 99, 7750, 7840, 180, 180, "G"),
    pe_rec("tD", 147, 7840, 7750, 180, -180, "G"),
    # tE: outside TE
    pe_rec("tE", 99, 5000, 5090, 180, 180, "T"),
    pe_rec("tE", 147, 5090, 5000, 180, -180, "T"),
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
        "--read-type", "PE",
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
    "total_alignments": "12",
    "secondary_alignments": "2",
    "total_templates": "5",
    "multimapped_templates": "1",
    "templates_in_te_unweighted": "4",
    "weighted_in_te": "3.500000",
    "multi_overlap_alignments": "0",
}
expected_families = {"Gypsy-1_DM": 2.5, "Gypsy-2_DM": 1.0}

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
print("PASS: PE TE quantification matches hand-computed expectations")
print("  weighted_in_te:", metrics["weighted_in_te"])
print("  families:", {k: v["weighted_count"] for k, v in families.items()})
