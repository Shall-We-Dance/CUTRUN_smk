#!/usr/bin/env python3
"""Order-independent replicate peak consensus and concordance metrics.

For each pair of replicates (i, j):
  - peak overlap in both directions (1 bp overlap, bedtools intersect -u)
  - base-pair Jaccard (bedtools jaccard)

Symmetric consensus (order-independent, works for >=2 replicates):
  a peak of replicate i is part of the consensus iff it overlaps (>=1 bp) at
  least one peak in EVERY other replicate. The consensus peak set is the
  union over all replicates i of such peaks.

Outputs:
  TSV rows (sample, metric, value):
    rep_peak_files, rep{i}_peaks_overlapping_rep{j} (i<j, both directions),
    consensus_peaks (symmetric, order-independent), union_peaks (merged
    interval count), pairwise bp Jaccard (rep{i}_vs_rep{j}_bp_jaccard)
  BED: symmetric consensus intervals
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

BEDTOOLS = os.environ.get("BEDTOOLS", "bedtools")


def run(cmd, log=None):
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if log:
        with open(log, "a") as fh:
            fh.write(proc.stdout)
            fh.write(proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed: {cmd}\n{proc.stderr}")
    return proc.stdout


def peak_count(path):
    with open(path) as fh:
        return sum(1 for line in fh if line.strip())


def main():
    sample = snakemake.wildcards.sample
    peaks = snakemake.input.peaks  # list of per-replicate narrowPeak paths
    out_tsv = snakemake.output.tsv
    out_cons = snakemake.output.consensus
    log = snakemake.log[0]
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    Path(out_cons).parent.mkdir(parents=True, exist_ok=True)
    with open(log, "w"):
        pass

    peaks = sorted(peaks)
    rows = [("sample", "metric", "value"), (sample, "rep_peak_files", str(len(peaks)))]
    for idx, p in enumerate(peaks, 1):
        rows.append((sample, f"rep{idx}_peaks", str(peak_count(p))))

    if len(peaks) < 2:
        # no biological replicate to assess consistency
        rows.append((sample, "symmetric_consensus_peaks", "NA"))
        rows.append((sample, "union_peaks", "NA"))
        rows.append((sample, "consensus_status", "single_replicate_only"))
        with open(out_tsv, "w") as out:
            for row in rows:
                out.write("\t".join(str(v) for v in row) + "\n")
        with open(out_cons, "w"):
            pass  # empty consensus file
        print(f"[replicate_consensus] {sample}: single replicate; consensus=NA")
        return

    tmpdir = tempfile.mkdtemp(prefix="repl_cons_")
    try:
        overlap_any = []  # per-replicate: peaks overlapping all others
        for i, pi in enumerate(peaks):
            with_others = pi
            for j, pj in enumerate(peaks):
                if i == j:
                    continue
                out = os.path.join(tmpdir, f"ov_{i}_{j}.bed")
                run(f"{BEDTOOLS} intersect -a '{with_others}' -b '{pj}' -u > '{out}' 2>> {log}")
                with_others = out
            overlap_any.append(with_others)

        # symmetric consensus = union over replicates of peaks present in all
        tmp_union = os.path.join(tmpdir, "union_cons.bed")
        run(
            f"cat {' '.join(overlap_any)} | {BEDTOOLS} sort -i - | {BEDTOOLS} merge -i - "
            f"> '{tmp_union}' 2>> {log}"
        )
        consensus_count = peak_count(tmp_union)
        # merged consensus intervals (one record per merged region)
        run(f"cp '{tmp_union}' '{out_cons}' 2>> {log}")

        # pairwise: bidirectional overlap counts + bp Jaccard
        n = len(peaks)
        for i in range(n):
            for j in range(i + 1, n):
                a, b = peaks[i], peaks[j]
                out_ij = os.path.join(tmpdir, f"pair_{i}_{j}.bed")
                run(f"{BEDTOOLS} intersect -a '{a}' -b '{b}' -u > '{out_ij}' 2>> {log}")
                rows.append(
                    (sample, f"rep{i+1}_peaks_overlapping_rep{j+1}", str(peak_count(out_ij)))
                )
                out_ji = os.path.join(tmpdir, f"pair_{j}_{i}.bed")
                run(f"{BEDTOOLS} intersect -a '{b}' -b '{a}' -u > '{out_ji}' 2>> {log}")
                rows.append(
                    (sample, f"rep{j+1}_peaks_overlapping_rep{i+1}", str(peak_count(out_ji)))
                )
                jac = run(
                    f"{BEDTOOLS} jaccard -a '{a}' -b '{b}' 2>> {log}"
                ).splitlines()[1].split("\t")
                rows.append((sample, f"rep{i+1}_vs_rep{j+1}_bp_jaccard", jac[2]))

        # union peaks (merged interval count across all replicates)
        tmp_all = os.path.join(tmpdir, "all_merged.bed")
        run(
            f"cat {' '.join(peaks)} | {BEDTOOLS} sort -i - | {BEDTOOLS} merge -i - "
            f"> '{tmp_all}' 2>> {log}"
        )
        rows.append((sample, "union_peaks", str(peak_count(tmp_all))))
        rows.append((sample, "symmetric_consensus_peaks", str(consensus_count)))
    finally:
        import shutil

        shutil.rmtree(tmpdir, ignore_errors=True)

    with open(out_tsv, "w") as out:
        for row in rows:
            out.write("\t".join(str(v) for v in row) + "\n")
    print(f"[replicate_consensus] {sample}: consensus={consensus_count}")


if __name__ == "__main__":
    main()
