"""Post-run assertions for all seven tiny fixtures; requires pyBigWig.

Run after the integration commands in README.md. Missing results fail rather
than silently skipping verification. No experimental output directory is read.
"""
from collections import Counter
import csv
import json
from pathlib import Path
import unittest

import pyBigWig

ROOT = Path(__file__).resolve().parent
FIXTURES = ["tiny", "tiny_se", "tiny_te", "tiny_replicates",
            "tiny_broad", "tiny_consensus", "tiny_spikein"]


def table(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fragments(path):
    with open(path) as handle:
        return [(r[0], int(r[1]), int(r[2])) for r in csv.reader(handle, delimiter="\t") if r]


class IntegrationOutputs(unittest.TestCase):
    def test_frip_qc_and_fragment_counts_agree(self):
        for fixture in FIXTURES:
            directory = ROOT / f"out_{fixture}"
            qc = {r["sample"]: r for r in table(directory / "qc/cutrun/cutrun_qc_summary.tsv")}
            for row in table(directory / "qc/frip/frip_summary.tsv"):
                sample = row["sample"]
                if row["level"] == "pooled":
                    bed = directory / f"fragments/pooled/{sample}.pooled.fragments.bed"
                else:
                    bed = directory / f"fragments/{sample}.fragments.bed"
                count = len(fragments(bed))
                with self.subTest(fixture=fixture, sample=sample, level=row["level"]):
                    self.assertGreater(count, 0)
                    self.assertEqual(int(row["usable_templates"]), count)
                    self.assertEqual(row["read_frip"], "NA")
                    self.assertLessEqual(int(row["templates_in_peaks"]), count)
                    self.assertAlmostEqual(float(row["fragment_frip"]),
                                           int(row["templates_in_peaks"]) / count, places=6)
                    if row["level"] != "pooled":
                        self.assertEqual(int(qc[sample]["usable_fragments_after_processing"]), count)
                        mates = 1 if qc[sample]["read_type"] == "SE" else 2
                        self.assertEqual(int(qc[sample]["usable_reads_after_processing"]), mates * count)

    def test_pool_preserves_all_independent_molecules(self):
        for fixture in ["tiny_replicates", "tiny_broad", "tiny_consensus"]:
            directory = ROOT / f"out_{fixture}/fragments"
            pools = list((directory / "pooled").glob("*.pooled.fragments.bed"))
            self.assertTrue(pools)
            for pool in pools:
                sample = pool.name.removesuffix(".pooled.fragments.bed")
                reps = list(directory.glob(f"{sample}.*.fragments.bed"))
                self.assertGreaterEqual(len(reps), 2)
                expected = Counter()
                for rep in reps:
                    expected.update(fragments(rep))
                self.assertEqual(Counter(fragments(pool)), expected)

    def test_actual_bigwig_values_match_fragments(self):
        checked = 0
        for fixture in FIXTURES:
            if fixture == "tiny_broad":  # explicitly disables all BigWigs
                continue
            directory = ROOT / f"out_{fixture}"
            beds = list((directory / "fragments").rglob("*.fragments.bed"))
            self.assertTrue(beds)
            for bed in beds:
                identifier = bed.name.removesuffix(".fragments.bed")
                if bed.parent.name == "pooled":
                    sample = identifier.removesuffix(".pooled")
                    base = directory / f"bigwig/pooled/{sample}/{identifier}"
                else:
                    base = directory / f"bigwig/{identifier}/{identifier}"
                records = fragments(bed)
                # The seven tiny configs use bin_size=5; consensus uses CPM.
                factor = (1e6 if fixture == "tiny_consensus" else 1e9 / 5) / len(records)
                tracks = [(Path(str(base) + ".normalized.bw"), factor)]
                raw = Path(str(base) + ".raw.bw")
                if raw.exists():
                    tracks.append((raw, 1.0))
                spike = Path(str(base) + ".spikein_normalized.bw")
                if spike.exists():
                    scale = table(directory / f"qc/spikein/{identifier}/{identifier}.spikein_scale.tsv")[0]
                    tracks.append((spike, float(scale["scale_factor"])))
                for path, multiplier in tracks:
                    with pyBigWig.open(str(path)) as bw:
                        for chrom, length in bw.chroms().items():
                            expected = [0] * length
                            for contig, start, end in records:
                                if contig == chrom:
                                    for pos in range(start // 5 * 5, min(length, ((end - 1) // 5 + 1) * 5)):
                                        expected[pos] += 1
                            for actual, value in zip(bw.values(chrom, 0, length), expected):
                                value *= multiplier
                                self.assertAlmostEqual(actual, value, delta=max(1e-5, abs(value) * 2e-6),
                                                       msg=str(path))
                    checked += 1
        self.assertGreater(checked, 20)

    def test_consensus_preserves_candidates_and_supporter_names(self):
        directory = ROOT / "out_tiny_consensus/macs3/consensus"
        rows = table(directory / "Tiny_H3K4me3.support.tsv")
        accepted = [r for r in rows if r["passes_support"] == "true"]
        self.assertTrue(accepted)
        expected = [(r["chrom"], int(r["start"]), int(r["end"])) for r in accepted]
        self.assertEqual(fragments(directory / "Tiny_H3K4me3.consensus.bed"), expected)
        for row in accepted:
            self.assertEqual(json.loads(row["support_replicates"]), ["rep1", "rep2"])
            self.assertEqual(row["match_rule"], "summit_distance")
            source = json.loads(row["source_peak_fields_json"])
            self.assertEqual(source[:3], [row["chrom"], row["start"], row["end"]])
        with (directory / "Tiny_H3K4me3.consensus.bed").open() as handle:
            self.assertTrue(all(len(line.rstrip().split("\t")) == 6 for line in handle))

    def test_broad_core_survives_rejected_large_candidate(self):
        directory = ROOT / "out_tiny_consensus/macs3/consensus"
        candidates = {r["candidate_id"]: r for r in table(directory / "Tiny_H3K27me3.support.tsv")}
        cores = table(directory / "Tiny_H3K27me3.core.support.tsv")
        self.assertTrue(cores)
        for core in cores:
            self.assertEqual(core["support_count"], "3")
            self.assertEqual(json.loads(core["support_replicates"]), ["rep1", "rep2", "rep3"])
            for identifier in json.loads(core["candidate_ids"]):
                parent = candidates[identifier]
                self.assertEqual(parent["passes_support"], "false")
                self.assertGreaterEqual(int(core["start"]), int(parent["start"]))
                self.assertLessEqual(int(core["end"]), int(parent["end"]))
        self.assertEqual(fragments(directory / "Tiny_H3K27me3.core.bed"),
                         [(r["chrom"], int(r["start"]), int(r["end"])) for r in cores])

    def test_tiny_idr_is_explicitly_insufficient(self):
        directory = ROOT / "out_tiny_consensus/macs3/idr/Tiny_H3K4me3"
        rows = table(directory / "summary.tsv")
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual((row["replicate_a"], row["replicate_b"]), ("rep1", "rep2"))
        self.assertEqual(row["status"], "insufficient_peaks")
        self.assertEqual(row["passing_peaks"], "0")
        self.assertEqual((directory / "pairs" / row["peak_file"]).stat().st_size, 0)

    def test_spikein_and_multiqc(self):
        directory = ROOT / "out_tiny_spikein"
        row = table(directory / "qc/spikein/spikein_scale_factors.tsv")[0]
        self.assertEqual(row["status"], "OK")
        self.assertEqual(int(row["spikein_reads"]), 2 * int(row["spikein_fragments"]))
        self.assertAlmostEqual(float(row["scale_factor"]),
                               float(row["scale_to"]) / int(row["spikein_fragments"]), places=8)
        self.assertGreater((directory / "qc/multiqc/multiqc_report.html").stat().st_size, 10000)


if __name__ == "__main__":
    unittest.main()
