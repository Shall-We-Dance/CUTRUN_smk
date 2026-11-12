"""Shared molecule selection, exact coverage, and cross-consumer regression tests."""
import csv
import importlib.util
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow/rules/scripts"
sys.path.insert(0, str(SCRIPTS))
from fragment_records import usable_fragment, reference_end, excluded_contig, bed_fragments
from fragment_coverage import chromosome_runs, normalization_factor, write_bigwig
from spikein_templates import count_templates

SAMTOOLS = os.environ.get("SAMTOOLS") or shutil.which("samtools")
HAS_BIGWIG = importlib.util.find_spec("pyBigWig") is not None


def record(name="read", flag=99, chrom="chr1", pos=11, cigar="10M", rg="a", mapq=42):
    return [name, str(flag), chrom, str(pos), str(mapq), cigar, "=", "31", "999",
            "A" * 10, "I" * 10, f"RG:Z:{rg}"]


def metric_table(path):
    with open(path) as handle:
        return dict(csv.reader(handle, delimiter="\t"))


class FragmentTests(unittest.TestCase):
    def setUp(self):
        self.sizes = {"chr1": 1000, "chr2": 500}
        self.pair = [record(), record(flag=147, pos=31)]

    def test_span_uses_both_alignments_not_tlen(self):
        self.assertEqual(usable_fragment(self.pair, "PE", 1804, self.sizes)[0], ("chr1", 10, 40, 2))

    def test_orphan_and_malformed_pair(self):
        for group in [self.pair[:1], self.pair + [record()], [record(flag=67), record(flag=131, chrom="chr2")],
                      [record(flag=73), record(flag=133)]]:
            self.assertIsNone(usable_fragment(group, "PE", 1804, self.sizes)[0])

    def test_secondary_and_supplementary_do_not_create_molecules(self):
        group = self.pair + [record(flag=355), record(flag=99 | 2048)]
        self.assertEqual(usable_fragment(group, "PE", 1804, self.sizes)[0], ("chr1", 10, 40, 2))

    def test_flags_apply_to_both_mates(self):
        for bit in [4, 8, 512, 1024]:
            group = [record(), record(flag=147 | bit, pos=31)]
            self.assertIsNone(usable_fragment(group, "PE", 1804, self.sizes)[0])
        self.assertIsNotNone(usable_fragment([record(flag=99 | 1024), record(flag=147 | 1024)],
                                            "PE", 4, self.sizes)[0])

    def test_se_extension_and_clipping(self):
        self.assertEqual(usable_fragment([record(flag=0, pos=991)], "SE", 1804, self.sizes, 100)[0], ("chr1", 990, 1000, 1))
        self.assertEqual(usable_fragment([record(flag=16, pos=11)], "SE", 1804, self.sizes, 100)[0], ("chr1", 0, 20, 1))

    def test_cigar_and_contig_filters(self):
        self.assertEqual(reference_end(record(cigar="2S5M3D5M")), 23)
        self.assertIsNone(usable_fragment(self.pair, "PE", 1804, {"chr2": 500})[0])
        for name in ["chrM", "MT", "M", "chr1_alt", "scaffold42"]:
            self.assertTrue(excluded_contig(name))
        self.assertFalse(excluded_contig("2L"))

    def test_normalization_denominators_are_fragments(self):
        self.assertEqual(normalization_factor("CPM", 2, 100, 12, 1000, 10), 500000)
        self.assertEqual(normalization_factor("RPKM", 2, 100, 12, 1000, 10), 50000000)
        self.assertEqual(normalization_factor("RPGC", 2, 100, 12, 1000, 10), 10)
        self.assertAlmostEqual(normalization_factor("BPM", 2, 100, 12, 1000, 10), 1e6 / 12)
        self.assertEqual(normalization_factor("CPM", 0, 0, 0, 1000, 10), 0)

    def test_streaming_coverage_matches_brute_force(self):
        rng = random.Random(13)
        for bin_size in [1, 7, 10, 37]:
            fragments = []
            for _ in range(300):
                start = rng.randrange(1000)
                end = rng.randrange(start + 1, 1002)
                fragments.append(("chr1", start, end))
            fragments.sort()
            expected = [0] * 1001
            for _, start, end in fragments:
                lo, hi = start // bin_size * bin_size, min(1001, ((end - 1) // bin_size + 1) * bin_size)
                for base in range(lo, hi):
                    expected[base] += 1
            actual = [None] * 1001
            for start, end, value in chromosome_runs(fragments, 1001, bin_size):
                actual[start:end] = [value] * (end - start)
            self.assertEqual(actual, expected)

    @unittest.skipUnless(HAS_BIGWIG, "Requires pinned pyBigWig environment")
    def test_bigwig_values_and_empty_track(self):
        import pyBigWig
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            (d/"fragments").write_text("chr1\t0\t20\nchr1\t10\t30\n")
            count, factor = write_bigwig(d/"fragments", {"chr1": 35, "chr2": 5}, d/"out.bw", 10, "CPM")
            self.assertEqual(count, 2)
            with pyBigWig.open(str(d/"out.bw")) as bw:
                self.assertEqual(bw.values("chr1", 0, 35), [500000.0]*10 + [1000000.0]*10 + [500000.0]*10 + [0.0]*5)
                self.assertEqual(bw.values("chr2", 0, 5), [0.0]*5)
            (d/"empty").write_text("")
            write_bigwig(d/"empty", {"chr1": 35}, d/"empty.bw", 10, "RPKM")
            with pyBigWig.open(str(d/"empty.bw")) as bw:
                self.assertEqual(bw.values("chr1", 0, 35), [0.0]*35)

    @unittest.skipUnless(SAMTOOLS, "Requires samtools")
    def test_bam_to_shared_fragments_and_frip(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            (d/"sizes").write_text("chr1\t1000\nchr2\t500\n")
            # Same QNAME in two independent libraries must remain two molecules.
            rows = self.pair + [record(rg="b"), record(flag=147, pos=31, rg="b")]
            rows += [record(name="orphan")]
            rows += [record(name="secondary", flag=355), record(name="secondary", flag=403, pos=31)]
            rows += [record(name="lowq"), record(name="lowq", flag=147, pos=31, mapq=0)]
            rows += [record(name="cross"), record(name="cross", flag=147, chrom="chr2", pos=31)]
            # Mates do not touch [55,65), but the actual fragment span does.
            rows += [record(name="black", pos=41), record(name="black", flag=147, pos=81)]
            (d/"input.sam").write_text("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n"
                                      "@RG\tID:a\tSM:a\n@RG\tID:b\tSM:b\n" +
                                      "".join("\t".join(r) + "\n" for r in rows))
            subprocess.run([SAMTOOLS, "view", "-b", "-q", "20", "-o", str(d/"input.bam"), str(d/"input.sam")], check=True)
            (d/"blacklist").write_text("chr1\t55\t65\n")
            env = dict(os.environ, SAMTOOLS=SAMTOOLS)
            old_samtools = os.environ.get("SAMTOOLS")
            os.environ["SAMTOOLS"] = SAMTOOLS
            try:
                # Two independent identical pairs plus the unmasked long pair.
                self.assertEqual(count_templates(d/"input.bam", "PE", 1804), (6, 3))
            finally:
                if old_samtools is None:
                    os.environ.pop("SAMTOOLS", None)
                else:
                    os.environ["SAMTOOLS"] = old_samtools
            subprocess.run([sys.executable, str(SCRIPTS/"prepare_fragments.py"), "--bam", str(d/"input.bam"),
                            "--chrom-sizes", str(d/"sizes"), "--blacklist", str(d/"blacklist"),
                            "--read-type", "PE", "--sample", "synthetic", "--out", str(d/"fragments"),
                            "--metrics", str(d/"metrics"), "--lengths", str(d/"lengths")], check=True, env=env)
            self.assertEqual(list(bed_fragments(d/"fragments")), [("chr1", 10, 40), ("chr1", 10, 40)])
            metrics = metric_table(d/"metrics")
            self.assertEqual(metrics["usable_fragments_after_processing"], "2")
            self.assertEqual(metrics["usable_reads_after_processing"], "4")
            self.assertEqual(metrics["rejected_blacklist_span"], "1")
            self.assertEqual((d/"lengths").read_text(), "fragment_length\tcount\n30\t2\n")
            (d/"peaks").write_text("chr1\t20\t21\n")  # neither sequenced mate covers this base
            subprocess.run([sys.executable, str(SCRIPTS/"frip.py"), "--fragments", str(d/"fragments"),
                            "--peaks", str(d/"peaks"), "--sample", "synthetic", "--out", str(d/"frip")], check=True)
            frip = metric_table(d/"frip")
            self.assertEqual(frip["usable_templates"], "2")
            self.assertEqual(frip["templates_in_peaks"], "2")
            self.assertEqual(frip["read_frip"], "NA")


if __name__ == "__main__":
    unittest.main()
