"""IDR thresholds, low-information status, and optional real-tool regression."""
import csv
import math
import os
from pathlib import Path
import random
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/rules/scripts"))
from replicate_idr import filter_idr, run_pairwise

IDR = os.environ.get("IDR") or shutil.which("idr")


class IDRTests(unittest.TestCase):
    def test_global_idr_and_ten_column_output(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            prefix = ["chr1", "0", "100", "peak", "0", ".", "100", "10", "99", "50"]
            (d/"raw").write_text("\t".join(prefix + ["10", "1"]) + "\n" + "\t".join(prefix + ["0", "2"]) + "\n")
            self.assertEqual(filter_idr(d/"raw", d/"out", 0.05), 1)
            self.assertEqual(len((d/"out").read_text().strip().split("\t")), 10)

    def test_insufficient_peaks_are_not_successful_reproducibility(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            (d/"empty").write_text("")
            run_pairwise({"a": d/"empty", "b": d/"empty", "c": d/"empty"}, d/"empty", d/"pairs", d/"summary",
                         executable="/definitely/not/a/tool")
            with (d/"summary").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 3)
            self.assertEqual({r["status"] for r in rows}, {"insufficient_peaks"})

    @unittest.skipUnless(IDR, "Set IDR or activate the pinned IDR environment for real-tool testing")
    def test_real_idr(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            rng = random.Random(42)
            signals = []
            for replicate in range(2):
                values = []
                for i in range(400):
                    values.append(1000 - i * 2 + rng.gauss(0, 15) if i < 200 else rng.uniform(1, 500))
                signals.append(values)
            for name, values in [("a", signals[0]), ("b", signals[1]), ("pooled", [(a+b)/2 for a,b in zip(*signals)])]:
                (d/name).write_text("".join(
                    f"chr1\t{i*1000}\t{i*1000+200}\tp{i}\t0\t.\t{value:.5f}\t{value:.5f}\t-1\t100\n"
                    for i, value in enumerate(values)))
            run_pairwise({"a": d/"a", "b": d/"b"}, d/"pooled", d/"pairs", d/"summary", executable=IDR)
            with (d/"summary").open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["status"], "ok")
            self.assertGreater(int(row["passing_peaks"]), 0)
            self.assertLessEqual(int(row["passing_peaks"]), 400)


if __name__ == "__main__":
    unittest.main()
