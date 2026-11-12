"""Hand-computed and adversarial replicate-consensus regressions."""
import csv
import json
from pathlib import Path
import random
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/rules/scripts"))
from consensus import (Peak, PeakIndex, candidate_support, core_segments, merged_intervals,
                       pairwise_metrics, required_support, supported_candidates, validate_options,
                       write_consensus)


def peak(start, end, summit=None, name="p", chrom="chr1"):
    fields = (chrom, str(start), str(end), name, "10", ".", "20", "5", "4", str(summit if summit is not None else (end-start)//2))
    return Peak(chrom, start, end, fields)


class ConsensusTests(unittest.TestCase):
    def test_thresholds(self):
        self.assertEqual([required_support("majority", n) for n in [2, 3, 4, 5]], [2, 2, 3, 3])
        self.assertEqual(required_support("all", 4), 4)
        for value in [0, 1, 4, 2.0, True, "2"]:
            with self.assertRaises(ValueError):
                required_support(value, 3)
        with self.assertRaises(ValueError):
            required_support("majority", 1)

    def test_one_bp_and_short_containment_are_not_reciprocal_support(self):
        candidate = peak(100, 200)
        for other in [peak(199, 300), peak(120, 130), peak(0, 1000)]:
            self.assertFalse(candidate_support(candidate, PeakIndex([other]), "narrow", {})[0])
        self.assertTrue(candidate_support(candidate, PeakIndex([peak(150, 250)]), "narrow", {})[0])

    def test_summit_matching_is_explicit(self):
        candidate = peak(0, 20, 10)
        options = {"match": "summit_distance", "summit_distance": 20}
        self.assertTrue(candidate_support(candidate, PeakIndex([peak(20, 40, 10)]), "narrow", options)[0])
        self.assertFalse(candidate_support(candidate, PeakIndex([peak(21, 41, 10)]), "narrow", options)[0])
        with self.assertRaises(ValueError):
            candidate_support(Peak("chr1", 0, 20), PeakIndex([]), "narrow", options)

    def test_broad_union_coverage_not_sum_of_overlaps(self):
        candidate = peak(0, 100)
        supported, fraction, _ = candidate_support(candidate, PeakIndex([peak(0, 30), peak(20, 40)]), "broad", {})
        self.assertFalse(supported)
        self.assertEqual(fraction, 0.4)
        supported, fraction, _ = candidate_support(candidate, PeakIndex([peak(0, 30), peak(70, 100)]), "broad", {})
        self.assertTrue(supported)
        self.assertEqual(fraction, 0.6)

    def test_pooled_boundaries_are_never_expanded(self):
        rows = supported_candidates([peak(100, 200)], {"r1": [peak(90, 210)], "r2": [peak(95, 205)]}, "narrow", {})
        self.assertTrue(rows[0]["passed"])
        self.assertEqual((rows[0]["peak"].start, rows[0]["peak"].end), (100, 200))

    def test_two_of_three(self):
        reps = {"good_b": [peak(0, 100)], "shallow": [], "good_a": [peak(0, 100)]}
        row = supported_candidates([peak(0, 100)], reps, "narrow", {})[0]
        self.assertEqual(row["supporters"], ["good_a", "good_b"])
        self.assertTrue(row["passed"])
        self.assertFalse(supported_candidates([peak(0, 100)], reps, "narrow", {"min_replicates": "all"})[0]["passed"])

    def test_long_bridge_has_no_three_way_core(self):
        reps = {"long": [peak(0, 1000)], "left": [peak(0, 100)], "right": [peak(900, 1000)]}
        options = {"min_replicates": "all", "broad": {"min_coverage_fraction": 0.1}}
        rows = supported_candidates([peak(0, 1000)], reps, "broad", options)
        self.assertTrue(rows[0]["passed"])
        self.assertEqual(core_segments(rows, reps, 3), [])
        cores = core_segments(rows, reps, 2)
        self.assertEqual([(c[1], c[2]) for c in cores], [(0, 100), (900, 1000)])

    def test_same_count_different_supporters_split(self):
        reps = {"a": [peak(0, 100)], "b": [peak(0, 50)], "c": [peak(50, 100)]}
        rows = supported_candidates([peak(0, 100)], reps, "broad", {})
        cores = core_segments(rows, reps, 2)
        self.assertEqual([(c[1], c[2], c[3]) for c in cores], [(0, 50, ("a", "b")), (50, 100, ("a", "c"))])

    def test_overlapping_peaks_in_one_replicate_count_once(self):
        reps = {"a": [peak(0, 100), peak(0, 100), peak(20, 90)], "b": []}
        rows = [{"peak": peak(0, 100), "id": "candidate", "passed": True}]
        self.assertEqual(core_segments(rows, reps, 2), [])

    def test_core_clipping_and_minimum_length(self):
        reps = {"a": [peak(0, 100)], "b": [peak(0, 100)]}
        rows = supported_candidates([peak(30, 50)], reps, "broad", {})
        self.assertEqual([(s[1], s[2]) for s in core_segments(rows, reps, 2)], [(30, 50)])
        self.assertEqual(core_segments(rows, reps, 2, 21), [])

    def test_core_survives_whole_candidate_coverage_failure(self):
        reps = {"a": [peak(40, 60)], "b": [peak(40, 60)]}
        rows = supported_candidates([peak(0, 100)], reps, "broad", {})
        self.assertFalse(rows[0]["passed"])
        self.assertEqual([(s[1], s[2]) for s in core_segments(rows, reps, 2)], [(40, 60)])

    def test_empty_metrics_and_order_independence(self):
        self.assertIn(("a_vs_b_bp_jaccard", "NA"), pairwise_metrics({"b": [], "a": []}))
        a = {"a": [peak(0, 100)], "b": [peak(50, 150)], "c": []}
        b = dict(reversed(list(a.items())))
        self.assertEqual(supported_candidates([peak(40, 110)], a, "broad", {}),
                         supported_candidates([peak(40, 110)], b, "broad", {}))

    def test_support_table_includes_rejected_and_source_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            for name, entries in {"pooled": [peak(0, 100), peak(500, 600)], "named_a": [peak(0, 100)], "named_b": [peak(0, 100)]}.items():
                (d / name).write_text("".join("\t".join(p.fields) + "\n" for p in entries))
            write_consensus(d/"pooled", {"named_b": d/"named_b", "named_a": d/"named_a"}, "narrow", {}, "test",
                            d/"out.bed", d/"support.tsv", d/"metrics.tsv")
            with (d/"support.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([r["passes_support"] for r in rows], ["true", "false"])
            self.assertEqual(json.loads(rows[0]["support_replicates"]), ["named_a", "named_b"])
            self.assertEqual(len(json.loads(rows[0]["source_peak_fields_json"])), 10)
            self.assertEqual(len((d/"out.bed").read_text().strip().split("\t")), 6)

    def test_empty_outputs_are_valid(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            (d/"empty").write_text("")
            write_consensus(d/"empty", {"a": d/"empty", "b": d/"empty"}, "broad", {}, "test",
                            d/"out.bed", d/"support.tsv", d/"metrics.tsv", d/"core.bed", d/"core.tsv")
            self.assertEqual((d/"out.bed").read_text(), "")
            self.assertEqual((d/"core.bed").read_text(), "")
            self.assertEqual(len((d/"support.tsv").read_text().splitlines()), 1)

    def test_nested_interval_index_matches_brute_force(self):
        rng = random.Random(4)
        peaks = [peak(s, s + rng.randint(1, 100)) for s in [rng.randrange(200) for _ in range(100)]]
        index = PeakIndex(peaks)
        for start in range(200):
            actual = set(index.query("chr1", start, start + 3))
            expected = {p for p in peaks if p.start < start + 3 and p.end > start}
            self.assertEqual(actual, expected)

    def test_invalid_options(self):
        for options in [{"narrow": {"match": "any"}}, {"broad": {"min_coverage_fraction": 0}},
                        {"narrow": {"min_overlap_fraction": 1.1}}, {"broad": {"min_core_length": 0}}]:
            with self.assertRaises(ValueError):
                validate_options(options, 3)


if __name__ == "__main__":
    unittest.main()
