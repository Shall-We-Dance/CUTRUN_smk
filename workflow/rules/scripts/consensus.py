"""Pooled-candidate support and basewise cores (0-based, half-open intervals)."""

import bisect
import csv
import json
from collections import defaultdict
from dataclasses import dataclass


@dataclass(frozen=True)
class Peak:
    chrom: str
    start: int
    end: int
    fields: tuple = ()

    @property
    def name(self):
        return self.fields[3] if len(self.fields) > 3 else "."

    @property
    def summit(self):
        if len(self.fields) < 10:
            raise ValueError("summit_distance matching requires narrowPeak with valid summits")
        offset = int(self.fields[9])
        if not 0 <= offset < self.end - self.start:
            raise ValueError(f"Invalid summit for {self.chrom}:{self.start}-{self.end}")
        return self.start + offset


def read_peaks(path):
    peaks = []
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = tuple(line.rstrip("\n").split("\t"))
            try:
                start, end = int(fields[1]), int(fields[2])
                if start < 0 or end <= start:
                    raise ValueError("invalid interval")
                peaks.append(Peak(fields[0], start, end, fields))
            except (IndexError, ValueError) as exc:
                raise ValueError(f"{path}:{number}: invalid BED/peak interval") from exc
    return sorted(peaks, key=lambda p: (p.chrom, p.start, p.end, p.fields))


class PeakIndex:
    def __init__(self, peaks):
        self.rows = defaultdict(list)
        for peak in peaks:
            self.rows[peak.chrom].append(peak)
        self.starts, self.max_ends = {}, {}
        for chrom, rows in self.rows.items():
            rows.sort(key=lambda p: (p.start, p.end))
            self.starts[chrom] = [p.start for p in rows]
            high, ends = 0, []
            for peak in rows:
                high = max(high, peak.end)
                ends.append(high)
            self.max_ends[chrom] = ends

    def query(self, chrom, start, end):
        rows = self.rows.get(chrom, [])
        if not rows:
            return []
        left = bisect.bisect_right(self.max_ends[chrom], start)
        right = bisect.bisect_left(self.starts[chrom], end)
        return [p for p in rows[left:right] if p.end > start]


def merged_intervals(peaks):
    result = []
    for peak in sorted(peaks, key=lambda p: (p.chrom, p.start, p.end)):
        if result and result[-1].chrom == peak.chrom and peak.start <= result[-1].end:
            previous = result[-1]
            result[-1] = Peak(peak.chrom, previous.start, max(previous.end, peak.end))
        else:
            result.append(Peak(peak.chrom, peak.start, peak.end))
    return result


def covered_bases(candidate, overlaps):
    clipped = [Peak(candidate.chrom, max(candidate.start, p.start),
                    min(candidate.end, p.end)) for p in overlaps]
    return sum(p.end - p.start for p in merged_intervals(clipped))


def required_support(value, n):
    if n < 2:
        raise ValueError("Consensus requires at least two independent replicates")
    if value == "majority":
        return n // 2 + 1
    if value == "all":
        return n
    if isinstance(value, bool) or not isinstance(value, int) or not 2 <= value <= n:
        raise ValueError(f"min_replicates must be majority, all, or an integer in [2, {n}]")
    return value


def validate_options(options, n):
    required_support(options.get("min_replicates", "majority"), n)
    narrow = options.get("narrow", {})
    broad = options.get("broad", {})
    if narrow.get("match", "reciprocal_overlap") not in {"reciprocal_overlap", "summit_distance"}:
        raise ValueError("narrow.match must be reciprocal_overlap or summit_distance")
    for label, value in [("narrow.min_overlap_fraction", narrow.get("min_overlap_fraction", 0.5)),
                         ("broad.min_coverage_fraction", broad.get("min_coverage_fraction", 0.5))]:
        if isinstance(value, bool) or not 0 < float(value) <= 1:
            raise ValueError(f"{label} must be in (0, 1]")
    for label, value, minimum in [("narrow.summit_distance", narrow.get("summit_distance", 100), 0),
                                  ("broad.min_core_length", broad.get("min_core_length", 1), 1)]:
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f"{label} must be an integer >= {minimum}")


def candidate_support(candidate, index, mode, options):
    overlaps = index.query(candidate.chrom, candidate.start, candidate.end)
    fraction = covered_bases(candidate, overlaps) / (candidate.end - candidate.start)
    if mode == "broad":
        return fraction >= options.get("min_coverage_fraction", 0.5), fraction, overlaps
    if options.get("match", "reciprocal_overlap") == "summit_distance":
        summit = candidate.summit
        distance = options.get("summit_distance", 100)
        nearby = index.query(candidate.chrom, max(0, summit - distance), summit + distance + 1)
        matched = [p for p in nearby if abs(p.summit - summit) <= distance]
    else:
        threshold = options.get("min_overlap_fraction", 0.5)
        matched = []
        for peak in overlaps:
            overlap = min(candidate.end, peak.end) - max(candidate.start, peak.start)
            if (overlap / (candidate.end - candidate.start) >= threshold
                    and overlap / (peak.end - peak.start) >= threshold):
                matched.append(peak)
    return bool(matched), fraction, matched


def supported_candidates(candidates, replicate_peaks, mode, options):
    validate_options(options, len(replicate_peaks))
    minimum = required_support(options.get("min_replicates", "majority"), len(replicate_peaks))
    indexes = {name: PeakIndex(peaks) for name, peaks in sorted(replicate_peaks.items())}
    rows = []
    for i, candidate in enumerate(candidates, 1):
        details, supporters = {}, []
        for name, index in indexes.items():
            supported, fraction, matches = candidate_support(candidate, index, mode, options.get(mode, {}))
            details[name] = {"supported": supported, "candidate_coverage": round(fraction, 8),
                             "matched_peaks": [p.name for p in matches]}
            if supported:
                supporters.append(name)
        rows.append({"peak": candidate, "id": f"candidate_{i:07d}", "supporters": supporters,
                     "required": minimum, "passed": len(supporters) >= minimum, "details": details})
    return rows


def core_segments(rows, replicate_peaks, minimum, min_length=1):
    """Clip >=k basewise support to pooled candidates; never bridge gaps.

    Adjacent segments are joined ONLY if both the replicate and candidate
    identity sets are unchanged. Overlaps within a replicate contribute once.
    Core eligibility is independent of whole-candidate coverage: pooling can
    broaden a candidate while a smaller region is genuinely reproducible.
    """
    events = defaultdict(lambda: defaultdict(list))
    for name, peaks in replicate_peaks.items():
        for peak in merged_intervals(peaks):
            events[peak.chrom][peak.start].append(("rep", name, 1))
            events[peak.chrom][peak.end].append(("rep", name, -1))
    for row in rows:
        peak = row["peak"]
        events[peak.chrom][peak.start].append(("candidate", row["id"], 1))
        events[peak.chrom][peak.end].append(("candidate", row["id"], -1))
    result = []
    for chrom, positions in sorted(events.items()):
        active = {"rep": set(), "candidate": set()}
        previous = None
        for pos, changes in sorted(positions.items()):
            if previous is not None and pos > previous and len(active["rep"]) >= minimum and active["candidate"]:
                reps, candidates = tuple(sorted(active["rep"])), tuple(sorted(active["candidate"]))
                segment = (chrom, previous, pos, reps, candidates)
                if result and result[-1][0] == chrom and result[-1][2] == previous and result[-1][3:] == segment[3:]:
                    old = result[-1]
                    result[-1] = (chrom, old[1], pos, reps, candidates)
                else:
                    result.append(segment)
            for kind, name, delta in changes:
                if delta > 0:
                    active[kind].add(name)
                else:
                    active[kind].discard(name)
            previous = pos
    return [seg for seg in result if seg[2] - seg[1] >= min_length]


def pairwise_metrics(replicate_peaks):
    names = sorted(replicate_peaks)
    metrics = []
    for i, a in enumerate(names):
        aa = replicate_peaks[a]
        for b in names[i + 1:]:
            bb = replicate_peaks[b]
            ia, ib = PeakIndex(aa), PeakIndex(bb)
            metrics.append((f"{a}_peaks_overlapping_{b}", sum(bool(ib.query(p.chrom, p.start, p.end)) for p in aa)))
            metrics.append((f"{b}_peaks_overlapping_{a}", sum(bool(ia.query(p.chrom, p.start, p.end)) for p in bb)))
            intersection = sum(covered_bases(p, ib.query(p.chrom, p.start, p.end)) for p in merged_intervals(aa))
            union = sum(p.end - p.start for p in merged_intervals(aa + bb))
            metrics.append((f"{a}_vs_{b}_bp_jaccard", intersection / union if union else "NA"))
    return metrics


def write_consensus(candidates_path, replicate_paths, mode, options, sample,
                    bed_path, support_path, metrics_path, core_path=None, core_support_path=None):
    if mode not in {"narrow", "broad"}:
        raise ValueError("Peak mode must be narrow or broad")
    replicates = {name: read_peaks(path) for name, path in replicate_paths.items()}
    candidates = read_peaks(candidates_path)
    rows = supported_candidates(candidates, replicates, mode, options)
    minimum = required_support(options.get("min_replicates", "majority"), len(replicates))
    with open(bed_path, "w") as bed, open(support_path, "w") as support:
        writer = csv.writer(support, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "candidate_id", "source_peak", "support_count",
                         "support_replicates", "required_support", "passes_support", "peak_mode",
                         "match_rule", "replicate_support_json", "source_peak_fields_json"])
        for row in rows:
            peak = row["peak"]
            if row["passed"]:
                bed.write(f"{peak.chrom}\t{peak.start}\t{peak.end}\t{row['id']}\t0\t.\n")
            writer.writerow([peak.chrom, peak.start, peak.end, row["id"], peak.name,
                             len(row["supporters"]), json.dumps(row["supporters"]), minimum,
                             str(row["passed"]).lower(), mode,
                             "candidate_coverage" if mode == "broad" else options.get("narrow", {}).get("match", "reciprocal_overlap"),
                             json.dumps(row["details"], sort_keys=True), json.dumps(peak.fields)])
    cores = []
    if mode == "broad":
        if not core_path or not core_support_path:
            raise ValueError("Broad consensus requires both core output paths")
        cores = core_segments(rows, replicates, minimum, options.get("broad", {}).get("min_core_length", 1))
        with open(core_path, "w") as bed, open(core_support_path, "w") as support:
            writer = csv.writer(support, delimiter="\t", lineterminator="\n")
            writer.writerow(["chrom", "start", "end", "core_id", "support_count", "support_replicates", "candidate_ids"])
            for i, (chrom, start, end, reps, ids) in enumerate(cores, 1):
                name = f"core_{i:07d}"
                bed.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t.\n")
                writer.writerow([chrom, start, end, name, len(reps), json.dumps(reps), json.dumps(ids)])
    metrics = [("status", "ok"), ("num_replicates", len(replicates)),
               ("required_support", minimum), ("pooled_candidates", len(candidates)),
               ("supported_candidates", sum(r["passed"] for r in rows))]
    metrics += [(f"{name}_peaks", len(peaks)) for name, peaks in sorted(replicates.items())]
    metrics += pairwise_metrics(replicates)
    if mode == "broad":
        metrics += [("core_segments", len(cores)), ("core_bases", sum(s[2] - s[1] for s in cores))]
    with open(metrics_path, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "metric", "value"])
        writer.writerows((sample, key, value) for key, value in metrics)
