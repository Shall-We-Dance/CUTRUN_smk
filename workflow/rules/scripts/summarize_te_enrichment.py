#!/usr/bin/env python3
"""Family-level TE enrichment of ChIP samples versus matched input.

Replicate matching:
  - Preferred: `control_replicates` maps every ChIP replicate to a named
    control replicate (e.g. `{rep1: input_a, rep2: input_b}`).
  - Legacy: a scalar `control_replicate` applies to every ChIP replicate.
  - If `control_replicate: pooled` (or the control has a single replicate),
    the pooled control is used for every ChIP replicate and the columns are
    labeled `vs_pooled_input`.
  - Otherwise, the ChIP replicate's own replicate name is matched to the
    control (rep1 vs rep1, rep2 vs rep2).

Family classification (from the TE annotation family column):
  curated  : name carries a Dfam-style accession (Name|DFxxxxx) or matches a
             known non-novel pattern
  novel    : EDTA novel IDs (`TE_#####` or `TE_#####_suffix`) or
             coordinate-style labels (2L:start..end, bare integers)
  unknown  : anything else (kept separate; never called curated)

Per family and per sample:
  - per-replicate 1/N weighted counts (ChIP and matched input)
  - pooled weighted counts and RPM (= count / (templates / 1e6))
  - pooled enrichment_ratio / log2_enrichment with pseudocount 1.0 RPM
  - per-replicate log2 enrichment (paired, matched input) and
    direction consistency (all replicates same sign beyond |log2|>0.1)

Outputs (all under {OUTDIR}/te/):
  te_family_enrichment.tsv             all families (pooled + per-rep)
  te_family_enrichment.curated.tsv     curated/Dfam-mapped families only
  te_family_enrichment.novel.tsv       EDTA novel families only
  te_family_enrichment.unknown.tsv     unclassified families
  te_family_enrichment.subclass.tsv    aggregation by subfamily column
"""

import math
import re
from collections import defaultdict
from pathlib import Path

PSEUDOCOUNT = 1.0
family_input = snakemake.input.family_counts
out_tsv = snakemake.output.tsv
curated_tsv = str(out_tsv).replace(".tsv", ".curated.tsv")
novel_tsv = str(out_tsv).replace(".tsv", ".novel.tsv")
unknown_tsv = str(out_tsv).replace(".tsv", ".unknown.tsv")
subclass_tsv = str(out_tsv).replace(".tsv", ".subclass.tsv")

NOVEL_RE = re.compile(
    r"^TE_\d+($|_)|^\d+[.:]|^[234XxYyMm]L?R?:|^2L:|^2R:|^3L:|^3R:|^X:|^Y:|^MT:|^chr",
    re.IGNORECASE,
)

# Dfam accession lookup: build {name: accession} from the TE library metadata
# (dfam4-dmel-TE-metadata.tsv). A family is "curated" only when its name is
# backed by a Dfam accession in this table (or carries an explicit DF accession
# suffix). Novel IDs and anything else are novel/unknown.
DFAM_NAMES = {}
try:
    _dfam_meta = snakemake.params.dfam_metadata
except Exception:
    _dfam_meta = None
if _dfam_meta and Path(_dfam_meta).exists():
    with open(_dfam_meta) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                DFAM_NAMES[f[1]] = f[0]


def classify_family(family):
    if NOVEL_RE.match(family):
        return "novel"
    base = family.split("|")[0].strip()
    if base in DFAM_NAMES or "|DF" in family:
        return "curated"
    # curated-looking Drosophila TE names without a Dfam entry are still
    # not verifiable -> unknown (not "curated")
    return "unknown"


def read_family_counts(path):
    """Return {family: weighted_count}, {family: subfamily}, total_templates."""
    counts = {}
    subfamilies = {}
    total_templates = None
    overlap_tsv = str(path).replace(".te_family_counts.tsv", ".te_overlap.tsv")
    if Path(overlap_tsv).exists():
        with open(overlap_tsv) as fh:
            for line in fh:
                metric, value = line.rstrip("\n").split("\t", 1)
                if metric == "total_templates":
                    total_templates = int(value)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            values = line.rstrip("\n").split("\t")
            row = dict(zip(header, values))
            counts[row["family"]] = float(row["weighted_count"])
            subfamilies[row["family"]] = row.get("subfamily", "")
    return counts, subfamilies, total_templates


def replicate_family_counts(sample):
    """Return dict replicate_name -> (counts, subfamilies, total_templates)."""
    reps = {}
    prefix = f"{sample}."
    implicit = bool(
        snakemake.config["samples"][sample].get("_implicit_single_replicate", False)
    )
    for path in family_input:
        name = Path(path).name
        if "family_counts" not in name:
            continue
        if implicit:
            if Path(path).parent.name != sample:
                continue
            counts, subfams, n = read_family_counts(path)
            reps["rep1"] = (counts, subfams, n or 0)
            continue
        if not name.startswith(prefix):
            continue
        rep = Path(path).parent.name
        if rep.startswith(prefix):
            rep = rep[len(prefix):]
        counts, subfams, n = read_family_counts(path)
        reps[rep] = (counts, subfams, n or 0)
    return reps


def rpm(counts, depth):
    return counts / (depth / 1e6) if depth else 0.0


def control_mapping(sample_cfg, chip_reps):
    """Return dict chip_rep -> control_rep (or "pooled")."""
    control = sample_cfg.get("control")
    if not control:
        return {}
    explicit_map = sample_cfg.get("control_replicates", {}) or {}
    explicit = sample_cfg.get("control_replicate")
    mapping = {}
    for rep in chip_reps:
        if rep in explicit_map:
            mapping[rep] = str(explicit_map[rep])
        elif explicit is not None:
            mapping[rep] = str(explicit)
        else:
            mapping[rep] = rep
    return mapping


rows = []
subclass_rows = []
for sample in snakemake.config["samples"]:
    sample_cfg = snakemake.config["samples"][sample]
    control = sample_cfg.get("control")
    if not control:
        continue
    chip_reps = replicate_family_counts(sample)
    input_reps = replicate_family_counts(control)
    if not chip_reps:
        continue
    mapping = control_mapping(sample_cfg, chip_reps)
    missing = sorted({ctrl for ctrl in mapping.values()
                      if ctrl != "pooled" and ctrl not in input_reps})
    if missing:
        raise ValueError(
            f"{sample}: mapped control replicate(s) not found in {control}: "
            + ", ".join(missing)
        )

    # pooled totals
    chip_depth = sum(r[2] for r in chip_reps.values())
    input_depth = sum(r[2] for r in input_reps.values())
    chip_pooled = defaultdict(float)
    for counts, _, _ in chip_reps.values():
        for fam, c in counts.items():
            chip_pooled[fam] += c
    input_pooled = defaultdict(float)
    for counts, _, _ in input_reps.values():
        for fam, c in counts.items():
            input_pooled[fam] += c

    if chip_depth <= 0 or input_depth <= 0:
        continue

    all_families = set(chip_pooled) | set(input_pooled)
    for family in sorted(all_families):
        chip_rpm = rpm(chip_pooled.get(family, 0.0), chip_depth)
        input_rpm = rpm(input_pooled.get(family, 0.0), input_depth)
        ratio = (chip_rpm + PSEUDOCOUNT) / (input_rpm + PSEUDOCOUNT)
        log2_enrichment = math.log2(ratio)

        # per-replicate paired enrichment
        rep_log2s = []
        rep_details = []
        for rep, (counts, _, depth) in sorted(chip_reps.items()):
            ctrl = mapping.get(rep, "pooled")
            if ctrl == "pooled":
                c_counts, c_depth = input_pooled, input_depth
                label = "vs_pooled"
            else:
                c_counts = input_reps.get(ctrl, (defaultdict(float), {}, 0))[0]
                c_depth = input_reps.get(ctrl, ({}, {}, 0))[2]
                label = f"vs_{ctrl}"
            rc = rpm(counts.get(family, 0.0), depth)
            ri = rpm(c_counts.get(family, 0.0), c_depth)
            rep_log2s.append(math.log2((rc + PSEUDOCOUNT) / (ri + PSEUDOCOUNT)))
            rep_details.append(label)

        signs = [1 if v > 0 else -1 for v in rep_log2s if abs(v) > 0.1]
        if len(chip_reps) < 2:
            direction_consistent = "NA"
        elif not signs:
            direction_consistent = "NA"
        else:
            direction_consistent = (
                "TRUE"
                if len(signs) == len(chip_reps) and len(set(signs)) == 1
                else "FALSE"
            )

        rows.append(
            {
                "sample": sample,
                "control": control,
                "family": family,
                "family_class": classify_family(family),
                "chip_weighted_count": f"{chip_pooled.get(family, 0.0):.4f}",
                "chip_templates": chip_depth,
                "chip_rpm": f"{chip_rpm:.6f}",
                "input_weighted_count": f"{input_pooled.get(family, 0.0):.4f}",
                "input_templates": input_depth,
                "input_rpm": f"{input_rpm:.6f}",
                "enrichment_ratio": f"{ratio:.4f}",
                "log2_enrichment": f"{log2_enrichment:.4f}",
                "rep_log2_enrichment": ";".join(f"{v:.4f}" for v in rep_log2s),
                "rep_input_mapping": ";".join(rep_details),
                "direction_consistent": direction_consistent,
            }
        )

    # subclass aggregation
    sub_chip = defaultdict(float)
    sub_input = defaultdict(float)
    for _, (counts, subfams, _) in chip_reps.items():
        for fam, c in counts.items():
            sub_chip[subfams.get(fam, "") or "unknown"] += c
    for _, (counts, subfams, _) in input_reps.items():
        for fam, c in counts.items():
            sub_input[subfams.get(fam, "") or "unknown"] += c
    for sf in sorted(set(sub_chip) | set(sub_input)):
        cr = rpm(sub_chip.get(sf, 0.0), chip_depth)
        ir = rpm(sub_input.get(sf, 0.0), input_depth)
        ratio = (cr + PSEUDOCOUNT) / (ir + PSEUDOCOUNT)
        subclass_rows.append(
            {
                "sample": sample, "control": control, "subclass": sf,
                "chip_weighted_count": f"{sub_chip.get(sf, 0.0):.4f}",
                "chip_rpm": f"{cr:.6f}",
                "input_weighted_count": f"{sub_input.get(sf, 0.0):.4f}",
                "input_rpm": f"{ir:.6f}",
                "enrichment_ratio": f"{ratio:.4f}",
                "log2_enrichment": f"{math.log2(ratio):.4f}",
            }
        )

Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
header = [
    "sample", "control", "family", "family_class",
    "chip_weighted_count", "chip_templates", "chip_rpm",
    "input_weighted_count", "input_templates", "input_rpm",
    "enrichment_ratio", "log2_enrichment",
    "rep_log2_enrichment", "rep_input_mapping", "direction_consistent",
]


def write_tsv(path, subset_rows):
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in subset_rows:
            out.write("\t".join(str(row[c]) for c in header) + "\n")


write_tsv(out_tsv, rows)
write_tsv(curated_tsv, [r for r in rows if r["family_class"] == "curated"])
write_tsv(novel_tsv, [r for r in rows if r["family_class"] == "novel"])
write_tsv(unknown_tsv, [r for r in rows if r["family_class"] == "unknown"])

with open(subclass_tsv, "w") as out:
    sub_header = ["sample", "control", "subclass", "chip_weighted_count",
                  "chip_rpm", "input_weighted_count", "input_rpm",
                  "enrichment_ratio", "log2_enrichment"]
    out.write("\t".join(sub_header) + "\n")
    for row in sorted(subclass_rows, key=lambda r: (r["sample"], -float(r["log2_enrichment"]))):
        out.write("\t".join(str(row[c]) for c in sub_header) + "\n")
