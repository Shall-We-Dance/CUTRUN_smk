"""Snakemake entry point for pooled-candidate consensus."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir)))
from consensus import write_consensus

for output in snakemake.output:
    Path(str(output)).parent.mkdir(parents=True, exist_ok=True)
Path(str(snakemake.log[0])).parent.mkdir(parents=True, exist_ok=True)
names = list(snakemake.params.replicates)
paths = list(snakemake.input.peaks)
if len(names) != len(paths) or len(set(names)) != len(names):
    raise ValueError("Replicate identities must map one-to-one to peak files")
write_consensus(
    str(snakemake.input.candidates), dict(zip(names, paths)),
    snakemake.params.mode, dict(snakemake.params.options), snakemake.wildcards.sample,
    str(snakemake.output.consensus), str(snakemake.output.support), str(snakemake.output.tsv),
    str(snakemake.output.core) if "core" in snakemake.output.keys() else None,
    str(snakemake.output.core_support) if "core_support" in snakemake.output.keys() else None,
)
Path(str(snakemake.log[0])).write_text("Pooled candidates evaluated; see support TSV.\n")
