# workflow/rules/common.smk
import os
from snakemake.io import temp

OUTDIR = config.get("output", {}).get("dir", "results")
SAMPLES = list(config.get("samples", {}).keys())
READ_TYPE = config.get("read_type", "PE").upper()

FILTER_BLACKLIST = bool(config.get("filter_blacklist", False))
REMOVE_DUPLICATES = bool(config.get("remove_duplicates", True))
KEEP_BAM = bool(config.get("keep_bam", False))
CUTRUN_QC_ENABLED = bool(config.get("cutrun_qc", {}).get("enabled", True))
PRESEQ_ENABLED = bool(config.get("cutrun_qc", {}).get("preseq", {}).get("enabled", True))
SPIKEIN_ENABLED = bool(config.get("spikein", {}).get("enabled", False))
SPIKEIN_BIGWIG_ENABLED = (
    SPIKEIN_ENABLED and bool(config.get("spikein", {}).get("generate_bigwig", True))
)
MACS3_CONFIG = config.get("macs3", {})
MACS3_ENABLED = bool(MACS3_CONFIG.get("enabled", True))


def blacklist_config():
    return config.get("blacklist", {})


def blacklist_path():
    blacklist = blacklist_config()
    if blacklist.get("path"):
        return blacklist["path"]
    if blacklist.get("url"):
        cache_dir = blacklist.get("cache_dir", "resources/blacklist")
        filename = os.path.basename(blacklist["url"])
        if filename.endswith(".gz"):
            filename = filename[:-3]
        return os.path.join(cache_dir, filename)
    return None


def unique_bam_path(sample):
    return f"{OUTDIR}/bowtie2/{sample}/{sample}.unique.bam"


def filtered_bam_path(sample):
    return f"{OUTDIR}/bowtie2/{sample}/{sample}.unique.filtered.bam"


def dedup_bam_path(sample):
    suffix = "unique.filtered" if FILTER_BLACKLIST else "unique"
    return f"{OUTDIR}/bowtie2/{sample}/{sample}.{suffix}.dedup.bam"


def duplicate_metrics_path(sample):
    return f"{OUTDIR}/bowtie2/{sample}/{sample}.dedup_metrics.txt"


def final_bam_path(sample):
    if REMOVE_DUPLICATES:
        return dedup_bam_path(sample)
    if FILTER_BLACKLIST:
        return filtered_bam_path(sample)
    return unique_bam_path(sample)


def final_bai_path(sample):
    return final_bam_path(sample) + ".bai"


def is_single_end():
    return READ_TYPE == "SE"


def aligned_bam_path(sample):
    suffix = "se" if is_single_end() else "pe"
    return f"{OUTDIR}/bowtie2/{sample}/{sample}_{suffix}.bam"


def maybe_temp(path):
    return path if KEEP_BAM else temp(path)


def aligned_flagstat_path(sample):
    suffix = "se" if is_single_end() else "pe"
    return f"{OUTDIR}/qc/bam_stats/{sample}/{sample}_{suffix}.flagstat.txt"


def unique_flagstat_path(sample):
    return f"{OUTDIR}/qc/bam_stats/{sample}/{sample}.unique.flagstat.txt"


def filtered_flagstat_path(sample):
    return f"{OUTDIR}/qc/bam_stats/{sample}/{sample}.unique.filtered.flagstat.txt"


def dedup_flagstat_path(sample):
    suffix = "unique.filtered" if FILTER_BLACKLIST else "unique"
    return f"{OUTDIR}/qc/bam_stats/{sample}/{sample}.{suffix}.dedup.flagstat.txt"


def bam_flagstat_paths(sample):
    paths = [aligned_flagstat_path(sample), unique_flagstat_path(sample)]
    if FILTER_BLACKLIST:
        paths.append(filtered_flagstat_path(sample))
    if REMOVE_DUPLICATES:
        paths.append(dedup_flagstat_path(sample))
    return paths


def bowtie2_log_path(sample):
    suffix = "se" if is_single_end() else "pe"
    return f"logs/bowtie2/{sample}_{suffix}.log"


def cutrun_usable_fragments_path(sample):
    return f"{OUTDIR}/qc/cutrun/{sample}/{sample}.usable_fragments.tsv"


def fragment_length_path(sample):
    return f"{OUTDIR}/qc/cutrun/{sample}/{sample}.fragment_lengths.tsv"


def preseq_path(sample):
    return f"{OUTDIR}/qc/cutrun/{sample}/{sample}.preseq.lc_extrap.txt"


def cutrun_qc_summary_path():
    return f"{OUTDIR}/qc/cutrun/cutrun_qc_summary.tsv"


def cutrun_qc_paths(sample):
    if not CUTRUN_QC_ENABLED:
        return []
    paths = [
        cutrun_usable_fragments_path(sample),
        fragment_length_path(sample),
    ]
    if PRESEQ_ENABLED:
        paths.append(preseq_path(sample))
    return paths


def spikein_bam_path(sample):
    return f"{OUTDIR}/spikein/bowtie2/{sample}/{sample}.spikein.bam"


def spikein_bai_path(sample):
    return spikein_bam_path(sample) + ".bai"


def spikein_flagstat_path(sample):
    return f"{OUTDIR}/qc/spikein/{sample}/{sample}.spikein.flagstat.txt"


def spikein_scale_factor_path(sample):
    return f"{OUTDIR}/qc/spikein/{sample}/{sample}.spikein_scale.tsv"


def spikein_scale_factors_path():
    return f"{OUTDIR}/qc/spikein/spikein_scale_factors.tsv"


def spikein_normalized_bigwig_path(sample):
    return f"{OUTDIR}/bigwig/{sample}/{sample}.spikein_normalized.bw"


def spikein_qc_paths(sample):
    if not SPIKEIN_ENABLED:
        return []
    return [
        spikein_flagstat_path(sample),
        spikein_scale_factor_path(sample),
    ]


def validate_workflow_config():
    if READ_TYPE not in {"PE", "SE"}:
        raise ValueError("read_type must be PE or SE.")
    if not SAMPLES:
        raise ValueError("config.samples must contain at least one sample.")

    for sample, sample_cfg in config.get("samples", {}).items():
        r1 = sample_cfg.get("R1", [])
        if not isinstance(r1, list) or not r1:
            raise ValueError(f"samples.{sample}.R1 must be a non-empty list.")
        if is_single_end():
            if sample_cfg.get("R2"):
                raise ValueError(f"samples.{sample}.R2 is set but read_type is SE.")
        else:
            r2 = sample_cfg.get("R2", [])
            if not isinstance(r2, list) or not r2:
                raise ValueError(f"samples.{sample}.R2 must be a non-empty list for PE data.")
            if len(r1) != len(r2):
                raise ValueError(
                    f"samples.{sample}.R1 and R2 must have the same number of FASTQs."
                )

    if not FILTER_BLACKLIST:
        pass
    elif not blacklist_path():
        raise ValueError(
            "filter_blacklist is enabled but no blacklist.path or blacklist.url is set."
        )

    if SPIKEIN_ENABLED:
        spikein = config.get("spikein", {})
        if not spikein.get("bowtie2_index"):
            raise ValueError("spikein.enabled is true but spikein.bowtie2_index is missing.")
        if float(spikein.get("scale_to", 10000)) <= 0:
            raise ValueError("spikein.scale_to must be greater than 0.")

    if MACS3_ENABLED:
        for key in ["gsize", "qvalue", "extsize"]:
            if key not in MACS3_CONFIG:
                raise ValueError(f"macs3.enabled is true but macs3.{key} is missing.")


rule faidx_reference:
    input:
        fa=config["reference"]["fasta"]
    output:
        fai=config["reference"]["fasta"] + ".fai"
    threads: 2
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools faidx {input.fa}"

validate_workflow_config()

if FILTER_BLACKLIST and blacklist_config().get("url"):
    rule get_blacklist:
        output:
            bed=blacklist_path()
        conda:
            "envs/samtools.yaml"
        script:
            "scripts/get_blacklist.py"
