# workflow/rules/common.smk
import os
from snakemake.io import temp

OUTDIR = config.get("output", {}).get("dir", "results")
SAMPLES = list(config.get("samples", {}).keys())
READ_TYPE = config.get("read_type", "PE").upper()

FILTER_BLACKLIST = bool(config.get("filter_blacklist", False))
REMOVE_DUPLICATES = bool(config.get("remove_duplicates", True))
KEEP_BAM = bool(config.get("keep_bam", False))


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

def validate_blacklist_config():
    if not FILTER_BLACKLIST:
        return
    if not blacklist_path():
        raise ValueError(
            "filter_blacklist is enabled but no blacklist.path or blacklist.url is set."
        )

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

validate_blacklist_config()

if FILTER_BLACKLIST and blacklist_config().get("url"):
    rule get_blacklist:
        output:
            bed=blacklist_path()
        conda:
            "envs/samtools.yaml"
        script:
            "scripts/get_blacklist.py"
