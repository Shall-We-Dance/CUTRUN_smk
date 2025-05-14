# rules/bam_to_bigwig.smk

import os
import hashlib
import urllib.request

def get_bam_path(sample):
    return f"results/star/{sample}/{sample}_Aligned.sortedByCoord.out.bam"

def get_blacklist_path(genome):
    """
    Download and cache ENCODE blacklist BED file for a given genome.
    """
    url_dict = {
        "mm10": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz",
        "hg19": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz",
        "hg38": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz",
        "ce11": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce11-blacklist.v2.bed.gz",
        "ce10": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce10-blacklist.v2.bed.gz",
        "dm3": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm3-blacklist.v2.bed.gz",
        "dm6": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.v2.bed.gz",
    }
    if genome not in url_dict:
        raise ValueError(f"Unsupported genome: {genome}")

    cache_dir = "resources/blacklist"
    os.makedirs(cache_dir, exist_ok=True)

    url = url_dict[genome]
    fname = os.path.join(cache_dir, os.path.basename(url))

    if not os.path.exists(fname):
        print(f"Downloading blacklist for {genome}...")
        urllib.request.urlretrieve(url, fname)
    return fname

rule bam_to_bigwig:
    input:
        bam = lambda wildcards: get_bam_path(wildcards.sample)
    params:
        chrom_sizes = config["chrom_sizes"],
        bin_size = config["bw_bin_size"]
    output:
        raw_bw = "results/bigwig/{sample}.bw",
        rpkm_bw = "results/bigwig/{sample}.rpkm.bw"
    log:
        "logs/bigwig/{sample}.log"
    threads: config["threads"]
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        # Unnormalized
        bamCoverage --bam {input.bam} \
                    --outFileName {output.raw_bw} \
                    --outFileFormat bigwig \
                    --numberOfProcessors {threads} \
                    --binSize {params.bin_size} \
                    --normalizeUsing None \
                    --samFlagExclude 1024 \
                    --effectiveGenomeSize $(awk '{{sum += $2}} END {{print sum}}' {params.chrom_sizes}) \
                    > {log} 2>&1
        # RPKM normalized
        bamCoverage --bam {input.bam} \
                    --outFileName {output.rpkm_bw} \
                    --outFileFormat bigwig \
                    --numberOfProcessors {threads} \
                    --binSize {params.bin_size} \
                    --normalizeUsing RPKM \
                    --samFlagExclude 1024 \
                    --effectiveGenomeSize $(awk '{{sum += $2}} END {{print sum}}' {params.chrom_sizes}) \
                    >> {log} 2>&1
        """

rule bam_to_bigwig_blacklist:
    input:
        bam = lambda wildcards: get_bam_path(wildcards.sample)
    params:
        chrom_sizes = config["chrom_sizes"],
        bin_size = config["bw_bin_size"],
        blacklist = lambda wildcards: get_blacklist_path(config["genome"])
    output:
        raw_bw = "results/bigwig/{sample}_blacklist.bw",
        rpkm_bw = "results/bigwig/{sample}_blacklist.rpkm.bw"
    log:
        "logs/bigwig/{sample}_blacklist.log"
    threads: config["threads"]
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        bamCoverage --bam {input.bam} \
                    --outFileName {output.raw_bw} \
                    --outFileFormat bigwig \
                    --numberOfProcessors {threads} \
                    --binSize {params.bin_size} \
                    --normalizeUsing None \
                    --samFlagExclude 1024 \
                    --blackListFileName {params.blacklist} \
                    --effectiveGenomeSize $(awk '{{sum += $2}} END {{print sum}}' {params.chrom_sizes}) \
                    > {log} 2>&1

        bamCoverage --bam {input.bam} \
                    --outFileName {output.rpkm_bw} \
                    --outFileFormat bigwig \
                    --numberOfProcessors {threads} \
                    --binSize {params.bin_size} \
                    --normalizeUsing RPKM \
                    --samFlagExclude 1024 \
                    --blackListFileName {params.blacklist} \
                    --effectiveGenomeSize $(awk '{{sum += $2}} END {{print sum}}' {params.chrom_sizes}) \
                    >> {log} 2>&1
        """