# rules/bam_to_bigwig.smk

import os
import hashlib
import urllib.request
import gzip
import shutil

rule bam_to_bigwig:
    input:
        bam = lambda wildcards: get_raw_bam(wildcards.sample)
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
        r"""
        # Unnormalized
        bamCoverage \
            --bam {input.bam} \
            --outFileName {output.raw_bw} \
            --outFileFormat bigwig \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --normalizeUsing None \
            --samFlagExclude 1024 \
            --effectiveGenomeSize $(awk '{{sum += $2}} END {{print sum}}' {params.chrom_sizes}) \
            > {log} 2>&1

        # RPKM normalized
        bamCoverage \
            --bam {input.bam} \
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
        bam = lambda wildcards: get_filtered_bam(wildcards.sample)
    params:
        chrom_sizes = config["chrom_sizes"],
        bin_size = config["bw_bin_size"]
    output:
        raw_bw = "results/bigwig/{sample}.filtered.bw",
        rpkm_bw = "results/bigwig/{sample}.filtered.rpkm.bw"
    log:
        "logs/bigwig/{sample}.filtered.log"
    threads: config["threads"]
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        # Unnormalized
        bamCoverage \
            --bam {input.bam} \
            --outFileName {output.raw_bw} \
            --outFileFormat bigwig \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --normalizeUsing None \
            --samFlagExclude 1024 \
            --effectiveGenomeSize $(awk '{{sum += $2}} END {{print sum}}' {params.chrom_sizes}) \
            > {log} 2>&1

        # RPKM normalized
        bamCoverage \
            --bam {input.bam} \
            --outFileName {output.rpkm_bw} \
            --outFileFormat bigwig \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --normalizeUsing RPKM \
            --samFlagExclude 1024 \
            --effectiveGenomeSize $(awk '{{sum += $2}} END {{print sum}}' {params.chrom_sizes}) \
            >> {log} 2>&1
        """
