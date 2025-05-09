# rules/fastp.smk

import os
import glob
def find_pe_read_file(sample, read):
    patterns = [
        f"{sample}_R{read}.fastq.gz",
        f"{sample}_{read}.fastq.gz",
        f"{sample}_R{read}.fq.gz",
        f"{sample}_{read}.fq.gz"
    ]
    recursive = config.get("auto_detect_subdirs", False)
    for pat in patterns:
        # Use recursive search with '**' if recursive is True
        path = os.path.join(config["raw_dir"], "**" if recursive else "", pat)
        matches = glob.glob(path, recursive=recursive)
        if matches:
            return matches[0]
    raise FileNotFoundError(f"No PE file found for {sample}, read {read}")

def find_se_read_file(sample):
    patterns = [
        f"{sample}.fastq.gz", f"{sample}.fq.gz",
        f"{sample}_R1.fastq.gz", f"{sample}_1.fastq.gz",
        f"{sample}_R1.fq.gz", f"{sample}_1.fq.gz"
    ]
    recursive = config.get("auto_detect_subdirs", False)
    for pat in patterns:
        # Use recursive search with '**' if recursive is True
        path = os.path.join(config["raw_dir"], "**" if recursive else "", pat)
        matches = glob.glob(path, recursive=recursive)
        if matches:
            return matches[0]
    raise FileNotFoundError(f"No SE file found for {sample}")

rule fastp_pe:
    input:
        r1 = lambda wildcards: find_pe_read_file(wildcards.sample, "1"),
        r2 = lambda wildcards: find_pe_read_file(wildcards.sample, "2")
    output:
        r1_out = temp("results/fastp/{sample}_fastp_R1.fastq.gz"),
        r2_out = temp("results/fastp/{sample}_fastp_R2.fastq.gz"),
        html = "results/fastp/{sample}_pe_fastp.html",
        json = "results/fastp/{sample}_pe_fastp.json"
    log:
        "logs/fastp/{sample}.log"
    threads: config["threads"]
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1_out} -O {output.r2_out} \
              -h {output.html} -j {output.json} \
              --thread {threads} > {log} 2>&1
        """

rule fastp_se:
    input:
        r1 = lambda wildcards: find_se_read_file(wildcards.sample)
    output:
        r1_out = temp("results/fastp/{sample}_fastp.fastq.gz"),
        html = "results/fastp/{sample}_se_fastp.html",
        json = "results/fastp/{sample}_se_fastp.json"
    log:
        "logs/fastp/{sample}.log"
    threads: config["threads"]
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        fastp -i {input.r1} \
              -o {output.r1_out} \
              -h {output.html} -j {output.json} \
              --thread {threads} > {log} 2>&1
        """