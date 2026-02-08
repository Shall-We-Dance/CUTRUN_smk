# workflow/rules/common.smk
import os

OUTDIR = config["output"]["dir"]

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

rule get_blacklist:
    output:
        bed=f"{config['blacklist']['cache_dir']}/{config['genome']}-blacklist.bed"
    conda:
        "envs/py_signal.yaml"
    script:
        "scripts/get_blacklist.py"
