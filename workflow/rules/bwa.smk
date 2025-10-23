# rules/bwa.smk

# This rule is for paired-end reads
rule bwa_pe:
    input:
        r1 = "results/fastp/{sample}_pe_fastp_R1.fastq.gz",
        r2 = "results/fastp/{sample}_pe_fastp_R2.fastq.gz"
    params:
        genome = config["bwa_index"]
    output:
        bam = "results/bwa/{sample}/{sample}_pe.bam",
    log:
        "logs/bwa/{sample}.log"
    threads: config["threads"]
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {params.genome} {input.r1} {input.r2} | \
        samtools view -bS - | \
        samtools sort -o {output.bam} -
        """

# This rule is for single-end reads
rule bwa_se:
    input:
        r1 = "results/fastp/{sample}_se_fastp.fastq.gz"
    params:
        genome = config["bwa_index"]
    output:
        bam = "results/bwa/{sample}/{sample}_se.bam",
    log:
        "logs/bwa/{sample}.log"
    threads: config["threads"]
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {params.genome} {input.r1} | \
        samtools view -bS - | \
        samtools sort -o {output.bam} -
        """
