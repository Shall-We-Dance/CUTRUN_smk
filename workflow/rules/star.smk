# rules/star.smk

ruleorder: star_pe > star_se

# This rule is for paired-end reads
rule star_pe:
    input:
        r1 = "results/fastp/{sample}_fastp_R1.fastq.gz",
        r2 = "results/fastp/{sample}_fastp_R2.fastq.gz"
    params:
        genome = config["STAR_index"]
    output:
        bam = "results/star/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        metrics = "results/star/{sample}/{sample}_metrics.txt"
    log:
        "logs/star/{sample}.log"
    threads: config["threads"]
    conda:
        "../envs/star.yaml"
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {params.genome} \
             --readFilesCommand zcat \
             --readFilesIn {input.r1} {input.r2} \
             --outFileNamePrefix results/star/{wildcards.sample}/{wildcards.sample}_ \
             --outFilterMultimapNmax 100 \
             --outSAMmultNmax 1 \
             --outSAMattributes NH HI AS nM \
             --winAnchorMultimapNmax 100 \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped None > {log} 2>&1

        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.metrics}
        """

# This rule is for single-end reads
rule star_se:
    input:
        r1 = "results/fastp/{sample}_fastp.fastq.gz"
    params:
        genome = config["STAR_index"]
    output:
        bam = "results/star/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        metrics = "results/star/{sample}/{sample}_metrics.txt"
    log:
        "logs/star/{sample}.log"
    threads: config["threads"]
    conda:
        "../envs/star.yaml"
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {params.genome} \
             --readFilesIn {input.r1} \
             --outFileNamePrefix results/star/{wildcards.sample}/{wildcards.sample}_ \
             --outFilterMultimapNmax 100 \
             --outSAMmultNmax 1 \
             --outSAMattributes NH HI AS nM \
             --winAnchorMultimapNmax 100 \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped None > {log} 2>&1

        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.metrics}
        """

# This rule is for masking repeats in the BAM file
if config.get("mask_repeats") and config.get("repeat_bed"):
    rule mask_repeats:
        input:
            bam = "results/star/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
        params:
            repeats = config["repeat_bed"]
        output:
            masked_bam = "results/star/{sample}/{sample}_masked.bam"
        log:
            "logs/star/{sample}_mask.log"
        threads: 1
        conda:
            "envs/bedtools.yaml"
        shell:
            """
            bedtools intersect -v -abam {input.bam} -b {params.repeats} > {output.masked_bam} 2> {log}
            samtools index {output.masked_bam}
            """
