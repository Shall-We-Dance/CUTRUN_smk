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
             --outSAMmultNmax 100 \
             --outSAMattributes NH HI AS nM \
             --winAnchorMultimapNmax 100 \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within > {log} 2>&1

        samtools index -@ {threads} {output.bam}
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
             --readFilesCommand zcat \
             --readFilesIn {input.r1} \
             --outFileNamePrefix results/star/{wildcards.sample}/{wildcards.sample}_ \
             --outFilterMultimapNmax 100 \
             --outSAMmultNmax 100 \
             --outSAMattributes NH HI AS nM \
             --winAnchorMultimapNmax 100 \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within > {log} 2>&1

        samtools index -@ {threads} {output.bam}
        samtools flagstat {output.bam} > {output.metrics}
        """

rule multimap_weight:
    input:
        bam = "results/star/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        bed = "results/star/{sample}/{sample}_multimap_weight.bed"
    log:
        "logs/star/{sample}_weight.log"
    threads: config["threads"]
    conda:
        "../envs/star.yaml"
    shell:
        """
        samtools view -@ {threads} {input.bam} | \
        awk '{{
            for(i=12;i<=NF;i++) {{
                if($i ~ /^NH:i:/) {{
                    split($i,a,":")
                    nh=a[3]
                    break
                }}
            }}
            if(nh) {{
                weight=1/nh
                print $3"\t"$4-1"\t"$4+length($10)"\t.\t"weight"\t"$2
            }}
        }}' > {output.bed} 2> {log}
        """


rule macs3_callpeak:
    input:
        bam = "results/star/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        narrowpeak = "results/macs3/{sample}/{sample}_peaks.narrowPeak"
    params:
        outdir = "results/macs3/{sample}",
        name = "{sample}",
        tempdir = "results/macs3/{sample}/temp",
        gsize = config["macs_gsize"]
    log:
        "logs/macs3/{sample}_macs3.log"
    conda:
        "envs/macs3.yaml"
    threads: 1
    shell:
        """
        mkdir -p {params.tempdir}
        macs3 callpeak -f AUTO \
                       -t {input.bam} \
                       -g {params.gsize} \
                       --outdir {params.outdir} \
                       -n {params.name} \
                       --tempdir {params.tempdir} \
                       --call-summits > {log} 2>&1
        """