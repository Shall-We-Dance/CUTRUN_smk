# rules/star.smk

ruleorder: star_pe > star_se

# This rule is for paired-end reads
rule star_pe:
    input:
        r1 = "results/fastp/{sample}_pe_fastp_R1.fastq.gz",
        r2 = "results/fastp/{sample}_pe_fastp_R2.fastq.gz"
    params:
        genome = config["STAR_index"]
    output:
        bam = "results/star/{sample}/{sample}_pe_Aligned.sortedByCoord.out.bam",
        link = "results/star/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
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
        r1 = "results/fastp/{sample}_se_fastp.fastq.gz"
    params:
        genome = config["STAR_index"]
    output:
        bam = "results/star/{sample}/{sample}_se_Aligned.sortedByCoord.out.bam",
        link = "results/star/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
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

def get_bam_path(sample_name):
    if sample_name in PE_SAMPLES:
        return f"results/star/{sample_name}/{sample_name}_pe_Aligned.sortedByCoord.out.bam"
    elif sample_name in SE_SAMPLES:
        return f"results/star/{sample_name}/{sample_name}_se_Aligned.sortedByCoord.out.bam"
    else:
        raise ValueError(f"Sample {sample_name} not found in PE_SAMPLES or SE_SAMPLES.")

rule multimap_weight:
    input:
        bam = lambda wildcards: get_bam_path(wildcards.sample)
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
