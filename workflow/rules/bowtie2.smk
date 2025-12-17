# rules/bowtie2.smk

# This rule is for paired-end reads
rule bowtie2_pe:
    input:
        r1 = "results/fastp/{sample}_pe_fastp_R1.fastq.gz",
        r2 = "results/fastp/{sample}_pe_fastp_R2.fastq.gz"
    params:
        genome = config["bowtie2_index"]
    output:
        bam = "results/bowtie2/{sample}/{sample}_pe.bam",
    log:
        "logs/bowtie2/{sample}.log"
    threads: config["threads"]
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        bowtie2 \
            -x {params.genome} \
            -1 {input.r1} \
            -2 {input.r2} \
            --end-to-end \
            --very-sensitive \
            -X 500 \
            --no-mixed \
            --no-discordant \
            -p {threads} | \
            samtools view -bS - | \
            samtools sort -@ {threads} -o {output.bam} -
        echo "Running samtools index for {wildcards.sample}..." >> {log} 2>&1
        samtools index -@ {threads} {output.bam} >> {log} 2>&1
        """

# This rule is for single-end reads
rule biwtie2_se:
    input:
        r1 = "results/fastp/{sample}_se_fastp.fastq.gz"
    params:
        genome = config["bowtie2_index"]
    output:
        bam = "results/bowtie2/{sample}/{sample}_se.bam",
    log:
        "logs/bowtie2/{sample}.log"
    threads: config["threads"]
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2 \
            -x {params.genome} \
            -U {input.r1} \
            --end-to-end \
            --very-sensitive \
            -p {threads} | \
            samtools view -bS - | \
            samtools sort -@ {threads} -o {output.bam} -
        echo "Running samtools index for {wildcards.sample}..." >> {log} 2>&1
        samtools index -@ {threads} {output.bam} >> {log} 2>&1
        """

def get_bam_path(sample_name):
    if sample_name in PE_SAMPLES:
        return f"results/bowtie2/{sample_name}/{sample_name}_pe.bam"
    elif sample_name in SE_SAMPLES:
        return f"results/bowtie2/{sample_name}/{sample_name}_se.bam"
    else:
        raise ValueError(f"Sample {sample_name} not found in PE_SAMPLES or SE_SAMPLES.")
