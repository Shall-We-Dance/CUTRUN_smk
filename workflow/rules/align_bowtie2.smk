# rules/bowtie2.smk

import os

# ----------------------------------------------------------------------------
# Alignment Rules
# ----------------------------------------------------------------------------

rule bowtie2_pe:
    input:
        r1=f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R2.fastq.gz"
    params:
        index=config["reference"]["bowtie2_index"],
        max_insert = config.get("bowtie2", {}).get("max_insert_size", 500),
        mode = config.get("bowtie2", {}).get("mode", "very-sensitive")
    output:
        bam = temp(f"{OUTDIR}/bowtie2/{{sample}}/{{sample}}_pe.bam"),
        bai = temp(f"{OUTDIR}/bowtie2/{{sample}}/{{sample}}_pe.bam.bai")
    log:
        "logs/bowtie2/{sample}_pe.log"
    threads: int(config["threads"]["bowtie2"])
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        bowtie2 \
            -x {params.index} \
            -1 {input.r1} \
            -2 {input.r2} \
            --end-to-end \
            --{params.mode} \
            -X {params.max_insert} \
            --no-mixed \
            --no-discordant \
            -p {threads} 2> {log} | \
        samtools view -bS - | \
        samtools sort -@ {threads} -o {output.bam} -
        
        echo "Running samtools index for {wildcards.sample}..." >> {log} 2>&1
        samtools index -@ {threads} {output.bam} >> {log} 2>&1
        # bowtie2 writes:
        # {OUTDIR}/tmp/bowtie2/
        """


# ----------------------------------------------------------------------------
# Filtering Rules
# ----------------------------------------------------------------------------

rule filter_unique_mappers:
    """
    Filter BAM to keep only uniquely mapped reads (MAPQ >= threshold).
    This removes multi-mapped reads.
    """
    input:
        bam = f"{OUTDIR}/bowtie2/{{sample}}/{{sample}}_pe.bam"
    output:
        bam = f"{OUTDIR}/bowtie2/{{sample}}/{{sample}}.unique.bam",
        bai = f"{OUTDIR}/bowtie2/{{sample}}/{{sample}}.unique.bam.bai"
    params:
        mapq = config.get("bowtie2", {}).get("mapq_threshold", 30),
        exclude_flags = config.get("samtools_exclude_flags", 1804)
    log:
        f"logs/bowtie2/{{sample}}.unique.log"
    threads: int(config["threads"]["samtools"])
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        # Filter for unique mappers: MAPQ >= {params.mapq}
        # Also exclude: unmapped, secondary, QC fail, duplicate
        samtools view -b -h -q {params.mapq} -F {params.exclude_flags} \
            -@ {threads} {input.bam} | \
        samtools sort -@ {threads} -o {output.bam} -
        
        samtools index -@ {threads} {output.bam}
        
        # Log filtering statistics
        echo "Original BAM:" > {log}
        samtools flagstat {input.bam} >> {log}
        echo -e "\nFiltered BAM (unique mappers):" >> {log}
        samtools flagstat {output.bam} >> {log}
        """

if FILTER_BLACKLIST:
    rule filter_blacklist_bam:
        """
        Remove reads overlapping blacklist regions.
        """
        input:
            bam = unique_bam_path("{sample}"),
            blacklist = blacklist_path()
        output:
            bam = filtered_bam_path("{sample}"),
            bai = filtered_bam_path("{sample}") + ".bai"
        params:
            nonamecheck = "--nonamecheck" if config.get("bedtools_nonamecheck", True) else ""
        threads: int(config["threads"]["samtools"])
        log:
            "logs/bowtie2/{sample}.filter_blacklist.log"
        conda:
            "envs/bowtie2.yaml"
        shell:
            """
            # Remove blacklist regions
            bedtools intersect -v \
                -abam {input.bam} \
                -b {input.blacklist} \
                {params.nonamecheck} | \
            samtools sort -@ {threads} -o {output.bam} -
            
            samtools index -@ {threads} {output.bam}
            
            # Log filtering statistics
            echo "Input BAM:" > {log}
            samtools flagstat {input.bam} >> {log}
            echo -e "\nAfter blacklist filtering:" >> {log}
            samtools flagstat {output.bam} >> {log}
            """


if REMOVE_DUPLICATES:
    rule remove_duplicates:
        """
        Remove PCR duplicates using Picard or samtools.
        """
        input:
            bam = lambda wc: filtered_bam_path(wc.sample) if FILTER_BLACKLIST else unique_bam_path(wc.sample)
        output:
            bam = dedup_bam_path("{sample}"),
            bai = dedup_bam_path("{sample}") + ".bai",
            metrics = f"{OUTDIR}/bowtie2/{{sample}}/{{sample}}.dedup_metrics.txt"
        params:
            method = config.get("dedup_method", "picard")
        threads: int(config["threads"]["samtools"])
        log:
            "logs/bowtie2/{sample}.dedup.log"
        conda:
            "envs/bowtie2.yaml"
        shell:
            """
            if [ "{params.method}" == "picard" ]; then
                picard MarkDuplicates \
                    I={input.bam} \
                    O={output.bam} \
                    M={output.metrics} \
                    REMOVE_DUPLICATES=true \
                    VALIDATION_STRINGENCY=LENIENT \
                    2> {log}
            else
                samtools markdup -r -@ {threads} {input.bam} {output.bam} 2> {log}
                samtools flagstat {output.bam} > {output.metrics}
            fi
            
            samtools index -@ {threads} {output.bam}
            """
