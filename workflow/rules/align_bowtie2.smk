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
            nonamecheck = config.get("bedtools_nonamecheck", True)
        threads: int(config["threads"]["samtools"])
        log:
            "logs/bowtie2/{sample}.filter_blacklist.log"
        conda:
            "envs/bowtie2.yaml"
        shell:
            """
            # Log filtering statistics
            echo "Input BAM:" > {log}
            samtools flagstat {input.bam} >> {log}

            nonamecheck_opt=""
            if [ "{params.nonamecheck}" = "True" ] && \
                bedtools intersect -h 2>&1 | grep -q -- '--nonamecheck'; then
                nonamecheck_opt="--nonamecheck"
            fi

            if [ ! -s "{input.blacklist}" ]; then
                echo -e "\nBlacklist BED is missing or empty; skipping blacklist filtering." >> {log}
                samtools view -b -@ {threads} {input.bam} -o {output.bam}
            else
                tmp_bam="$(mktemp --suffix=.bam)"
                # Remove blacklist regions
                if bedtools intersect -v \
                    -abam {input.bam} \
                    -b {input.blacklist} \
                    $nonamecheck_opt 2>> {log} | \
                samtools sort -@ {threads} -o "$tmp_bam" - 2>> {log}; then
                    mv "$tmp_bam" {output.bam}
                else
                    echo -e "\nBlacklist filtering failed; passing through input BAM." >> {log}
                    rm -f "$tmp_bam"
                    samtools view -b -@ {threads} {input.bam} -o {output.bam}
                fi
            fi

            samtools index -@ {threads} {output.bam} 2>> {log}

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
            method = config.get("dedup_method", "picard"),
            picard_java_opts = config.get("picard_java_opts", "-Xmx4g"),
            picard_tmpdir = config.get("picard_tmpdir", f"{OUTDIR}/tmp/picard")
        threads: int(config["threads"]["samtools"])
        log:
            "logs/bowtie2/{sample}.dedup.log"
        conda:
            "envs/bowtie2.yaml"
        shell:
            """
            mkdir -p $(dirname {output.bam}) "{params.picard_tmpdir}"
            if [ "{params.method}" == "picard" ]; then
                if ! command -v picard >/dev/null 2>&1; then
                    echo "picard not found in PATH; falling back to samtools markdup." > {log}
                    if samtools markdup -r -@ {threads} {input.bam} {output.bam} 2>> {log}; then
                        samtools flagstat {output.bam} > {output.metrics}
                    else
                        echo "samtools markdup failed; passing through input BAM." >> {log}
                        cp {input.bam} {output.bam}
                        samtools flagstat {output.bam} > {output.metrics}
                    fi
                else
                    if picard --java-options "{params.picard_java_opts}" MarkDuplicates \
                        I={input.bam} \
                        O={output.bam} \
                        M={output.metrics} \
                        TMP_DIR={params.picard_tmpdir} \
                        REMOVE_DUPLICATES=true \
                        VALIDATION_STRINGENCY=LENIENT \
                        2> {log}; then
                        :
                    else
                        echo "picard MarkDuplicates failed; falling back to samtools markdup." >> {log}
                        if samtools markdup -r -@ {threads} {input.bam} {output.bam} 2>> {log}; then
                            samtools flagstat {output.bam} > {output.metrics}
                        else
                            echo "samtools markdup failed; passing through input BAM." >> {log}
                            cp {input.bam} {output.bam}
                            samtools flagstat {output.bam} > {output.metrics}
                        fi
                    fi
                fi
            else
                if samtools markdup -r -@ {threads} {input.bam} {output.bam} 2> {log}; then
                    samtools flagstat {output.bam} > {output.metrics}
                else
                    echo "samtools markdup failed; passing through input BAM." >> {log}
                    cp {input.bam} {output.bam}
                    samtools flagstat {output.bam} > {output.metrics}
                fi
            fi

            samtools index -@ {threads} {output.bam}
            """
