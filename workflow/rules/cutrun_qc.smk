# workflow/rules/cutrun_qc.smk

if CUTRUN_QC_ENABLED:
    rule cutrun_usable_fragments:
        input:
            bam=lambda wc: final_bam_path(wc.sample),
            bai=lambda wc: final_bai_path(wc.sample)
        output:
            tsv=cutrun_usable_fragments_path("{sample}")
        params:
            exclude_flags=config.get("samtools_exclude_flags", 1804),
            read_type=READ_TYPE
        log:
            "logs/cutrun_qc/{sample}.usable_fragments.log"
        threads: int(config["threads"]["samtools"])
        conda:
            "envs/cutrun_qc.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log})

            total_reads=$(samtools view -c {input.bam} 2>> {log})
            mapped_reads=$(samtools view -c -F 4 {input.bam} 2>> {log})
            usable_reads=$(samtools view -c -F {params.exclude_flags} {input.bam} 2>> {log})

            if [ "{params.read_type}" = "PE" ]; then
                proper_pair_reads=$(samtools view -c -f 2 -F {params.exclude_flags} {input.bam} 2>> {log})
                usable_fragments=$((proper_pair_reads / 2))
            else
                proper_pair_reads="NA"
                usable_fragments="$usable_reads"
            fi

            {{
                printf "metric\tvalue\n"
                printf "sample\t{wildcards.sample}\n"
                printf "read_type\t{params.read_type}\n"
                printf "total_reads_after_processing\t%s\n" "$total_reads"
                printf "mapped_reads_after_processing\t%s\n" "$mapped_reads"
                printf "usable_reads_after_processing\t%s\n" "$usable_reads"
                printf "proper_pair_reads_after_processing\t%s\n" "$proper_pair_reads"
                printf "usable_fragments_after_processing\t%s\n" "$usable_fragments"
            }} > {output.tsv}
            """

    rule cutrun_fragment_lengths:
        input:
            bam=lambda wc: final_bam_path(wc.sample),
            bai=lambda wc: final_bai_path(wc.sample)
        output:
            tsv=fragment_length_path("{sample}")
        params:
            exclude_flags=config.get("samtools_exclude_flags", 1804),
            read_type=READ_TYPE
        log:
            "logs/cutrun_qc/{sample}.fragment_lengths.log"
        threads: int(config["threads"]["samtools"])
        conda:
            "envs/cutrun_qc.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log})

            if [ "{params.read_type}" = "PE" ]; then
                samtools view -f 2 -F {params.exclude_flags} {input.bam} 2>> {log} | \
                awk 'BEGIN {{OFS="\t"; print "fragment_length","count"}} $9 > 0 {{counts[$9]++}} END {{for (flen in counts) print flen, counts[flen]}}' | \
                sort -k1,1n > {output.tsv}
            else
                samtools view -F {params.exclude_flags} {input.bam} 2>> {log} | \
                awk 'BEGIN {{OFS="\t"; print "fragment_length","count"}} {{counts[length($10)]++}} END {{for (flen in counts) print flen, counts[flen]}}' | \
                sort -k1,1n > {output.tsv}
            fi
            """

    if PRESEQ_ENABLED:
        if not is_single_end():
            rule preseq_fastp_sample_level:
                input:
                    r1=f"{OUTDIR}/tmp/merged_raw/{{sample}}_R1.fastq.gz",
                    r2=f"{OUTDIR}/tmp/merged_raw/{{sample}}_R2.fastq.gz"
                output:
                    clean_r1=temp(preseq_clean_r1_path("{sample}")),
                    clean_r2=temp(preseq_clean_r2_path("{sample}")),
                    html=temp(f"{OUTDIR}/qc/preseq/{{sample}}/{{sample}}.preseq_fastp.html"),
                    json=temp(f"{OUTDIR}/qc/preseq/{{sample}}/{{sample}}.preseq_fastp.json")
                log:
                    "logs/preseq/{sample}.fastp.log"
                threads: int(config["threads"]["fastp"])
                conda:
                    "envs/qc.yaml"
                params:
                    extra_args=fastp_extra_args(force_dedup=False)
                shell:
                    r"""
                    set -euo pipefail
                    mkdir -p $(dirname {output.clean_r1}) $(dirname {output.html}) $(dirname {log})

                    fastp \
                      -i {input.r1} -I {input.r2} \
                      -o {output.clean_r1} -O {output.clean_r2} \
                      --thread {threads} \
                      {params.extra_args} \
                      --html {output.html} --json {output.json} \
                      > {log} 2>&1
                    """

            rule preseq_bowtie2:
                input:
                    r1=preseq_clean_r1_path("{sample}"),
                    r2=preseq_clean_r2_path("{sample}")
                output:
                    bam=temp(preseq_aligned_bam_path("{sample}")),
                    bai=temp(preseq_aligned_bam_path("{sample}") + ".bai"),
                    flagstat=preseq_aligned_flagstat_path("{sample}")
                params:
                    index=config["reference"]["bowtie2_index"],
                    max_insert=config.get("bowtie2", {}).get("max_insert_size", 500),
                    mode=config.get("bowtie2", {}).get("mode", "very-sensitive")
                log:
                    "logs/preseq/{sample}.bowtie2_pe.log"
                threads: int(config["threads"]["bowtie2"])
                conda:
                    "envs/bowtie2.yaml"
                shell:
                    r"""
                    set -euo pipefail
                    mkdir -p $(dirname {output.bam}) $(dirname {output.flagstat}) $(dirname {log})

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

                    samtools index -@ {threads} {output.bam} >> {log} 2>&1
                    samtools flagstat {output.bam} > {output.flagstat}
                    """
        else:
            rule preseq_fastp_sample_level:
                input:
                    r1=f"{OUTDIR}/tmp/merged_raw/{{sample}}_R1.fastq.gz"
                output:
                    clean_r1=temp(preseq_clean_r1_path("{sample}")),
                    html=temp(f"{OUTDIR}/qc/preseq/{{sample}}/{{sample}}.preseq_fastp.html"),
                    json=temp(f"{OUTDIR}/qc/preseq/{{sample}}/{{sample}}.preseq_fastp.json")
                log:
                    "logs/preseq/{sample}.fastp.log"
                threads: int(config["threads"]["fastp"])
                conda:
                    "envs/qc.yaml"
                params:
                    extra_args=fastp_extra_args(force_dedup=False)
                shell:
                    r"""
                    set -euo pipefail
                    mkdir -p $(dirname {output.clean_r1}) $(dirname {output.html}) $(dirname {log})

                    fastp \
                      -i {input.r1} \
                      -o {output.clean_r1} \
                      --thread {threads} \
                      {params.extra_args} \
                      --html {output.html} --json {output.json} \
                      > {log} 2>&1
                    """

            rule preseq_bowtie2:
                input:
                    r1=preseq_clean_r1_path("{sample}")
                output:
                    bam=temp(preseq_aligned_bam_path("{sample}")),
                    bai=temp(preseq_aligned_bam_path("{sample}") + ".bai"),
                    flagstat=preseq_aligned_flagstat_path("{sample}")
                params:
                    index=config["reference"]["bowtie2_index"],
                    mode=config.get("bowtie2", {}).get("mode", "very-sensitive")
                log:
                    "logs/preseq/{sample}.bowtie2_se.log"
                threads: int(config["threads"]["bowtie2"])
                conda:
                    "envs/bowtie2.yaml"
                shell:
                    r"""
                    set -euo pipefail
                    mkdir -p $(dirname {output.bam}) $(dirname {output.flagstat}) $(dirname {log})

                    bowtie2 \
                        -x {params.index} \
                        -U {input.r1} \
                        --end-to-end \
                        --{params.mode} \
                        -p {threads} 2> {log} | \
                    samtools view -bS - | \
                    samtools sort -@ {threads} -o {output.bam} -

                    samtools index -@ {threads} {output.bam} >> {log} 2>&1
                    samtools flagstat {output.bam} > {output.flagstat}
                    """

        rule preseq_filter_unique_mappers:
            input:
                bam=preseq_aligned_bam_path("{sample}"),
                bai=preseq_aligned_bam_path("{sample}") + ".bai"
            output:
                bam=temp(preseq_unique_bam_path("{sample}")),
                bai=temp(preseq_unique_bam_path("{sample}") + ".bai"),
                flagstat=preseq_unique_flagstat_path("{sample}")
            params:
                mapq=config.get("bowtie2", {}).get("mapq_threshold", 30),
                exclude_flags=config.get("samtools_exclude_flags", 1804)
            log:
                "logs/preseq/{sample}.unique.log"
            threads: int(config["threads"]["samtools"])
            conda:
                "envs/bowtie2.yaml"
            shell:
                r"""
                set -euo pipefail
                mkdir -p $(dirname {output.bam}) $(dirname {output.flagstat}) $(dirname {log})

                samtools view -b -h -q {params.mapq} -F {params.exclude_flags} \
                    -@ {threads} {input.bam} | \
                samtools sort -@ {threads} -o {output.bam} -

                samtools index -@ {threads} {output.bam}
                samtools flagstat {output.bam} > {output.flagstat}
                """

        if FILTER_BLACKLIST:
            rule preseq_filter_blacklist_bam:
                input:
                    bam=preseq_unique_bam_path("{sample}"),
                    bai=preseq_unique_bam_path("{sample}") + ".bai",
                    blacklist=blacklist_path()
                output:
                    bam=temp(preseq_filtered_bam_path("{sample}")),
                    bai=temp(preseq_filtered_bam_path("{sample}") + ".bai"),
                    flagstat=preseq_filtered_flagstat_path("{sample}")
                params:
                    nonamecheck=config.get("bedtools_nonamecheck", True)
                threads: int(config["threads"]["samtools"])
                log:
                    "logs/preseq/{sample}.filter_blacklist.log"
                conda:
                    "envs/bowtie2.yaml"
                shell:
                    r"""
                    set -euo pipefail
                    mkdir -p $(dirname {output.bam}) $(dirname {output.flagstat}) $(dirname {log})

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
                    samtools flagstat {output.bam} > {output.flagstat}
                    """

        rule preseq_library_complexity:
            input:
                bam=lambda wc: preseq_input_bam_path(wc.sample),
                bai=lambda wc: preseq_input_bai_path(wc.sample)
            output:
                txt=preseq_path("{sample}")
            params:
                pe_arg="" if is_single_end() else "-pe"
            log:
                "logs/preseq/{sample}.preseq.log"
            threads: 1
            conda:
                "envs/cutrun_qc.yaml"
            shell:
                r"""
                set -euo pipefail
                mkdir -p $(dirname {output.txt}) $(dirname {log})
                preseq lc_extrap {params.pe_arg} -output {output.txt} {input.bam} > {log} 2>&1
                """

    rule summarize_cutrun_qc:
        input:
            usable=[cutrun_usable_fragments_path(sample) for sample in SAMPLES],
            fragments=[fragment_length_path(sample) for sample in SAMPLES],
            preseq=[preseq_path(sample) for sample in SAMPLES] if PRESEQ_ENABLED else [],
            spikein=[spikein_scale_factors_path()] if SPIKEIN_ENABLED else []
        output:
            tsv=cutrun_qc_summary_path()
        conda:
            "envs/cutrun_qc.yaml"
        script:
            "scripts/summarize_cutrun_qc.py"
