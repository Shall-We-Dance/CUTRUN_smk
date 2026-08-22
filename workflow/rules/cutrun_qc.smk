# workflow/rules/cutrun_qc.smk

if CUTRUN_QC_ENABLED:
    rule cutrun_usable_fragments:
        input:
            bam=lambda wc: final_bam_path(wc.run_sample),
            bai=lambda wc: final_bai_path(wc.run_sample)
        output:
            tsv=cutrun_usable_fragments_path("{run_sample}")
        params:
            exclude_flags=config.get("samtools_exclude_flags", 1804),
            read_type=READ_TYPE
        log:
            f"{LOG_DIR}/cutrun_qc/{{run_sample}}.usable_fragments.log"
        threads: int(config["threads"]["samtools"])
        resources:
            mem_mb=4000,
            disk_mb=20000
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
                printf "sample\t{wildcards.run_sample}\n"
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
            bam=lambda wc: final_bam_path(wc.run_sample),
            bai=lambda wc: final_bai_path(wc.run_sample)
        output:
            tsv=fragment_length_path("{run_sample}")
        params:
            exclude_flags=config.get("samtools_exclude_flags", 1804),
            read_type=READ_TYPE
        log:
            f"{LOG_DIR}/cutrun_qc/{{run_sample}}.fragment_lengths.log"
        threads: int(config["threads"]["samtools"])
        resources:
            mem_mb=4000,
            disk_mb=20000
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

    rule summarize_cutrun_qc:
        input:
            usable=[cutrun_usable_fragments_path(run_sample) for run_sample in RUN_SAMPLES],
            fragments=[fragment_length_path(run_sample) for run_sample in RUN_SAMPLES],
            spikein=[spikein_scale_factors_path()] if SPIKEIN_ENABLED else []
        output:
            tsv=cutrun_qc_summary_path(),
            fragments_mqc=cutrun_fragment_lengths_multiqc_path(),
            summary_mqc=cutrun_summary_multiqc_path()
        conda:
            "envs/cutrun_qc.yaml"
        script:
            "scripts/summarize_cutrun_qc.py"

if MACS3_ENABLED:
    rule frip:
        """Fragment-level FRiP for a replicate's final BAM."""
        input:
            bam=lambda wc: final_bam_path(wc.run_sample),
            peaks=lambda wc: macs3_narrow_path(wc.run_sample),
        output:
            tsv=frip_path("{run_sample}")
        params:
            read_type=READ_TYPE
        log:
            f"{LOG_DIR}/cutrun_qc/{{run_sample}}.frip.log"
        threads: int(config["threads"]["samtools"])
        resources:
            mem_mb=4000,
            disk_mb=20000
        conda:
            "envs/bowtie2.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log})
            python workflow/rules/scripts/frip.py \
                --bam {input.bam} \
                --peaks {input.peaks} \
                --read-type {params.read_type} \
                --sample {wildcards.run_sample} \
                --threads {threads} \
                --out {output.tsv} > {log} 2>&1
            """

    rule frip_pooled:
        """Fragment-level FRiP for pooled BAM against pooled peaks."""
        wildcard_constraints:
            sample=POOLED_PEAK_SAMPLE_PATTERN
        input:
            bam=lambda wc: pooled_final_bam_path(wc.sample),
            peaks=lambda wc: macs3_narrow_path(wc.sample, pooled=True),
        output:
            tsv=frip_path("{sample}", pooled=True)
        params:
            read_type=READ_TYPE
        log:
            f"{LOG_DIR}/cutrun_qc/{{sample}}.pooled.frip.log"
        threads: int(config["threads"]["samtools"])
        resources:
            mem_mb=4000,
            disk_mb=20000
        conda:
            "envs/bowtie2.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log})
            python workflow/rules/scripts/frip.py \
                --bam {input.bam} \
                --peaks {input.peaks} \
                --read-type {params.read_type} \
                --sample {wildcards.sample} \
                --threads {threads} \
                --out {output.tsv} > {log} 2>&1
            """

    rule summarize_frip:
        input:
            per_rep=[frip_path(rs) for rs in RUN_SAMPLES if rs in PEAK_RUN_SAMPLES],
            pooled=[frip_path(sample, pooled=True) for sample in POOLED_PEAK_SAMPLES],
        output:
            tsv=frip_summary_path(),
            mqc=frip_multiqc_path()
        script:
            "scripts/summarize_frip.py"
