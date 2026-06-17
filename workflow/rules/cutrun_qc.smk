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

    rule summarize_cutrun_qc:
        input:
            usable=[cutrun_usable_fragments_path(sample) for sample in SAMPLES],
            fragments=[fragment_length_path(sample) for sample in SAMPLES],
            spikein=[spikein_scale_factors_path()] if SPIKEIN_ENABLED else []
        output:
            tsv=cutrun_qc_summary_path(),
            mqc=cutrun_fragment_lengths_multiqc_path()
        conda:
            "envs/cutrun_qc.yaml"
        script:
            "scripts/summarize_cutrun_qc.py"
