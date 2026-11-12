# workflow/rules/spikein.smk

if SPIKEIN_ENABLED:
    if not is_single_end():
        rule bowtie2_spikein_pe:
            input:
                r1=fastp_sample_r1("{run_sample}"),
                r2=fastp_sample_r2("{run_sample}")
            output:
                bam=spikein_bam_path("{run_sample}"),
                bai=spikein_bai_path("{run_sample}"),
                flagstat=spikein_flagstat_path("{run_sample}")
            params:
                index=config["spikein"]["bowtie2_index"],
                max_insert=config.get("bowtie2", {}).get("max_insert_size", 500),
                mode=config.get("bowtie2", {}).get("mode", "very-sensitive"),
                sort_mem=config.get("samtools_sort_memory", "2G"),
                pipe_split=pipe_thread_split()
            log:
                f"{LOG_DIR}/spikein/{{run_sample}}_pe.log"
            threads: int(config["threads"]["bowtie2"])
            conda:
                "envs/bowtie2.yaml"
            shell:
                r"""
                set -euo pipefail
                mkdir -p $(dirname {output.bam}) $(dirname {output.flagstat}) $(dirname {log})
                THREADS={threads}
            {params.pipe_split}
                echo "granted_threads={threads} bowtie2_p=$MAIN_T sort_at=$SORT_T serial=$SERIAL" > {log}
                tmp_sam="$(mktemp --suffix=.sam)"
                trap 'rm -f "$tmp_sam"' EXIT
                if [ "$SERIAL" = "1" ]; then
                    bowtie2 \
                        -x {params.index} \
                        -1 {input.r1} \
                        -2 {input.r2} \
                        --end-to-end \
                        --{params.mode} \
                        -X {params.max_insert} \
                        --no-mixed \
                        --no-discordant \
                        -p 1 -S "$tmp_sam" 2>> {log}
                    samtools sort -@ 0 -m {params.sort_mem} -o {output.bam} "$tmp_sam" 2>> {log}
                else
                    bowtie2 \
                        -x {params.index} \
                        -1 {input.r1} \
                        -2 {input.r2} \
                        --end-to-end \
                        --{params.mode} \
                        -X {params.max_insert} \
                        --no-mixed \
                        --no-discordant \
                        -p $MAIN_T 2>> {log} | \
                    samtools view -bS - | \
                    samtools sort -@ $SORT_T -m {params.sort_mem} -o {output.bam} -
                fi
                samtools index -@ $SORT_T {output.bam} >> {log} 2>&1
                samtools flagstat {output.bam} > {output.flagstat}
                """
    else:
        rule bowtie2_spikein_se:
            input:
                r1=fastp_sample_r1("{run_sample}")
            output:
                bam=spikein_bam_path("{run_sample}"),
                bai=spikein_bai_path("{run_sample}"),
                flagstat=spikein_flagstat_path("{run_sample}")
            params:
                index=config["spikein"]["bowtie2_index"],
                mode=config.get("bowtie2", {}).get("mode", "very-sensitive"),
                sort_mem=config.get("samtools_sort_memory", "2G"),
                pipe_split=pipe_thread_split()
            log:
                f"{LOG_DIR}/spikein/{{run_sample}}_se.log"
            threads: int(config["threads"]["bowtie2"])
            conda:
                "envs/bowtie2.yaml"
            shell:
                r"""
                set -euo pipefail
                mkdir -p $(dirname {output.bam}) $(dirname {output.flagstat}) $(dirname {log})
                THREADS={threads}
            {params.pipe_split}
                echo "granted_threads={threads} bowtie2_p=$MAIN_T sort_at=$SORT_T serial=$SERIAL" > {log}
                tmp_sam="$(mktemp --suffix=.sam)"
                trap 'rm -f "$tmp_sam"' EXIT
                if [ "$SERIAL" = "1" ]; then
                    bowtie2 \
                        -x {params.index} \
                        -U {input.r1} \
                        --end-to-end \
                        --{params.mode} \
                        -p 1 -S "$tmp_sam" 2>> {log}
                    samtools sort -@ 0 -m {params.sort_mem} -o {output.bam} "$tmp_sam" 2>> {log}
                else
                    bowtie2 \
                        -x {params.index} \
                        -U {input.r1} \
                        --end-to-end \
                        --{params.mode} \
                        -p $MAIN_T 2>> {log} | \
                    samtools view -bS - | \
                    samtools sort -@ $SORT_T -m {params.sort_mem} -o {output.bam} -
                fi
                samtools index -@ $SORT_T {output.bam} >> {log} 2>&1
                samtools flagstat {output.bam} > {output.flagstat}
                """

    rule spikein_scale_factor:
        input:
            bam=lambda wc: spikein_bam_path(wc.run_sample),
            bai=lambda wc: spikein_bai_path(wc.run_sample),
            implementation=["workflow/rules/scripts/spikein_templates.py", "workflow/rules/scripts/fragment_records.py"]
        output:
            tsv=spikein_scale_factor_path("{run_sample}")
        params:
            read_type=READ_TYPE,
            scale_to=config.get("spikein", {}).get("scale_to", 10000),
            min_fragments=config.get("spikein", {}).get("min_fragments", 1000),
            exclude_flags=FRAGMENT_EXCLUDE_FLAGS
        log:
            f"{LOG_DIR}/spikein/{{run_sample}}.scale_factor.log"
        threads: int(config["threads"]["samtools"])
        conda:
            "envs/bowtie2.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log})

            python workflow/rules/scripts/spikein_templates.py \
                --bam {input.bam:q} --read-type {params.read_type} \
                --exclude-flags {params.exclude_flags} --threads {threads} \
                --sample {wildcards.run_sample:q} --scale-to {params.scale_to} \
                --min-fragments {params.min_fragments} --out {output.tsv:q} > {log:q} 2>&1
            """

    rule summarize_spikein_scale_factors:
        input:
            [spikein_scale_factor_path(run_sample) for run_sample in RUN_SAMPLES]
        output:
            tsv=spikein_scale_factors_path()
        conda:
            "envs/bowtie2.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.tsv})
            awk 'FNR == 1 && NR != 1 {{next}} {{print}}' {input} > {output.tsv}
            """
