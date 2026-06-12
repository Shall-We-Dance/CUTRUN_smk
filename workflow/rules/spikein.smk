# workflow/rules/spikein.smk

if SPIKEIN_ENABLED:
    if not is_single_end():
        rule bowtie2_spikein_pe:
            input:
                r1=f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R1.fastq.gz",
                r2=f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R2.fastq.gz"
            output:
                bam=spikein_bam_path("{sample}"),
                bai=spikein_bai_path("{sample}"),
                flagstat=spikein_flagstat_path("{sample}")
            params:
                index=config["spikein"]["bowtie2_index"],
                max_insert=config.get("bowtie2", {}).get("max_insert_size", 500),
                mode=config.get("bowtie2", {}).get("mode", "very-sensitive")
            log:
                "logs/spikein/{sample}_pe.log"
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
        rule bowtie2_spikein_se:
            input:
                r1=f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R1.fastq.gz"
            output:
                bam=spikein_bam_path("{sample}"),
                bai=spikein_bai_path("{sample}"),
                flagstat=spikein_flagstat_path("{sample}")
            params:
                index=config["spikein"]["bowtie2_index"],
                mode=config.get("bowtie2", {}).get("mode", "very-sensitive")
            log:
                "logs/spikein/{sample}_se.log"
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

    rule spikein_scale_factor:
        input:
            bam=lambda wc: spikein_bam_path(wc.sample),
            bai=lambda wc: spikein_bai_path(wc.sample)
        output:
            tsv=spikein_scale_factor_path("{sample}")
        params:
            read_type=READ_TYPE,
            scale_to=config.get("spikein", {}).get("scale_to", 10000),
            min_fragments=config.get("spikein", {}).get("min_fragments", 1000)
        log:
            "logs/spikein/{sample}.scale_factor.log"
        threads: int(config["threads"]["samtools"])
        conda:
            "envs/bowtie2.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log})

            if [ "{params.read_type}" = "PE" ]; then
                spikein_reads=$(samtools view -c -f 2 -F 4 {input.bam} 2>> {log})
                spikein_fragments=$((spikein_reads / 2))
            else
                spikein_reads=$(samtools view -c -F 4 {input.bam} 2>> {log})
                spikein_fragments="$spikein_reads"
            fi

            scale_factor=$(awk -v scale_to="{params.scale_to}" -v fragments="$spikein_fragments" 'BEGIN {if (fragments > 0) printf "%.10f", scale_to / fragments; else printf "1"}')
            status=$(awk -v fragments="$spikein_fragments" -v min_fragments="{params.min_fragments}" 'BEGIN {if (fragments == 0) print "NO_SPIKEIN"; else if (fragments < min_fragments) print "LOW_SPIKEIN"; else print "OK"}')

            {
                printf "sample\tspikein_reads\tspikein_fragments\tscale_to\tscale_factor\tstatus\n"
                printf "{wildcards.sample}\t%s\t%s\t{params.scale_to}\t%s\t%s\n" "$spikein_reads" "$spikein_fragments" "$scale_factor" "$status"
            } > {output.tsv}
            """

    rule summarize_spikein_scale_factors:
        input:
            [spikein_scale_factor_path(sample) for sample in SAMPLES]
        output:
            tsv=spikein_scale_factors_path()
        conda:
            "envs/bowtie2.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.tsv})
            awk 'FNR == 1 && NR != 1 {next} {print}' {input} > {output.tsv}
            """
