# workflow/rules/cutrun_qc.smk
# Counts/lengths are emitted by prepare_fragments, even when summary QC is off.

if CUTRUN_QC_ENABLED:
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
        """Fragment-level FRiP for a replicate's shared fragment table."""
        input:
            fragments=lambda wc: fragments_path(wc.run_sample),
            peaks=lambda wc: macs3_peak_path(wc.run_sample),
            implementation=["workflow/rules/scripts/frip.py", "workflow/rules/scripts/fragment_records.py"]
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
                --fragments {input.fragments:q} \
                --peaks {input.peaks} \
                --read-type {params.read_type} \
                --sample {wildcards.run_sample} \
                --threads {threads} \
                --out {output.tsv} > {log} 2>&1
            """

    rule frip_pooled:
        """Fragment-level FRiP for pooled fragments against pooled peaks."""
        wildcard_constraints:
            sample=POOLED_PEAK_SAMPLE_PATTERN
        input:
            fragments=lambda wc: pooled_fragments_path(wc.sample),
            peaks=lambda wc: macs3_peak_path(wc.sample, pooled=True),
            implementation=["workflow/rules/scripts/frip.py", "workflow/rules/scripts/fragment_records.py"]
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
                --fragments {input.fragments:q} \
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
