# workflow/rules/qc_fastp.smk
import os

OUTDIR = config["output"]["dir"]


def units(sample):
    return list(range(len(config["samples"][sample]["R1"])))


def fastp_extra_args():
    dedup_adapter = config.get("fastp", {}).get("dedup_adapter", {})
    if not bool(dedup_adapter.get("enable", True)):
        return "--disable_adapter_trimming"

    args = []
    if bool(dedup_adapter.get("dedup", True)):
        args.append("--dedup")
    if dedup_adapter.get("extra"):
        args.append(str(dedup_adapter["extra"]))
    return " ".join(args)


def multiqc_search_paths():
    paths = [
        f"{OUTDIR}/qc/fastp",
        f"{OUTDIR}/qc/bam_stats",
        f"{OUTDIR}/bowtie2",
        "logs/fastp",
        "logs/bowtie2",
    ]
    if CUTRUN_QC_ENABLED:
        paths.extend([cutrun_qc_summary_path(), cutrun_fragment_lengths_multiqc_path(), "logs/cutrun_qc"])
    if SPIKEIN_ENABLED:
        paths.extend([f"{OUTDIR}/qc/spikein", "logs/spikein"])
    return " ".join(paths)


# ----------------------------
# Merge raw FASTQ per sample first
# ----------------------------
if not is_single_end():
    rule merge_raw_fastq_per_sample:
        input:
            r1=lambda wc: [config["samples"][wc.sample]["R1"][i] for i in units(wc.sample)],
            r2=lambda wc: [config["samples"][wc.sample]["R2"][i] for i in units(wc.sample)]
        output:
            merged_r1=temp(f"{OUTDIR}/tmp/merged_raw/{{sample}}_R1.fastq.gz"),
            merged_r2=temp(f"{OUTDIR}/tmp/merged_raw/{{sample}}_R2.fastq.gz")
        threads: 2
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.merged_r1})
            cat {input.r1} > {output.merged_r1}
            cat {input.r2} > {output.merged_r2}
            """
else:
    rule merge_raw_fastq_per_sample:
        input:
            r1=lambda wc: [config["samples"][wc.sample]["R1"][i] for i in units(wc.sample)]
        output:
            merged_r1=temp(f"{OUTDIR}/tmp/merged_raw/{{sample}}_R1.fastq.gz")
        threads: 2
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.merged_r1})
            cat {input.r1} > {output.merged_r1}
            """


# ----------------------------
# Run fastp at sample level (after merge)
# ----------------------------
if not is_single_end():
    rule fastp_sample_level:
        input:
            r1=f"{OUTDIR}/tmp/merged_raw/{{sample}}_R1.fastq.gz",
            r2=f"{OUTDIR}/tmp/merged_raw/{{sample}}_R2.fastq.gz"
        output:
            clean_r1=temp(f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R1.fastq.gz"),
            clean_r2=temp(f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R2.fastq.gz"),
            html=f"{OUTDIR}/qc/fastp/{{sample}}/fastp.html",
            json=f"{OUTDIR}/qc/fastp/{{sample}}/fastp.json"
        log:
            f"logs/fastp/{{sample}}.log"
        threads: int(config["threads"]["fastp"])
        conda:
            "envs/qc.yaml"
        params:
            extra_args=fastp_extra_args()
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
else:
    rule fastp_sample_level:
        input:
            r1=f"{OUTDIR}/tmp/merged_raw/{{sample}}_R1.fastq.gz"
        output:
            clean_r1=temp(f"{OUTDIR}/tmp/fastp_sample/{{sample}}_R1.fastq.gz"),
            html=f"{OUTDIR}/qc/fastp/{{sample}}/fastp.html",
            json=f"{OUTDIR}/qc/fastp/{{sample}}/fastp.json"
        log:
            f"logs/fastp/{{sample}}.log"
        threads: int(config["threads"]["fastp"])
        conda:
            "envs/qc.yaml"
        params:
            extra_args=fastp_extra_args()
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


# ----------------------------
# MultiQC summary
# ----------------------------
rule multiqc:
    input:
        expand(f"{OUTDIR}/qc/fastp/{{sample}}/fastp.html", sample=SAMPLES),
        [bowtie2_log_path(sample) for sample in SAMPLES],
        [path for sample in SAMPLES for path in bam_flagstat_paths(sample)],
        [duplicate_metrics_path(sample) for sample in SAMPLES] if REMOVE_DUPLICATES else [],
        [path for sample in SAMPLES for path in cutrun_qc_paths(sample)],
        [cutrun_qc_summary_path(), cutrun_fragment_lengths_multiqc_path()] if CUTRUN_QC_ENABLED else [],
        [path for sample in SAMPLES for path in spikein_qc_paths(sample)],
        [spikein_scale_factors_path()] if SPIKEIN_ENABLED else []
    output:
        html=f"{OUTDIR}/qc/multiqc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    threads: 1
    conda:
        "envs/multiqc.yaml"
    params:
        search_paths=multiqc_search_paths()
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.html}) $(dirname {log})
        multiqc --force -o {OUTDIR}/qc/multiqc {params.search_paths} > {log} 2>&1
        """
