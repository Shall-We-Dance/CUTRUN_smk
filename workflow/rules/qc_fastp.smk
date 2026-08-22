# workflow/rules/qc_fastp.smk
import os


def run_sample_parts_impl(run_sample):
    return run_sample_parts(run_sample)


def units(run_sample):
    sample, rep = run_sample_parts_impl(run_sample)
    return list(range(len(config["samples"][sample]["replicates"][rep]["R1"])))


def rep_fastqs(run_sample, key):
    sample, rep = run_sample_parts_impl(run_sample)
    return config["samples"][sample]["replicates"][rep][key]


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


def fastp_html(run_sample):
    return f"{OUTDIR}/qc/fastp/{run_sample}/{run_sample}.fastp.html"


def fastp_json(run_sample):
    return f"{OUTDIR}/qc/fastp/{run_sample}/{run_sample}.fastp.json"


def merged_raw_r1(run_sample):
    return f"{OUTDIR}/tmp/merged_raw/{run_sample}_R1.fastq.gz"


def merged_raw_r2(run_sample):
    return f"{OUTDIR}/tmp/merged_raw/{run_sample}_R2.fastq.gz"


def fastp_sample_r1(run_sample):
    return f"{OUTDIR}/tmp/fastp_sample/{run_sample}_R1.fastq.gz"


def fastp_sample_r2(run_sample):
    return f"{OUTDIR}/tmp/fastp_sample/{run_sample}_R2.fastq.gz"


def multiqc_search_paths():
    # Give MultiQC exact report/metrics files. Scanning the complete bowtie2
    # output directory also scans large BAMs and can make report generation
    # needlessly slow, especially on networked storage.
    paths = [fastp_json(run_sample) for run_sample in RUN_SAMPLES]
    paths.extend(bowtie2_log_path(run_sample) for run_sample in RUN_SAMPLES)
    paths.extend(
        path
        for run_sample in RUN_SAMPLES
        for path in bam_flagstat_paths(run_sample)
    )
    if TE_ENABLED:
        paths.extend(te_flagstat_path(run_sample) for run_sample in RUN_SAMPLES)
    if REMOVE_DUPLICATES:
        paths.extend(duplicate_metrics_path(run_sample) for run_sample in RUN_SAMPLES)
    if CUTRUN_QC_ENABLED:
        paths.extend(
            [cutrun_summary_multiqc_path(), cutrun_fragment_lengths_multiqc_path()]
        )
    if SPIKEIN_ENABLED:
        paths.extend(spikein_flagstat_path(run_sample) for run_sample in RUN_SAMPLES)
    if FRIP_QC_ENABLED:
        paths.append(frip_multiqc_path())
    if TE_OVERLAP_QC_ENABLED:
        paths.append(te_overlap_mqc_path())
    return paths


# ----------------------------
# Merge raw FASTQ per replicate (multi-lane support; a single biological
# replicate may consist of several sequencing lanes)
# ----------------------------
if not is_single_end():
    rule merge_raw_fastq_per_replicate:
        input:
            r1=lambda wc: rep_fastqs(wc.run_sample, "R1"),
            r2=lambda wc: rep_fastqs(wc.run_sample, "R2")
        output:
            merged_r1=temp(merged_raw_r1("{run_sample}")),
            merged_r2=temp(merged_raw_r2("{run_sample}"))
        threads: 2
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.merged_r1})
            cat {input.r1} > {output.merged_r1}
            cat {input.r2} > {output.merged_r2}
            """
else:
    rule merge_raw_fastq_per_replicate:
        input:
            r1=lambda wc: rep_fastqs(wc.run_sample, "R1")
        output:
            merged_r1=temp(merged_raw_r1("{run_sample}"))
        threads: 2
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.merged_r1})
            cat {input.r1} > {output.merged_r1}
            """


# ----------------------------
# fastp per replicate
# ----------------------------
if not is_single_end():
    rule fastp_sample_level:
        input:
            r1=merged_raw_r1("{run_sample}"),
            r2=merged_raw_r2("{run_sample}")
        output:
            clean_r1=temp(fastp_sample_r1("{run_sample}")),
            clean_r2=temp(fastp_sample_r2("{run_sample}")),
            html=fastp_html("{run_sample}"),
            json=fastp_json("{run_sample}")
        log:
            f"{LOG_DIR}/fastp/{{run_sample}}.log"
        threads: int(config["threads"]["fastp"])
        resources:
            mem_mb=4000,
            disk_mb=50000
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
            r1=merged_raw_r1("{run_sample}")
        output:
            clean_r1=temp(fastp_sample_r1("{run_sample}")),
            html=fastp_html("{run_sample}"),
            json=fastp_json("{run_sample}")
        log:
            f"{LOG_DIR}/fastp/{{run_sample}}.log"
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
        expand(fastp_json("{run_sample}"), run_sample=RUN_SAMPLES),
        [bowtie2_log_path(run_sample) for run_sample in RUN_SAMPLES],
        [path for run_sample in RUN_SAMPLES for path in bam_flagstat_paths(run_sample)],
        [te_flagstat_path(run_sample) for run_sample in RUN_SAMPLES] if TE_ENABLED else [],
        [duplicate_metrics_path(run_sample) for run_sample in RUN_SAMPLES] if REMOVE_DUPLICATES else [],
        [path for run_sample in RUN_SAMPLES for path in cutrun_qc_paths(run_sample)],
        [cutrun_qc_summary_path(), cutrun_summary_multiqc_path(), cutrun_fragment_lengths_multiqc_path()] if CUTRUN_QC_ENABLED else [],
        [path for run_sample in RUN_SAMPLES for path in spikein_qc_paths(run_sample)],
        [spikein_scale_factors_path()] if SPIKEIN_ENABLED else [],
        [frip_summary_path(), frip_multiqc_path()] if FRIP_QC_ENABLED else [],
        [te_overlap_summary_path(), te_overlap_mqc_path()] if TE_OVERLAP_QC_ENABLED else []
    output:
        html=f"{OUTDIR}/qc/multiqc/multiqc_report.html"
    log:
        f"{LOG_DIR}/multiqc.log"
    threads: 1
    conda:
        "envs/multiqc.yaml"
    params:
        search_paths=multiqc_search_paths()
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.html}) $(dirname {log})
        multiqc --force -o {OUTDIR}/qc/multiqc {params.search_paths:q} > {log} 2>&1

        # Fail-fast: a broken custom-content module makes MultiQC exit 0 but
        # silently drops the section. Check for tracebacks and required IDs.
        if grep -qE "Oops! The 'custom_content'|AssertionError|Traceback \(most recent" {log}; then
            echo "MultiQC custom-content module failed; see {log}" >&2
            exit 1
        fi
        if [ ! -d {OUTDIR}/qc/multiqc/multiqc_data ]; then
            echo "MultiQC multiqc_data directory missing" >&2
            exit 1
        fi
        if [ "{CUTRUN_QC_ENABLED}" = "True" ]; then
            for section_id in cutrun_qc_table cutrun_fragment_lengths_plot; do
                if ! grep -rq "\"$section_id\"" {OUTDIR}/qc/multiqc/multiqc_data/; then
                    echo "MultiQC report missing required '$section_id' section" >&2
                    exit 1
                fi
            done
        fi
        if [ "{FRIP_QC_ENABLED}" = "True" ]; then
            if ! grep -rq '"frip_table"' {OUTDIR}/qc/multiqc/multiqc_data/; then
                echo "MultiQC report missing required 'frip_table' section" >&2
                exit 1
            fi
        fi
        if [ "{TE_OVERLAP_QC_ENABLED}" = "True" ]; then
            if ! grep -rq '"te_overlap_table"' {OUTDIR}/qc/multiqc/multiqc_data/; then
                echo "MultiQC report missing required 'te_overlap_table' section" >&2
                exit 1
            fi
        fi
        """
