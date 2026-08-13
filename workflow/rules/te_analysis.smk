# workflow/rules/te_analysis.smk
#
# TE-mode analysis rules:
#  - te_quant: fractional (1/N) family-weighted quantification + overlap QC of
#    the TE-mode all-alignment BAM against the TE annotation BED
#  - summarize: cross-sample overlap summary (TSV + MultiQC JSON)
#  - enrichment: family-level ChIP-vs-matched-input enrichment (primary
#    quantitative TE result)

if TE_ENABLED and TE_ANNOTATION:

    rule te_overlap_qc:
        input:
            bam=te_bam_path("{run_sample}"),
            bai=te_bai_path("{run_sample}"),
            te_bed=TE_ANNOTATION
        output:
            tsv=te_overlap_path("{run_sample}"),
            family=te_family_counts_path("{run_sample}")
        params:
            exclude_flags=TE_EXCLUDE_FLAGS,
            read_type=READ_TYPE,
            family_column=TE_FAMILY_COLUMN,
            subfamily_column=TE_SUBFAMILY_COLUMN,
            class_column=TE_CLASS_COLUMN,
            min_overlap_frac=TE_MIN_OVERLAP_FRAC,
            min_overlap_bp=int(TE_CONFIG.get("min_overlap_bp", 1))
        log:
            f"{LOG_DIR}/te/{{run_sample}}.overlap.log"
        threads: int(config["threads"].get("samtools", 4))
        resources:
            mem_mb=16000,
            disk_mb=20000
        benchmark:
            f"{OUTDIR}/qc/benchmarks/{ANALYSIS_MODE}/te_overlap_qc/{{run_sample}}.benchmark.tsv"
        conda:
            "envs/cutrun_qc.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log})
            SORT_T=$(( {threads} - 1 )); [ "$SORT_T" -lt 1 ] && SORT_T=1
            python workflow/rules/scripts/te_quant.py \
                --bam {input.bam} \
                --te-bed {input.te_bed} \
                --exclude-flags {params.exclude_flags} \
                --read-type {params.read_type} \
                --family-column {params.family_column} \
                --subfamily-column {params.subfamily_column} \
                --class-column {params.class_column} \
                --min-overlap-frac {params.min_overlap_frac} \
                --min-overlap-bp {params.min_overlap_bp} \
                --sample {wildcards.run_sample} \
                --threads {threads} \
                --out-tsv {output.tsv} \
                --out-family {output.family} > {log} 2>&1
            """

    rule summarize_te_overlap:
        input:
            expand(te_overlap_path("{run_sample}"), run_sample=RUN_SAMPLES)
        output:
            tsv=te_overlap_summary_path(),
            mqc=te_overlap_mqc_path()
        conda:
            "envs/cutrun_qc.yaml"
        script:
            "scripts/summarize_te_overlap.py"

    rule te_family_enrichment:
        """Family-level weighted enrichment of each ChIP sample versus its
        matched input, normalized by template depth (RPM). Outputs pooled +
        per-replicate metrics, curated/novel family splits, and subclass
        aggregation."""
        input:
            family_counts=[te_family_counts_path(run_sample) for run_sample in RUN_SAMPLES]
        params:
            dfam_metadata=TE_CONFIG.get("dfam_metadata")
        output:
            tsv=te_enrichment_path(),
            curated=te_enrichment_path().replace(".tsv", ".curated.tsv"),
            novel=te_enrichment_path().replace(".tsv", ".novel.tsv"),
            unknown=te_enrichment_path().replace(".tsv", ".unknown.tsv"),
            subclass=te_enrichment_path().replace(".tsv", ".subclass.tsv")
        conda:
            "envs/cutrun_qc.yaml"
        script:
            "scripts/summarize_te_enrichment.py"

    if TE_BLACKLIST_MODE == "full" and FILTER_BLACKLIST:
        rule te_blacklist_lowconf:
            """Blacklist intervals overlapping annotated TE loci. These
            regions are NOT rescued automatically; the file marks them as
            low-confidence loci for sensitivity/mappability inspection."""
            input:
                blacklist=blacklist_path(),
                te=TE_ANNOTATION
            output:
                bed=te_blacklist_lowconf_path()
            log:
                f"{LOG_DIR}/te/te_blacklist_lowconf.log"
            conda:
                "envs/bowtie2.yaml"
            shell:
                """
                set -euo pipefail
                mkdir -p $(dirname {output.bed})
                bedtools intersect -wa -a {input.blacklist} -b {input.te} \
                    > {output.bed} 2> {log}
                echo "low-confidence (blacklist x TE) regions: $(wc -l < {output.bed})" >> {log}
                """
