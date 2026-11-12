# Quantitative tracks consume the same fragments as MACS3 and FRiP.
# The all-alignment TE diagnostic remains intentionally read-based.
def bigwig_params():
    return {
        "chrom_sizes": config["reference"]["chrom_sizes"],
        "bin_size": config.get("bigwig", {}).get("bin_size", 10),
        "normalization": config.get("bigwig", {}).get("normalization", "RPKM"),
        "exclude_flags": config.get("bigwig", {}).get("samtools_exclude_flags", config.get("samtools_exclude_flags", 1804)),
        "exclude_bed_cmd": exclude_bed_commands(),
        "remove_chrm": "--remove-chrm" if FRAGMENT_REMOVE_CHRM else "",
    }


def fragment_bigwig_shell(normalized=True, spikein=False):
    scale = "--scale-table {input.scale:q}" if spikein else ""
    normalization = "{params.normalization}" if normalized else "None"
    return r"""
        set -euo pipefail
        mkdir -p $(dirname {log:q})
        python workflow/rules/scripts/fragment_coverage.py \
            --fragments {input.fragments:q} --chrom-sizes {input.chrom_sizes:q} \
            --out {output.bw:q} --bin-size {params.bin_size} \
            --normalization """ + normalization + " " + scale + r""" \
            {params.remove_chrm} > {log:q} 2>&1
        """


if GENERATE_RAW_BIGWIG:
    rule bam_to_raw_bigwig:
        input:
            fragments=lambda wc: fragments_path(wc.run_sample),
            chrom_sizes=config["reference"]["chrom_sizes"],
            implementation=["workflow/rules/scripts/fragment_coverage.py", "workflow/rules/scripts/fragment_records.py"]
        params:
            **bigwig_params()
        output:
            bw=raw_bigwig_path("{run_sample}")
        log:
            f"{LOG_DIR}/bigwig/{{run_sample}}.fragments.raw.log"
        threads: 1
        resources:
            mem_mb=4000
        conda:
            "envs/deeptools.yaml"
        shell: fragment_bigwig_shell(normalized=False)


if GENERATE_NORMALIZED_BIGWIG:
    rule bam_to_normalized_bigwig:
        input:
            fragments=lambda wc: fragments_path(wc.run_sample),
            chrom_sizes=config["reference"]["chrom_sizes"],
            implementation=["workflow/rules/scripts/fragment_coverage.py", "workflow/rules/scripts/fragment_records.py"]
        params:
            **bigwig_params()
        output:
            bw=normalized_bigwig_path("{run_sample}")
        log:
            f"{LOG_DIR}/bigwig/{{run_sample}}.fragments.normalized.log"
        threads: 1
        resources:
            mem_mb=4000
        conda:
            "envs/deeptools.yaml"
        shell: fragment_bigwig_shell()

    rule pooled_bam_to_bigwig:
        wildcard_constraints:
            sample=POOLED_SAMPLE_PATTERN
        input:
            fragments=lambda wc: pooled_fragments_path(wc.sample),
            chrom_sizes=config["reference"]["chrom_sizes"],
            implementation=["workflow/rules/scripts/fragment_coverage.py", "workflow/rules/scripts/fragment_records.py"]
        params:
            **bigwig_params()
        output:
            bw=pooled_normalized_bigwig_path("{sample}")
        log:
            f"{LOG_DIR}/bigwig/{{sample}}.pooled.fragments.log"
        threads: 1
        resources:
            mem_mb=4000
        conda:
            "envs/deeptools.yaml"
        shell: fragment_bigwig_shell()

if TE_ENABLED and TE_GENERATE_BIGWIG:
    rule te_diagnostic_bigwig:
        """DIAGNOSTIC ONLY: unweighted coverage of the all-alignment TE BAM
        (multimappers counted once per alignment). Not a quantitative track;
        quantitative TE results come from te_family_enrichment (1/N weighted)."""
        input:
            bam=te_bam_path("{run_sample}"),
            bai=te_bai_path("{run_sample}")
        params:
            **bigwig_params()
        output:
            bw=te_diagnostic_bigwig_path("{run_sample}")
        log:
            f"{LOG_DIR}/bigwig/{{run_sample}}.te_diagnostic.log"
        threads: int(config["threads"]["deeptools"])
        resources:
            mem_mb=8000,
            disk_mb=20000
        conda:
            "envs/deeptools.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.bw}) $(dirname {log})
            EXCLUDE_BED=""
            {params.exclude_bed_cmd}
            EXTRA_BAMCOVERAGE_ARGS=""
            if [ -n "$EXCLUDE_BED" ]; then
                EXTRA_BAMCOVERAGE_ARGS="--blackListFileName $EXCLUDE_BED"
            fi
            bamCoverage \
                --bam {input.bam} \
                --outFileName {output.bw} \
                --outFileFormat bigwig \
                --numberOfProcessors {threads} \
                --binSize {params.bin_size} \
                --normalizeUsing None \
                --samFlagExclude {params.exclude_flags} \
                $EXTRA_BAMCOVERAGE_ARGS \
                > {log} 2>&1
            if [ -n "$EXCLUDE_BED" ]; then rm -f "$EXCLUDE_BED"; fi
            """



if SPIKEIN_BIGWIG_ENABLED:
    rule bam_to_spikein_normalized_bigwig:
        input:
            fragments=lambda wc: fragments_path(wc.run_sample),
            chrom_sizes=config["reference"]["chrom_sizes"],
            scale=spikein_scale_factor_path("{run_sample}"),
            implementation=["workflow/rules/scripts/fragment_coverage.py", "workflow/rules/scripts/fragment_records.py"]
        params:
            **bigwig_params()
        output:
            bw=spikein_normalized_bigwig_path("{run_sample}")
        log:
            f"{LOG_DIR}/bigwig/{{run_sample}}.spikein_fragments.log"
        threads: 1
        conda:
            "envs/deeptools.yaml"
        shell: fragment_bigwig_shell(normalized=False, spikein=True)
