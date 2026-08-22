# rules/bam_to_bigwig.smk

BIGWIG_EXCLUDE_FLAGS = int(
    config.get("bigwig", {}).get(
        "samtools_exclude_flags",
        config.get("samtools_exclude_flags", 1804),
    )
)
BIGWIG_REMOVE_CHRM = bool(
    config.get("bigwig", {}).get("remove_chrM_and_scaffolds", True)
)


rule compute_genome_size:
    """Write the effective genome size over the contigs included in tracks."""
    output:
        txt=genome_size_path()
    conda:
        "envs/samtools.yaml"
    shell:
        effective_genome_size_cmd()


def bigwig_params():
    return {
        "chrom_sizes": config["reference"]["chrom_sizes"],
        "bin_size": config.get("bigwig", {}).get("bin_size", 10),
        "normalization": config.get("bigwig", {}).get("normalization", "RPKM"),
        "exclude_flags": BIGWIG_EXCLUDE_FLAGS,
        "remove_chrM_and_scaffolds": BIGWIG_REMOVE_CHRM,
        "exclude_bed_cmd": exclude_bed_commands(),
    }


def bigwig_shell():
    return r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bw}) $(dirname {log})

        GENOME_SIZE=$(cat {input.genome_size})
        if [ -z "$GENOME_SIZE" ] || [ "$GENOME_SIZE" = "0" ]; then
            echo "effective genome size is empty; falling back to sum of chrom sizes" >&2
            GENOME_SIZE=$(awk '{{sum += $2}} END {{print sum}}' {params.chrom_sizes})
        fi

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
            --normalizeUsing {params.normalization} \
            --samFlagExclude {params.exclude_flags} \
            --effectiveGenomeSize $GENOME_SIZE \
            $EXTRA_BAMCOVERAGE_ARGS \
            > {log} 2>&1

        if [ -n "$EXCLUDE_BED" ]; then rm -f "$EXCLUDE_BED"; fi
        """


if GENERATE_RAW_BIGWIG:
    rule bam_to_raw_bigwig:
        """Unnormalized coverage (read counts per bin)."""
        input:
            bam=lambda wc: final_bam_path(wc.run_sample),
            bai=lambda wc: final_bai_path(wc.run_sample),
            genome_size=genome_size_path()
        params:
            **bigwig_params()
        output:
            bw=raw_bigwig_path("{run_sample}")
        log:
            f"{LOG_DIR}/bigwig/{{run_sample}}.bamCoverage.raw.log"
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


if GENERATE_NORMALIZED_BIGWIG:
    rule bam_to_normalized_bigwig:
        """RPKG/RPKM-normalized coverage using the reference-correct effective
        genome size (computed over the contigs kept in the track)."""
        input:
            bam=lambda wc: final_bam_path(wc.run_sample),
            bai=lambda wc: final_bai_path(wc.run_sample),
            genome_size=genome_size_path()
        params:
            **bigwig_params()
        output:
            bw=normalized_bigwig_path("{run_sample}")
        log:
            f"{LOG_DIR}/bigwig/{{run_sample}}.bamCoverage.log"
        threads: int(config["threads"]["deeptools"])
        resources:
            mem_mb=8000,
            disk_mb=20000
        conda:
            "envs/deeptools.yaml"
        shell: bigwig_shell()


if GENERATE_NORMALIZED_BIGWIG and len(RUN_SAMPLES) > 0:
    rule pooled_bam_to_bigwig:
        """Normalized coverage of the pooled (replicate-merged) final BAM."""
        wildcard_constraints:
            sample=POOLED_SAMPLE_PATTERN
        input:
            bam=lambda wc: pooled_final_bam_path(wc.sample),
            bai=lambda wc: pooled_final_bai_path(wc.sample),
            genome_size=genome_size_path()
        params:
            **bigwig_params()
        output:
            bw=pooled_normalized_bigwig_path("{sample}")
        log:
            f"{LOG_DIR}/bigwig/{{sample}}.pooled.bamCoverage.log"
        threads: int(config["threads"]["deeptools"])
        resources:
            mem_mb=8000,
            disk_mb=20000
        conda:
            "envs/deeptools.yaml"
        shell: bigwig_shell()


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
        """Generate target-genome bigWig tracks scaled by spike-in alignment depth."""
        input:
            bam=lambda wc: final_bam_path(wc.run_sample),
            bai=lambda wc: final_bai_path(wc.run_sample),
            scale=spikein_scale_factor_path("{run_sample}")
        params:
            **bigwig_params()
        output:
            bw=spikein_normalized_bigwig_path("{run_sample}")
        log:
            f"{LOG_DIR}/bigwig/{{run_sample}}.spikein_bamCoverage.log"
        threads: int(config["threads"]["deeptools"])
        conda:
            "envs/deeptools.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.bw}) $(dirname {log})

            scale_factor=$(awk 'NR == 2 {{print $5}}' {input.scale})
            if [ -z "$scale_factor" ]; then
                scale_factor=1
            fi

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
                --scaleFactor "$scale_factor" \
                --samFlagExclude {params.exclude_flags} \
                $EXTRA_BAMCOVERAGE_ARGS \
                > {log} 2>&1

            if [ -n "$EXCLUDE_BED" ]; then
                rm -f "$EXCLUDE_BED"
            fi
            """
