# rules/bam_to_bigwig.smk

rule bam_to_bigwig:
    """Generate bigWig coverage tracks from raw BAM files."""
    input:
        bam = lambda wc: final_bam_path(wc.sample),
        bai = lambda wc: final_bai_path(wc.sample)
    params:
        chrom_sizes = config["reference"]["chrom_sizes"],
        bin_size = config.get("bigwig", {}).get("bin_size", 10),
        normalization = config.get("bigwig", {}).get("normalization", "RPKM"),
        exclude_flags = config.get("samtools_exclude_flags", 1804)
    output:
        raw_bw = f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.raw.bw",
        norm_bw = f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.normalized.bw"
    log:
        f"logs/bigwig/{{sample}}.bamCoverage.log"
    threads: int(config["threads"]["deeptools"])
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        GENOME_SIZE=$(awk '{{sum += $2}} END {{print sum}}' {params.chrom_sizes})
        mkdir -p $(dirname {output.raw_bw})
        
        # Unnormalized coverage
        bamCoverage \
            --bam {input.bam} \
            --outFileName {output.raw_bw} \
            --outFileFormat bigwig \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --normalizeUsing None \
            --samFlagExclude {params.exclude_flags} \
            --effectiveGenomeSize $GENOME_SIZE \
            > {log} 2>&1
        
        # Normalized coverage
        bamCoverage \
            --bam {input.bam} \
            --outFileName {output.norm_bw} \
            --outFileFormat bigwig \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --normalizeUsing {params.normalization} \
            --samFlagExclude {params.exclude_flags} \
            --effectiveGenomeSize $GENOME_SIZE \
            >> {log} 2>&1
        """
