# rules/bam_to_bigwig.smk

def get_bam_path(sample):
    if config.get("mask_repeats", True):
        return f"results/star/{sample}/{sample}_masked.bam"
    else:
        return f"results/star/{sample}/{sample}_Aligned.sortedByCoord.out.bam"

rule bam_to_bigwig:
    input:
        bam = lambda wildcards: get_bam_path(wildcards.sample)
    params:
        chrom_sizes = config["chrom_sizes"],
        bin_size = config["bw_bin_size"]
    output:
        raw_bw = "results/bigwig/{sample}.bw",
        rpkm_bw = "results/bigwig/{sample}.rpkm.bw"
    log:
        "logs/bigwig/{sample}.log"
    threads: config["threads"]
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        # Unnormalized
        bamCoverage --bam {input.bam} \
                    --outFileName {output.raw_bw} \
                    --outFileFormat bigwig \
                    --numberOfProcessors {threads} \
                    --binSize {params.bin_size} \
                    --normalizeUsing None \
                    --effectiveGenomeSize $(awk '{{sum += $2}} END {{print sum}}' {params.chrom_sizes}) \
                    > {log} 2>&1

        # RPKM normalized
        bamCoverage --bam {input.bam} \
                    --outFileName {output.rpkm_bw} \
                    --outFileFormat bigwig \
                    --numberOfProcessors {threads} \
                    --binSize {params.bin_size} \
                    --normalizeUsing RPKM \
                    --effectiveGenomeSize $(awk '{{sum += $2}} END {{print sum}}' {params.chrom_sizes}) \
                    > {log} 2>&1
        """
