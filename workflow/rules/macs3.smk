# rules/macs3.smk

def get_macs3_input_file_type(sample_name):
    """Determine input file type for MACS3 based on sample type."""
    if sample_name in PE_SAMPLES:
        return "BAMPE"
    elif sample_name in SE_SAMPLES:
        return "BAM"
    else:
        raise ValueError(f"Sample {sample_name} not found in PE_SAMPLES or SE_SAMPLES.")


def get_macs3_input_bam(sample_name):
    """Get the appropriate BAM file for MACS3 peak calling."""
    from rules.bowtie2 import get_analysis_bam
    return get_analysis_bam(sample_name)


rule macs3_callpeak_unique:
    """Call peaks on uniquely mapped reads (unfiltered)."""
    input:
        bam = lambda wc: get_raw_bam(wc.sample)
    output:
        xls = "results/macs3/{sample}/{sample}_peaks.xls",
        narrow = "results/macs3/{sample}/{sample}_peaks.narrowPeak",
        summits = "results/macs3/{sample}/{sample}_summits.bed"
    params:
        outdir = "results/macs3/{sample}",
        name = "{sample}",
        gsize = config["macs"]["gsize"],
        qvalue = config["macs"]["qvalue"],
        extsize = config["macs"]["extsize"],
        input_type = lambda wc: get_macs3_input_file_type(wc.sample),
        broad = "--broad --broad-cutoff " + str(config.get("macs_broad_cutoff", 0.1)) \
                if config.get("macs_broad", False) else ""
    log:
        "logs/macs3/{sample}.log"
    threads: 1
    conda:
        "../envs/macs3.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        macs3 callpeak \
            -t {input.bam} \
            -f {params.input_type} \
            -g {params.gsize} \
            -n {params.name} \
            --outdir {params.outdir} \
            --nomodel \
            --extsize {params.extsize} \
            --keep-dup 1 \
            -q {params.qvalue} \
            --call-summits \
            {params.broad} >> {log} 2>&1
        """


rule macs3_callpeak_filtered:
    """Call peaks on filtered/processed reads."""
    input:
        bam = lambda wc: get_macs3_input_bam(wc.sample)
    output:
        xls = "results/macs3/{sample}/{sample}.filtered_peaks.xls",
        narrow = "results/macs3/{sample}/{sample}.filtered_peaks.narrowPeak",
        summits = "results/macs3/{sample}/{sample}.filtered_summits.bed"
    params:
        outdir = "results/macs3/{sample}",
        name = "{sample}.filtered",
        gsize = config["macs"]["gsize"],
        qvalue = config["macs"]["qvalue"],
        extsize = config["macs"]["extsize"],
        input_type = lambda wc: get_macs3_input_file_type(wc.sample),
        broad = "--broad --broad-cutoff " + str(config.get("macs_broad_cutoff", 0.1)) \
                if config.get("macs_broad", False) else ""
    log:
        "logs/macs3/{sample}.filtered.log"
    threads: 1
    conda:
        "../envs/macs3.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        macs3 callpeak \
            -t {input.bam} \
            -f {params.input_type} \
            -g {params.gsize} \
            -n {params.name} \
            --outdir {params.outdir} \
            --nomodel \
            --extsize {params.extsize} \
            --keep-dup 1 \
            -q {params.qvalue} \
            --call-summits \
            {params.broad} >> {log} 2>&1
        """


rule annotate_peaks:
    """Annotate peaks with nearby genes (optional, requires ChIPseeker or HOMER)."""
    input:
        peaks = "results/macs3/{sample}/{sample}_peaks.narrowPeak"
    output:
        annotated = "results/macs3/{sample}/{sample}_peaks.annotated.txt"
    params:
        genome = config.get("genome", "hg38")
    log:
        "logs/macs3/{sample}.annotate.log"
    conda:
        "../envs/chipseeker.yaml"
    script:
        "../scripts/annotate_peaks.R"