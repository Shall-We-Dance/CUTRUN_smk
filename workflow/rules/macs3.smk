# rules/macs3.smk

def get_macs3_input_file_type(sample_name):
    if sample_name in PE_SAMPLES:
        return "BAMPE"
    elif sample_name in SE_SAMPLES:
        return "BAM"
    else:
        raise ValueError(f"Sample {sample_name} not found in PE_SAMPLES or SE_SAMPLES.")

rule macs3_callpeak_unique:
    input:
        bam = lambda wildcards: get_bam_path(wildcards.sample)
    output:
        macs3_xls = "results/macs3/{sample}/{sample}_peaks.xls",
        narrow_peaks = "results/macs3/{sample}/{sample}_peaks.narrowPeak",
        summits = "results/macs3/{sample}/{sample}_summits.bed"
    params:
        outdir = "results/macs3/{sample}",
        name = "{sample}",
        tempdir = temp("results/macs3/{sample}/temp"),
        gsize = config["macs_gsize"],
        qvalue = config["macs_qvalue"],
        extsize = config["macs_extsize"],
        input_type = lambda wildcards: get_macs3_input_file_type(wildcards.sample)
    log:
        "logs/macs3/{sample}_macs3.log"
    conda:
        "../envs/macs3.yaml"
    threads: 1
    shell:
        """
        mkdir -p {params.tempdir}
        macs3 callpeak -f {params.input_type} \
                       -t {input.bam} \
                       -g {params.gsize} \
                       --outdir {params.outdir} \
                       -n {params.name} \
                       --nomodel \
                       --nolambda \
                       --extsize {params.extsize} \
                       --keep-dup 1 \
                       -q {params.qvalue} \
                       --tempdir {params.tempdir} \
                       --call-summits > {log} 2>&1
        rm -rf {params.tempdir}
        echo "MACS3 peak calling completed for {wildcards.sample}." >> {log}
        """

rule macs3_filter_peaks:
    input:
        narrow_peaks = "results/macs3/{sample}/{sample}_peaks.narrowPeak",
        summits = "results/macs3/{sample}/{sample}_summits.bed"
    output:
        filtered_peaks = "results/macs3/{sample}/{sample}_peaks_filtered.narrowPeak",
        filtered_summits = "results/macs3/{sample}/{sample}_summits_filtered.bed"
    params:
        blacklist = lambda wildcards: get_blacklist_path(config["genome"])
    log:
        "logs/macs3/{sample}_macs3_filter.log"
    shell:
        """
        bedtools intersect -v -a {input.summits} -b {params.blacklist} -nonamecheck > {output.filtered_summits}
        bedtools intersect -v -a {input.narrow_peaks} -b {params.blacklist} -nonamecheck > {output.filtered_peaks}
        echo "Filtered peaks and summits saved to {output.filtered_peaks} and {output.filtered_summits}." >> {log}
        """
