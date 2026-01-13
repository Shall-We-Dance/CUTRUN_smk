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
        bam = lambda wc: get_raw_bam(wc.sample)
    output:
        xls = "results/macs3/{sample}/{sample}_peaks.xls",
        narrow = "results/macs3/{sample}/{sample}_peaks.narrowPeak",
        summits = "results/macs3/{sample}/{sample}_summits.bed"
    params:
        outdir = "results/macs3/{sample}",
        name = "{sample}",
        gsize = config["macs_gsize"],
        qvalue = config["macs_qvalue"],
        extsize = config["macs_extsize"],
        input_type = lambda wc: get_macs3_input_file_type(wc.sample)
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
            --call-summits >> {log} 2>&1
        """

rule macs3_callpeak_unique_blacklist:
    input:
        bam = lambda wc: get_filtered_bam(wc.sample)
    output:
        xls = "results/macs3/{sample}/{sample}.filtered_peaks.xls",
        narrow = "results/macs3/{sample}/{sample}.filtered_peaks.narrowPeak",
        summits = "results/macs3/{sample}/{sample}.filtered_summits.bed"
    params:
        outdir = "results/macs3/{sample}",
        name = "{sample}.filtered",
        gsize = config["macs_gsize"],
        qvalue = config["macs_qvalue"],
        extsize = config["macs_extsize"],
        input_type = lambda wc: get_macs3_input_file_type(wc.sample)
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
            --call-summits >> {log} 2>&1
        """
