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
        macs3_xls = "results/macs3/{sample}/{sample}_peaks.xls"
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
        """