# rules/macs3.smk


rule macs3_callpeak:
    input:
        bam = lambda wildcards: get_bam_path(wildcards.sample)
    output:
        narrowpeak = "results/macs3/{sample}/{sample}_peaks.narrowPeak"
    params:
        outdir = "results/macs3/{sample}",
        name = "{sample}",
        tempdir = "results/macs3/{sample}/temp",
        gsize = config["macs_gsize"]
    log:
        "logs/macs3/{sample}_macs3.log"
    conda:
        "envs/macs3.yaml"
    threads: 1
    shell:
        """
        mkdir -p {params.tempdir}
        macs3 callpeak -f AUTO \
                       -t {input.bam} \
                       -g {params.gsize} \
                       --outdir {params.outdir} \
                       -n {params.name} \
                       --tempdir {params.tempdir} \
                       --call-summits > {log} 2>&1
        """