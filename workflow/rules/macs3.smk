# rules/macs3.smk

rule macs3_callpeak:
    """Call peaks on final BAM (post-filtering)."""
    input:
        bam = lambda wc: final_bam_path(wc.sample)
    output:
        xls = f"{OUTDIR}/macs3/{{sample}}/{{sample}}_peaks.xls",
        narrow = f"{OUTDIR}/macs3/{{sample}}/{{sample}}_peaks.narrowPeak",
        summits = f"{OUTDIR}/macs3/{{sample}}/{{sample}}_summits.bed"
    params:
        outdir = f"{OUTDIR}/macs3/{{sample}}",
        name = "{sample}",
        gsize = config["macs"]["gsize"],
        qvalue = config["macs"]["qvalue"],
        extsize = config["macs"]["extsize"],
        fmt = "BAM" if is_single_end() else "BAMPE",
        broad = "--broad --broad-cutoff " + str(config["macs"].get("broad_cutoff", 0.1))
                if config["macs"].get("broad", False) else ""
    log:
        f"logs/macs3/{{sample}}.log"
    threads: 1
    conda:
        "envs/macs3.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        macs3 callpeak \
            -t {input.bam} \
            -f {params.fmt} \
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
