# rules/macs3.smk

if MACS3_ENABLED:
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
            gsize = MACS3_CONFIG["gsize"],
            qvalue = MACS3_CONFIG["qvalue"],
            extsize = MACS3_CONFIG["extsize"],
            fmt = "BAM" if is_single_end() else "BAMPE",
            broad = "--broad --broad-cutoff " + str(MACS3_CONFIG.get("broad_cutoff", 0.1))
                    if MACS3_CONFIG.get("broad", False) else ""
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
