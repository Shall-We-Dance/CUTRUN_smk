# rules/macs3.smk

if MACS3_ENABLED:
    def macs3_call_peaks_run_samples():
        return [
            run_sample
            for run_sample in RUN_SAMPLES
            if config["samples"][sample_of(run_sample)].get("call_peaks", True)
        ]

    def macs3_sample_params(sample, run_sample=None):
        """Per-sample MACS3 parameters (sample-level overrides on top of the
        global macs3 config).

        ``keep_dup`` default: the final BAM is already coordinate-deduped
        upstream (Picard/samtools) when ``remove_duplicates: true``, so MACS3
        must not apply a second coordinate-level dedup (--keep-dup 1) -- that
        would also wrongly collapse fragments from different biological
        replicates in pooled calls. Default is therefore "all"; keep the
        option for sensitivity runs.
        """
        merged = {}
        merged.update(MACS3_CONFIG)
        merged.update(config["samples"][sample].get("macs3", {}))
        if "keep_dup" not in merged:
            merged["keep_dup"] = "all" if REMOVE_DUPLICATES else "1"
        return merged

    def macs3_control_bam(wc, pooled=False):
        sample = sample_of(wc.run_sample) if not pooled else wc.sample
        control = config["samples"].get(sample, {}).get("control")
        if not control:
            return []
        if pooled:
            return pooled_final_bam_path(control)
        rep = replicate_of(wc.run_sample)
        control_rep = control_replicate_for(sample, rep)
        if control_rep == "pooled":
            return pooled_final_bam_path(control)
        return final_bam_path(f"{control}.{control_rep}")

    PEAK_RUN_SAMPLES = macs3_call_peaks_run_samples()
    PEAK_SAMPLES = sorted({sample_of(rs) for rs in PEAK_RUN_SAMPLES})

    rule macs3_callpeak:
        """Call peaks on the final BAM of one replicate, optionally with the
        matched INPUT/control BAM of the same replicate."""
        input:
            bam=lambda wc: final_bam_path(wc.run_sample),
            control=macs3_control_bam
        output:
            xls=macs3_xls_path("{run_sample}"),
            narrow=macs3_narrow_path("{run_sample}"),
            summits=macs3_summits_path("{run_sample}")
        params:
            outdir=lambda wc: f"{OUTDIR}/macs3/{wc.run_sample}",
            name="{run_sample}",
            gsize=lambda wc: macs3_sample_params(sample_of(wc.run_sample))["gsize"],
            qvalue=lambda wc: macs3_sample_params(sample_of(wc.run_sample))["qvalue"],
            extsize=lambda wc: macs3_sample_params(sample_of(wc.run_sample))["extsize"],
            fmt="BAM" if is_single_end() else "BAMPE",
            keep_dup=lambda wc: str(macs3_sample_params(sample_of(wc.run_sample)).get("keep_dup", "all")),
            nomodel=lambda wc: macs3_sample_params(sample_of(wc.run_sample)).get("nomodel", True),
            broad=lambda wc: "--broad --broad-cutoff " + str(macs3_sample_params(sample_of(wc.run_sample)).get("broad_cutoff", 0.1))
                    if macs3_sample_params(sample_of(wc.run_sample)).get("broad", False) else ""
        log:
            f"{LOG_DIR}/macs3/{{run_sample}}.log"
        threads: 2
        resources:
            mem_mb=8000,
            disk_mb=20000
        conda:
            "envs/macs3.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p {params.outdir} $(dirname {log})
            echo "resolved keep_dup={params.keep_dup} (upstream remove_duplicates={REMOVE_DUPLICATES})" > {log}
            control_args=""
            if [ -n "{input.control}" ]; then
                control_args="-c {input.control}"
            fi
            nomodel_opt=""
            if [ "{params.nomodel}" = "True" ]; then
                nomodel_opt="--nomodel --extsize {params.extsize}"
            else
                nomodel_opt="--extsize {params.extsize}"
            fi
            macs3 callpeak \
                -t {input.bam} \
                $control_args \
                -f {params.fmt} \
                -g {params.gsize} \
                -n {params.name} \
                --outdir {params.outdir} \
                $nomodel_opt \
                --keep-dup {params.keep_dup} \
                -q {params.qvalue} \
                --call-summits \
                {params.broad} >> {log} 2>&1
            """

    rule macs3_pooled_peaks:
        """Call peaks on the pooled (replicate-merged) BAM with pooled input."""
        input:
            bam=lambda wc: pooled_final_bam_path(wc.sample),
            control=lambda wc: pooled_control_bam(wc.sample) or [],
            peaks=expand(macs3_narrow_path("{run_sample}"), run_sample=PEAK_RUN_SAMPLES)
        output:
            xls=macs3_xls_path("{sample}", pooled=True),
            narrow=macs3_narrow_path("{sample}", pooled=True),
            summits=macs3_summits_path("{sample}", pooled=True)
        params:
            outdir=lambda wc: f"{OUTDIR}/macs3/pooled/{wc.sample}",
            name=lambda wc: f"{wc.sample}.pooled",
            gsize=lambda wc: macs3_sample_params(wc.sample)["gsize"],
            qvalue=lambda wc: macs3_sample_params(wc.sample)["qvalue"],
            extsize=lambda wc: macs3_sample_params(wc.sample)["extsize"],
            fmt="BAM" if is_single_end() else "BAMPE",
            keep_dup=lambda wc: str(macs3_sample_params(wc.sample).get("keep_dup", "all")),
            nomodel=lambda wc: macs3_sample_params(wc.sample).get("nomodel", True),
            broad=lambda wc: "--broad --broad-cutoff " + str(macs3_sample_params(wc.sample).get("broad_cutoff", 0.1))
                    if macs3_sample_params(wc.sample).get("broad", False) else ""
        log:
            f"{LOG_DIR}/macs3/{{sample}}.pooled.log"
        threads: 2
        resources:
            mem_mb=8000,
            disk_mb=20000
        conda:
            "envs/macs3.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p {params.outdir} $(dirname {log})
            echo "resolved keep_dup={params.keep_dup} (upstream remove_duplicates={REMOVE_DUPLICATES})" > {log}
            control_args=""
            if [ -n "{input.control}" ]; then
                control_args="-c {input.control}"
            fi
            nomodel_opt=""
            if [ "{params.nomodel}" = "True" ]; then
                nomodel_opt="--nomodel --extsize {params.extsize}"
            else
                nomodel_opt="--extsize {params.extsize}"
            fi
            macs3 callpeak \
                -t {input.bam} \
                $control_args \
                -f {params.fmt} \
                -g {params.gsize} \
                -n {params.name} \
                --outdir {params.outdir} \
                $nomodel_opt \
                --keep-dup {params.keep_dup} \
                -q {params.qvalue} \
                --call-summits \
                {params.broad} >> {log} 2>&1
            """

    rule pooled_control_bam_merge:
        """Merge replicate final BAMs of a sample into the pooled BAM."""
        input:
            bams=lambda wc: [final_bam_path(run_sample_id(wc.sample, rep)) for rep in REPLICATES[wc.sample]],
            bais=lambda wc: [final_bai_path(run_sample_id(wc.sample, rep)) for rep in REPLICATES[wc.sample]]
        output:
            bam=pooled_final_bam_path("{sample}"),
            bai=pooled_final_bai_path("{sample}")
        log:
            f"{LOG_DIR}/bowtie2/{{sample}}.pooled_merge.log"
        threads: int(config["threads"]["samtools"])
        resources:
            mem_mb=8000,
            disk_mb=20000
        conda:
            "envs/bowtie2.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.bam}) $(dirname {log})
            echo "granted_threads={threads}" > {log}
            bams=({input.bams})
            M_T=$(( {threads} - 1 )); [ "$M_T" -lt 0 ] && M_T=0
            if [ "${{#bams[@]}}" -eq 1 ]; then
                cp "${{bams[0]}}" {output.bam} 2>> {log}
            else
                samtools merge -@ $M_T -f {output.bam} "${{bams[@]}}" 2>> {log}
            fi
            samtools index -@ $M_T {output.bam} 2>> {log}
            samtools quickcheck {output.bam} || {{ echo "quickcheck failed on pooled BAM" >&2; exit 1; }}
            """

    def concordance_run_samples(sample):
        return [rs for rs in sample_replicates(sample) if rs in PEAK_RUN_SAMPLES]

    rule replicate_concordance:
        """Order-independent replicate peak consensus + concordance metrics.

        Uses replicate_consensus.py (symmetric consensus, bidirectional
        overlap, bp Jaccard); NOT the previous order-dependent intersect."""
        input:
            peaks=lambda wc: [
                macs3_narrow_path(run_sample_id(wc.sample, rep))
                for rep in REPLICATES[wc.sample]
                if run_sample_id(wc.sample, rep) in PEAK_RUN_SAMPLES
            ],
            pooled=macs3_narrow_path("{sample}", pooled=True)
        output:
            tsv=replicate_concordance_path("{sample}"),
            consensus=macs3_consensus_path("{sample}")
        log:
            f"{LOG_DIR}/qc/{{sample}}.replicate_concordance.log"
        threads: 1
        conda:
            "envs/bowtie2.yaml"
        script:
            "scripts/replicate_consensus.py"
