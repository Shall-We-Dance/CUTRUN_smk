# rules/macs3.smk

if MACS3_ENABLED:
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
        merged = macs3_config_for_sample(sample)
        if "keep_dup" not in merged:
            merged["keep_dup"] = "all" if REMOVE_DUPLICATES else "1"
        return merged

    def macs3_control_bam(wc, pooled=False):
        sample = sample_of(wc.run_sample) if not pooled else wc.sample
        control = config["samples"].get(sample, {}).get("control")
        if not control:
            return []
        if pooled:
            return pooled_control_bam(sample)
        rep = replicate_of(wc.run_sample)
        control_rep = control_replicate_for(sample, rep)
        if control_rep == "pooled":
            return pooled_control_bam(sample)
        return final_bam_path(run_sample_id(control, control_rep))

    def macs3_callpeak_shell():
        """Shared MACS3 command; mode_args is mutually exclusive by rule."""
        return r"""
            set -euo pipefail
            mkdir -p {params.outdir} $(dirname {log})
            echo "resolved keep_dup={params.keep_dup} (upstream remove_duplicates={REMOVE_DUPLICATES})" > {log}
            echo "peak_mode={params.peak_mode}" >> {log}
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
                {params.mode_args} >> {log} 2>&1
            """

    rule macs3_callpeak:
        """Call narrow peaks and summits for one replicate."""
        wildcard_constraints:
            run_sample=NARROW_PEAK_RUN_SAMPLE_PATTERN
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
            peak_mode="narrow",
            mode_args="--call-summits"
        log:
            f"{LOG_DIR}/macs3/{{run_sample}}.log"
        threads: 2
        resources:
            mem_mb=8000,
            disk_mb=20000
        conda:
            "envs/macs3.yaml"
        shell: macs3_callpeak_shell()

    rule macs3_callpeak_broad:
        """Call broad and gapped peaks for one replicate; MACS3 broad mode has no summits."""
        wildcard_constraints:
            run_sample=BROAD_PEAK_RUN_SAMPLE_PATTERN
        input:
            bam=lambda wc: final_bam_path(wc.run_sample),
            control=macs3_control_bam
        output:
            xls=macs3_xls_path("{run_sample}"),
            broad=macs3_broad_path("{run_sample}"),
            gapped=macs3_gapped_path("{run_sample}")
        params:
            outdir=lambda wc: f"{OUTDIR}/macs3/{wc.run_sample}",
            name="{run_sample}",
            gsize=lambda wc: macs3_sample_params(sample_of(wc.run_sample))["gsize"],
            qvalue=lambda wc: macs3_sample_params(sample_of(wc.run_sample))["qvalue"],
            extsize=lambda wc: macs3_sample_params(sample_of(wc.run_sample))["extsize"],
            fmt="BAM" if is_single_end() else "BAMPE",
            keep_dup=lambda wc: str(macs3_sample_params(sample_of(wc.run_sample)).get("keep_dup", "all")),
            nomodel=lambda wc: macs3_sample_params(sample_of(wc.run_sample)).get("nomodel", True),
            peak_mode="broad",
            mode_args=lambda wc: "--broad --broad-cutoff " + str(
                macs3_sample_params(sample_of(wc.run_sample)).get("broad_cutoff", 0.1)
            )
        log:
            f"{LOG_DIR}/macs3/{{run_sample}}.log"
        threads: 2
        resources:
            mem_mb=8000,
            disk_mb=20000
        conda:
            "envs/macs3.yaml"
        shell: macs3_callpeak_shell()

    rule macs3_pooled_peaks:
        """Call narrow peaks and summits on a pooled replicate BAM."""
        wildcard_constraints:
            sample=POOLED_NARROW_PEAK_SAMPLE_PATTERN
        input:
            bam=lambda wc: pooled_final_bam_path(wc.sample),
            control=lambda wc: pooled_control_bam(wc.sample) or []
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
            peak_mode="narrow",
            mode_args="--call-summits"
        log:
            f"{LOG_DIR}/macs3/{{sample}}.pooled.log"
        threads: 2
        resources:
            mem_mb=8000,
            disk_mb=20000
        conda:
            "envs/macs3.yaml"
        shell: macs3_callpeak_shell()

    rule macs3_pooled_peaks_broad:
        """Call broad and gapped peaks on a pooled replicate BAM; no summits."""
        wildcard_constraints:
            sample=POOLED_BROAD_PEAK_SAMPLE_PATTERN
        input:
            bam=lambda wc: pooled_final_bam_path(wc.sample),
            control=lambda wc: pooled_control_bam(wc.sample) or []
        output:
            xls=macs3_xls_path("{sample}", pooled=True),
            broad=macs3_broad_path("{sample}", pooled=True),
            gapped=macs3_gapped_path("{sample}", pooled=True)
        params:
            outdir=lambda wc: f"{OUTDIR}/macs3/pooled/{wc.sample}",
            name=lambda wc: f"{wc.sample}.pooled",
            gsize=lambda wc: macs3_sample_params(wc.sample)["gsize"],
            qvalue=lambda wc: macs3_sample_params(wc.sample)["qvalue"],
            extsize=lambda wc: macs3_sample_params(wc.sample)["extsize"],
            fmt="BAM" if is_single_end() else "BAMPE",
            keep_dup=lambda wc: str(macs3_sample_params(wc.sample).get("keep_dup", "all")),
            nomodel=lambda wc: macs3_sample_params(wc.sample).get("nomodel", True),
            peak_mode="broad",
            mode_args=lambda wc: "--broad --broad-cutoff " + str(
                macs3_sample_params(wc.sample).get("broad_cutoff", 0.1)
            )
        log:
            f"{LOG_DIR}/macs3/{{sample}}.pooled.log"
        threads: 2
        resources:
            mem_mb=8000,
            disk_mb=20000
        conda:
            "envs/macs3.yaml"
        shell: macs3_callpeak_shell()

    def replicate_peak_inputs(wc):
        return [
            macs3_peak_path(run_sample_id(wc.sample, rep))
            for rep in REPLICATES[wc.sample]
            if run_sample_id(wc.sample, rep) in PEAK_RUN_SAMPLES
        ]

    rule replicate_concordance:
        """Order-independent replicate peak consensus + concordance metrics.

        Uses replicate_consensus.py (symmetric consensus, bidirectional
        overlap, bp Jaccard); NOT the previous order-dependent intersect."""
        wildcard_constraints:
            sample=POOLED_NARROW_PEAK_SAMPLE_PATTERN
        input:
            peaks=replicate_peak_inputs
        output:
            tsv=replicate_concordance_path("{sample}"),
            consensus=f"{OUTDIR}/macs3/consensus/{{sample}}.consensus.narrowPeak"
        log:
            f"{LOG_DIR}/qc/{{sample}}.replicate_concordance.log"
        threads: 1
        conda:
            "envs/bowtie2.yaml"
        script:
            "scripts/replicate_consensus.py"

    rule replicate_concordance_broad:
        """Replicate concordance for broadPeak inputs, with a generic BED consensus."""
        wildcard_constraints:
            sample=POOLED_BROAD_PEAK_SAMPLE_PATTERN
        input:
            peaks=replicate_peak_inputs
        output:
            tsv=replicate_concordance_path("{sample}"),
            consensus=macs3_broad_consensus_path("{sample}")
        log:
            f"{LOG_DIR}/qc/{{sample}}.replicate_concordance.log"
        threads: 1
        conda:
            "envs/bowtie2.yaml"
        script:
            "scripts/replicate_consensus.py"
