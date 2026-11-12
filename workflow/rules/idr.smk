# IDR uses separate relaxed calls, never the stringent consensus BED.
if IDR_SAMPLES:
    def idr_call_shell():
        return r"""
            set -euo pipefail
            mkdir -p {params.outdir:q} $(dirname {log:q})
            if [ ! -s {input.fragments:q} ]; then
                : > {output.peaks:q}
                echo "No fragments: no relaxed peaks." > {log:q}
                exit 0
            fi
            control_args=()
            for control in {input.control:q}; do
                test -s "$control"
                control_args+=(-c "$control")
            done
            macs3 callpeak -t {input.fragments:q} "${{control_args[@]}}" \
                -f BEDPE --keep-dup all --call-summits \
                -g {params.gsize} -p {params.pvalue} \
                -n {params.name:q} --outdir {params.outdir:q} > {log:q} 2>&1
        """

    rule idr_relaxed_peaks:
        wildcard_constraints:
            run_sample=IDR_RUN_SAMPLE_PATTERN
        input:
            fragments=lambda wc: fragments_path(wc.run_sample),
            control=macs3_control_fragments
        output:
            peaks=idr_relaxed_path("{run_sample}")
        params:
            outdir=lambda wc: os.path.dirname(idr_relaxed_prefix(wc.run_sample)),
            name="{run_sample}",
            gsize=lambda wc: macs3_config_for_sample(sample_of(wc.run_sample))["gsize"],
            pvalue=lambda wc: idr_options(sample_of(wc.run_sample)).get("relaxed_pvalue", 0.1)
        log:
            f"{LOG_DIR}/idr/{{run_sample}}.relaxed.log"
        threads: 2
        resources:
            mem_mb=8000
        conda:
            "envs/macs3.yaml"
        shell: idr_call_shell()

    rule idr_relaxed_pooled:
        wildcard_constraints:
            sample=IDR_SAMPLE_PATTERN
        input:
            fragments=lambda wc: pooled_fragments_path(wc.sample),
            control=lambda wc: pooled_control_fragments(wc.sample) or []
        output:
            peaks=idr_relaxed_path("{sample}", pooled=True)
        params:
            outdir=lambda wc: os.path.dirname(idr_relaxed_prefix(wc.sample, pooled=True)),
            name="{sample}",
            gsize=lambda wc: macs3_config_for_sample(wc.sample)["gsize"],
            pvalue=lambda wc: idr_options(wc.sample).get("relaxed_pvalue", 0.1)
        log:
            f"{LOG_DIR}/idr/{{sample}}.pooled.relaxed.log"
        threads: 2
        resources:
            mem_mb=8000
        conda:
            "envs/macs3.yaml"
        shell: idr_call_shell()

    rule replicate_idr:
        wildcard_constraints:
            sample=IDR_SAMPLE_PATTERN
        input:
            peaks=lambda wc: [idr_relaxed_path(rs) for rs in sample_replicates(wc.sample)],
            pooled=lambda wc: idr_relaxed_path(wc.sample, pooled=True),
            implementation=["workflow/rules/scripts/replicate_idr.py", "workflow/rules/scripts/consensus.py"]
        output:
            pairs=directory(idr_pairs_path("{sample}")),
            summary=idr_summary_path("{sample}")
        params:
            replicates=lambda wc: REPLICATES[wc.sample],
            threshold=lambda wc: idr_options(wc.sample).get("threshold", 0.05),
            rank=lambda wc: idr_options(wc.sample).get("rank", "p.value"),
            min_peaks=lambda wc: idr_options(wc.sample).get("min_peaks", 20)
        log:
            f"{LOG_DIR}/idr/{{sample}}.idr.log"
        threads: 1
        resources:
            mem_mb=4000
        conda:
            "envs/idr.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {log:q})
            python workflow/rules/scripts/replicate_idr.py \
                --replicates {params.replicates:q} --peaks {input.peaks:q} \
                --pooled {input.pooled:q} --outdir {output.pairs:q} \
                --summary {output.summary:q} --threshold {params.threshold} \
                --rank {params.rank} --min-peaks {params.min_peaks} > {log:q} 2>&1
            """
