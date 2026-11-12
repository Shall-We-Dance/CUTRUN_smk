# One shared molecule population. PE pairs are validated before downstream use.
rule prepare_fragments:
    wildcard_constraints:
        run_sample=wildcard_pattern(RUN_SAMPLES)
    input:
        bam=lambda wc: final_bam_path(wc.run_sample),
        chrom_sizes=config["reference"]["chrom_sizes"],
        blacklist=blacklist_path() if FILTER_BLACKLIST else [],
        implementation=[
            "workflow/rules/scripts/prepare_fragments.py",
            "workflow/rules/scripts/fragment_records.py",
            "workflow/rules/scripts/consensus.py",
        ]
    output:
        bed=fragments_path("{run_sample}"),
        metrics=cutrun_usable_fragments_path("{run_sample}"),
        lengths=fragment_length_path("{run_sample}")
    params:
        read_type=READ_TYPE,
        exclude_flags=FRAGMENT_EXCLUDE_FLAGS,
        remove_chrm="--remove-chrm" if FRAGMENT_REMOVE_CHRM else "",
        se_extension=lambda wc: fragment_se_extension(sample_of(wc.run_sample))
    log:
        f"{LOG_DIR}/fragments/{{run_sample}}.log"
    threads: int(config["threads"]["samtools"])
    resources:
        mem_mb=4000,
        disk_mb=20000
    conda:
        "envs/cutrun_qc.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log:q})
        blacklist_args=()
        for path in {input.blacklist:q}; do blacklist_args+=(--blacklist "$path"); done
        python workflow/rules/scripts/prepare_fragments.py \
            --bam {input.bam:q} --chrom-sizes {input.chrom_sizes:q} \
            --read-type {params.read_type} --exclude-flags {params.exclude_flags} \
            --se-extension {params.se_extension} {params.remove_chrm} "${{blacklist_args[@]}}" \
            --sample {wildcards.run_sample:q} --threads {threads} \
            --out {output.bed:q} --metrics {output.metrics:q} \
            --lengths {output.lengths:q} > {log:q} 2>&1
        """

rule pool_fragments:
    wildcard_constraints:
        sample=POOLED_SAMPLE_PATTERN
    input:
        lambda wc: [fragments_path(rs) for rs in sample_replicates(wc.sample)]
    output:
        bed=pooled_fragments_path("{sample}")
    log:
        f"{LOG_DIR}/fragments/{{sample}}.pooled.log"
    threads: 1
    resources:
        mem_mb=4000,
        disk_mb=20000
    conda:
        "envs/cutrun_qc.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bed:q}) $(dirname {log:q})
        LC_ALL=C sort -k1,1 -k2,2n -k3,3n {input:q} > {output.bed:q} 2> {log:q}
        """
