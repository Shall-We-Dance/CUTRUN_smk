# rules/align_bowtie2.smk

import os

# ----------------------------------------------------------------------------
# Alignment Rules
#
# CPU budget contract: every rule computes its internal thread split from the
# *actually granted* Snakemake ``threads`` value (which may be scaled down by
# ``--cores``). A pipe of two multi-threaded tools (e.g. bowtie2 -p A |
# samtools sort -@ B) peaks at A + (B + 1) threads, so the split guarantees
# A + B + 1 <= {threads}. For tiny allocations (threads <= 2) tools run
# sequentially (SERIAL=1) so peak concurrent threads never exceed the grant.
# Multi-stage pipes (dedup) budget (threads - N) extra threads across N
# stages. Rules fail fast: no silent pass-through of tool failures.
# ----------------------------------------------------------------------------


def _bowtie2_params(wc, for_spikein=False):
    """Shared params for the bowtie2 alignment rules."""
    base = {
        "max_insert": config.get("bowtie2", {}).get("max_insert_size", 500),
        "min_insert": config.get("bowtie2", {}).get("min_insert_size", 10),
        "mode": config.get("bowtie2", {}).get("mode", "very-sensitive"),
        "k_arg": bowtie2_k_arg(),
        "extra": config.get("bowtie2", {}).get("extra", ""),
        "sort_mem": config.get("samtools_sort_memory", "2G"),
    }
    if for_spikein:
        base["index"] = config["spikein"]["bowtie2_index"]
    else:
        base["index"] = config["reference"]["bowtie2_index"]
    return base


if not is_single_end():
    rule bowtie2_pe:
        input:
            r1=fastp_sample_r1("{run_sample}"),
            r2=fastp_sample_r2("{run_sample}")
        params:
            **_bowtie2_params({}),
            pipe_split=pipe_thread_split()
        output:
            bam=maybe_temp(aligned_bam_path("{run_sample}")),
            bai=maybe_temp(aligned_bam_path("{run_sample}") + ".bai"),
            flagstat=aligned_flagstat_path("{run_sample}")
        log:
            bowtie2_log_path("{run_sample}")
        threads: int(config["threads"]["bowtie2"])
        resources:
            mem_mb=8000,
            disk_mb=50000
        benchmark:
            f"{OUTDIR}/qc/benchmarks/{ANALYSIS_MODE}/bowtie2/{{run_sample}}.align.benchmark.tsv"
        conda:
            "envs/bowtie2.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.bam}) $(dirname {output.flagstat})
            THREADS={threads}
            {params.pipe_split}
            echo "granted_threads={threads} bowtie2_p=$MAIN_T sort_at=$SORT_T serial=$SERIAL" > {log}
            tmp_sam="$(mktemp --suffix=.sam)"
            trap 'rm -f "$tmp_sam"' EXIT
            if [ "$SERIAL" = "1" ]; then
                bowtie2 \
                    -x {params.index} \
                    -1 {input.r1} \
                    -2 {input.r2} \
                    --end-to-end \
                    --{params.mode} \
                    -I {params.min_insert} \
                    -X {params.max_insert} \
                    --no-mixed \
                    --no-discordant \
                    {params.k_arg} \
                    {params.extra} \
                    --rg-id {wildcards.run_sample} --rg "SM:{wildcards.run_sample}" --rg PL:ILLUMINA \
                    -p 1 -S "$tmp_sam" 2>> {log}
                samtools sort -@ 0 -m {params.sort_mem} -o {output.bam} "$tmp_sam" 2>> {log}
            else
                bowtie2 \
                    -x {params.index} \
                    -1 {input.r1} \
                    -2 {input.r2} \
                    --end-to-end \
                    --{params.mode} \
                    -I {params.min_insert} \
                    -X {params.max_insert} \
                    --no-mixed \
                    --no-discordant \
                    {params.k_arg} \
                    {params.extra} \
                    --rg-id {wildcards.run_sample} --rg "SM:{wildcards.run_sample}" --rg PL:ILLUMINA \
                    -p $MAIN_T 2>> {log} | \
                samtools view -bS -@ 0 - | \
                samtools sort -@ $SORT_T -m {params.sort_mem} -o {output.bam} -
            fi

            echo "Running samtools index for {wildcards.run_sample}..." >> {log} 2>&1
            samtools index -@ $SORT_T {output.bam} >> {log} 2>&1
            samtools flagstat {output.bam} > {output.flagstat}
            """


else:
    rule bowtie2_se:
        input:
            r1=fastp_sample_r1("{run_sample}")
        params:
            **_bowtie2_params({}),
            pipe_split=pipe_thread_split()
        output:
            bam=maybe_temp(aligned_bam_path("{run_sample}")),
            bai=maybe_temp(aligned_bam_path("{run_sample}") + ".bai"),
            flagstat=aligned_flagstat_path("{run_sample}")
        log:
            bowtie2_log_path("{run_sample}")
        threads: int(config["threads"]["bowtie2"])
        resources:
            mem_mb=8000,
            disk_mb=50000
        benchmark:
            f"{OUTDIR}/qc/benchmarks/{ANALYSIS_MODE}/bowtie2/{{run_sample}}.align.benchmark.tsv"
        conda:
            "envs/bowtie2.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.bam}) $(dirname {output.flagstat})
            THREADS={threads}
            {params.pipe_split}
            echo "granted_threads={threads} bowtie2_p=$MAIN_T sort_at=$SORT_T serial=$SERIAL" > {log}
            tmp_sam="$(mktemp --suffix=.sam)"
            trap 'rm -f "$tmp_sam"' EXIT
            if [ "$SERIAL" = "1" ]; then
                bowtie2 \
                    -x {params.index} \
                    -U {input.r1} \
                    --end-to-end \
                    --{params.mode} \
                    {params.k_arg} \
                    {params.extra} \
                    --rg-id {wildcards.run_sample} --rg "SM:{wildcards.run_sample}" --rg PL:ILLUMINA \
                    -p 1 -S "$tmp_sam" 2>> {log}
                samtools sort -@ 0 -m {params.sort_mem} -o {output.bam} "$tmp_sam" 2>> {log}
            else
                bowtie2 \
                    -x {params.index} \
                    -U {input.r1} \
                    --end-to-end \
                    --{params.mode} \
                    {params.k_arg} \
                    {params.extra} \
                    --rg-id {wildcards.run_sample} --rg "SM:{wildcards.run_sample}" --rg PL:ILLUMINA \
                    -p $MAIN_T 2>> {log} | \
                samtools view -bS -@ 0 - | \
                samtools sort -@ $SORT_T -m {params.sort_mem} -o {output.bam} -
            fi

            echo "Running samtools index for {wildcards.run_sample}..." >> {log} 2>&1
            samtools index -@ $SORT_T {output.bam} >> {log} 2>&1
            samtools flagstat {output.bam} > {output.flagstat}
            """


# Pooled final BAMs are shared by pooled tracks and MACS3. Keep this rule
# independent of MACS3 so replicated samples can still produce pooled BigWigs
# when peak calling is disabled.
rule pooled_bam_merge:
    """Merge final BAMs for samples with at least two biological replicates."""
    wildcard_constraints:
        sample=POOLED_SAMPLE_PATTERN
    input:
        bams=lambda wc: [
            final_bam_path(run_sample_id(wc.sample, rep))
            for rep in REPLICATES[wc.sample]
        ],
        bais=lambda wc: [
            final_bai_path(run_sample_id(wc.sample, rep))
            for rep in REPLICATES[wc.sample]
        ]
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
        samtools merge -@ $M_T -f {output.bam} "${{bams[@]}}" 2>> {log}
        samtools index -@ $M_T {output.bam} 2>> {log}
        samtools quickcheck {output.bam} || {{ echo "quickcheck failed on pooled BAM" >&2; exit 1; }}
        """


# ----------------------------------------------------------------------------
# Filtering Rules
# ----------------------------------------------------------------------------

rule filter_unique_mappers:
    """
    Filter BAM to keep only uniquely mapped reads (MAPQ >= threshold).
    This removes multi-mapped reads. Input BAM is coordinate sorted.
    """
    input:
        bam=lambda wc: aligned_bam_path(wc.run_sample)
    output:
        bam=maybe_temp(unique_bam_path("{run_sample}")),
        bai=maybe_temp(unique_bam_path("{run_sample}") + ".bai"),
        flagstat=unique_flagstat_path("{run_sample}")
    params:
        mapq=config.get("bowtie2", {}).get("mapq_threshold", 30),
        exclude_flags=config.get("samtools_exclude_flags", 1804),
        sort_mem=config.get("samtools_sort_memory", "2G"),
        pipe_split=pipe_thread_split()
    log:
        f"{LOG_DIR}/bowtie2/{{run_sample}}.unique.log"
    threads: int(config["threads"]["samtools"])
    resources:
        mem_mb=8000,
        disk_mb=50000
    benchmark:
        f"{OUTDIR}/qc/benchmarks/{ANALYSIS_MODE}/bowtie2/{{run_sample}}.unique.benchmark.tsv"
    conda:
        "envs/bowtie2.yaml"
    shell:
        r"""
        set -euo pipefail
        THREADS={threads}
            {params.pipe_split}
        echo "granted_threads={threads} view_at=$MAIN_T sort_at=$SORT_T serial=$SERIAL" > {log}
        tmp_bam="$(mktemp --suffix=.bam)"
        trap 'rm -f "$tmp_bam"' EXIT
        if [ "$SERIAL" = "1" ]; then
            samtools view -b -h -q {params.mapq} -F {params.exclude_flags} \
                -@ 0 {input.bam} -o "$tmp_bam" 2>> {log}
            samtools sort -@ 0 -m {params.sort_mem} -o {output.bam} "$tmp_bam" 2>> {log}
        else
            samtools view -b -h -q {params.mapq} -F {params.exclude_flags} \
                -@ $MAIN_T {input.bam} | \
            samtools sort -@ $SORT_T -m {params.sort_mem} -o {output.bam} -
        fi

        samtools index -@ $SORT_T {output.bam}

        echo "Original BAM:" >> {log}
        samtools flagstat {input.bam} >> {log}
        echo -e "\nFiltered BAM (unique mappers):" >> {log}
        samtools flagstat {output.bam} >> {log}
        samtools flagstat {output.bam} > {output.flagstat}
        """


# ----------------------------------------------------------------------------
# Blacklist filtering (template/fragment-aware)
#
# For PE data, a template is removed if ANY of its alignments overlaps a
# blacklist region: both mates are dropped together, so no orphan mates or
# stale proper-pair flags remain. For SE data each alignment is a template.
# ----------------------------------------------------------------------------
def blacklist_shell():
    return r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {output.flagstat})
        THREADS={threads}
            {params.pipe_split}
        STEP_T=$(( THREADS - 1 )); [ "$STEP_T" -lt 0 ] && STEP_T=0
        echo "granted_threads={threads} step_at=$STEP_T serial=$SERIAL" > {log}
        echo "Input BAM:" >> {log}
        samtools flagstat {input.bam} >> {log}

        if [ ! -s "{input.blacklist}" ]; then
            echo -e "\nBlacklist BED is missing or empty; FAILING (no silent pass-through)." >> {log}
            exit 1
        fi

        tmp_bam="$(mktemp --suffix=.bam)"
        tmp_ns="$(mktemp --suffix=.bam)"
        tmp_fix="$(mktemp --suffix=.bam)"
        tmp_blbam="$(mktemp --suffix=.bam)"
        tmp_blbed="$(mktemp --suffix=.bed)"
        tmp_blhit="$(mktemp --suffix=.bed)"
        tmp_blnames="$(mktemp)"
        trap 'rm -f "$tmp_bam" "$tmp_ns" "$tmp_fix" "$tmp_blbam" "$tmp_blbed" "$tmp_blhit" "$tmp_blnames"' EXIT

        nonamecheck_opt=""
        if [ "{params.nonamecheck}" = "True" ] && \
            bedtools intersect -h 2>&1 | grep -q -- '--nonamecheck'; then
            nonamecheck_opt="--nonamecheck"
        fi

        if [ "{params.read_type}" = "PE" ]; then
            # 1) name-sort and add MC tags (needed to keep pairs consistent)
            samtools sort -n -@ $STEP_T -m {params.sort_mem} -o "$tmp_ns" {input.bam} 2>> {log}
            samtools fixmate -m -@ $STEP_T -O bam "$tmp_ns" "$tmp_fix" 2>> {log}

            # 2) collect template names whose ANY alignment overlaps blacklist.
            # Staging these commands avoids a 5-process pipe exceeding the
            # Snakemake CPU grant at low core counts.
            samtools view -b -F {params.exclude_flags} -@ $STEP_T \
                -o "$tmp_blbam" "$tmp_fix" 2>> {log}
            bedtools bamtobed -i "$tmp_blbam" > "$tmp_blbed" 2>> {log}
            bedtools intersect -a "$tmp_blbed" -b {input.blacklist} \
                $nonamecheck_opt -wa > "$tmp_blhit" 2>> {log}
            awk '{{ seen[$4] = 1 }} END {{ for (q in seen) print q }}' \
                "$tmp_blhit" > "$tmp_blnames"

            # 3) drop entire templates, re-sort to coordinate order
            samtools view -b -N ^"$tmp_blnames" -@ $STEP_T \
                -o "$tmp_bam" "$tmp_fix" 2>> {log}
        else
            bedtools intersect -v -abam {input.bam} -b {input.blacklist} \
                $nonamecheck_opt > "$tmp_bam" 2>> {log}
        fi

        samtools sort -@ $STEP_T -m {params.sort_mem} -o {output.bam} "$tmp_bam" 2>> {log}
        samtools index -@ $STEP_T {output.bam} 2>> {log}

        echo -e "\nAfter blacklist filtering:" >> {log}
        samtools flagstat {output.bam} >> {log}
        samtools flagstat {output.bam} > {output.flagstat}
        """


if FILTER_BLACKLIST:
    rule filter_blacklist_bam:
        """
        Remove templates overlapping blacklist regions (fragment-aware for PE).
        """
        input:
            bam=unique_bam_path("{run_sample}"),
            blacklist=blacklist_path()
        output:
            bam=maybe_temp(filtered_bam_path("{run_sample}")),
            bai=maybe_temp(filtered_bam_path("{run_sample}") + ".bai"),
            flagstat=filtered_flagstat_path("{run_sample}")
        params:
            nonamecheck=config.get("bedtools_nonamecheck", True),
            exclude_flags=config.get("samtools_exclude_flags", 1804),
            read_type=READ_TYPE,
            sort_mem=config.get("samtools_sort_memory", "2G"),
            pipe_split=pipe_thread_split()
        threads: int(config["threads"]["samtools"])
        resources:
            mem_mb=8000,
            disk_mb=50000
        benchmark:
            f"{OUTDIR}/qc/benchmarks/{ANALYSIS_MODE}/bowtie2/{{run_sample}}.blacklist.benchmark.tsv"
        log:
            f"{LOG_DIR}/bowtie2/{{run_sample}}.filter_blacklist.log"
        conda:
            "envs/bowtie2.yaml"
        shell: blacklist_shell()


# ----------------------------------------------------------------------------
# Duplicate removal
#
# dedup_method = picard  : `picard MarkDuplicates` with the PicardCommandLine
#                          wrapper syntax (no `--java-options` subcommand;
#                          JVM memory via JAVA_TOOL_OPTIONS). Real Picard
#                          metrics are written to the metrics output.
# dedup_method = samtools : collate/name-sort -> fixmate -m -> coordinate sort
#                          -> markdup -r, with markdup stats as metrics.
# Both modes fail the job on error (no pass-through of the input BAM).
# ----------------------------------------------------------------------------
if REMOVE_DUPLICATES:
    rule remove_duplicates:
        input:
            bam=lambda wc: filtered_bam_path(wc.run_sample) if FILTER_BLACKLIST else unique_bam_path(wc.run_sample)
        output:
            bam=maybe_temp(dedup_bam_path("{run_sample}")),
            bai=maybe_temp(dedup_bam_path("{run_sample}") + ".bai"),
            metrics=duplicate_metrics_path("{run_sample}"),
            flagstat=dedup_flagstat_path("{run_sample}")
        params:
            method=config.get("dedup_method", "picard"),
            picard_mem=config.get("picard_java_mem", "4g"),
            sort_mem=config.get("samtools_sort_memory", "2G")
        threads: int(config["threads"].get("dedup", config["threads"].get("samtools", 8)))
        resources:
            mem_mb=8000,
            disk_mb=50000
        benchmark:
            f"{OUTDIR}/qc/benchmarks/{ANALYSIS_MODE}/bowtie2/{{run_sample}}.dedup.benchmark.tsv"
        log:
            f"{LOG_DIR}/bowtie2/{{run_sample}}.dedup.log"
        conda:
            "envs/bowtie2.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.bam}) $(dirname {output.flagstat}) $(dirname {output.metrics})
            echo "granted_threads={threads}" > {log}

            if [ "{params.method}" == "picard" ]; then
                command -v picard >/dev/null 2>&1 || {{ echo "picard not found in PATH (dedup_method=picard)" >&2; exit 1; }}
                picard_tmpdir="$(mktemp -d)"
                trap 'rm -rf "$picard_tmpdir"' EXIT
                export JAVA_TOOL_OPTIONS="-Xmx{params.picard_mem}"
                picard MarkDuplicates \
                    I={input.bam} \
                    O={output.bam} \
                    M={output.metrics} \
                    TMP_DIR="$picard_tmpdir" \
                    REMOVE_DUPLICATES=true \
                    ASSUME_SORTED=true \
                    VALIDATION_STRINGENCY=LENIENT \
                    2>> {log}
                samtools index -@ 0 {output.bam} 2>> {log}
            else
                # 4-stage pipeline: budget (threads - 4) extra threads, one
                # per stage, so peak concurrent threads <= granted threads.
                if [ {threads} -le 4 ]; then
                    D_T=0
                else
                    D_T=$(( ({threads} - 4) / 4 ))
                fi
                tmp_md="$(mktemp --suffix=.bam)"
                trap 'rm -f "$tmp_md"' EXIT
                if [ {threads} -le 4 ]; then
                    tmp_c="$(mktemp --suffix=.bam)"
                    tmp_f="$(mktemp --suffix=.bam)"
                    trap 'rm -f "$tmp_md" "$tmp_c" "$tmp_f"' EXIT
                    samtools collate -@ 0 -u -o "$tmp_c" {input.bam} 2>> {log}
                    samtools fixmate -m -@ 0 -O bam -o "$tmp_f" "$tmp_c" 2>> {log}
                    samtools sort -@ 0 -m {params.sort_mem} -o "$tmp_md" "$tmp_f" 2>> {log}
                else
                    samtools collate -@ $D_T -O -u {input.bam} 2>> {log} | \
                    samtools fixmate -m -@ $D_T - - 2>> {log} | \
                    samtools sort -@ $D_T -m {params.sort_mem} -o /dev/stdout - 2>> {log} | \
                    samtools markdup -r -@ $D_T -f {output.metrics} - "$tmp_md" 2>> {log}
                fi
                if [ {threads} -le 4 ]; then
                    samtools markdup -r -@ 0 -f {output.metrics} "$tmp_md" {output.bam} 2>> {log}
                    rm -f "$tmp_c" "$tmp_f"
                else
                    mv "$tmp_md" {output.bam}
                fi
                samtools index -@ 0 {output.bam} 2>> {log}
            fi

            samtools quickcheck {output.bam} || {{ echo "samtools quickcheck failed on {output.bam}" >&2; exit 1; }}
            samtools flagstat {output.bam} > {output.flagstat}
            echo -e "\nAfter deduplication:" >> {log}
            cat {output.flagstat} >> {log}
            """


# ----------------------------------------------------------------------------
# TE-mode all-alignment BAM (multimappers retained, secondary kept)
# ----------------------------------------------------------------------------
if TE_ENABLED:
    rule filter_te_multimappers:
        """
        Build the TE-mode BAM that RETAINS multi-mapped alignments
        (bowtie2 -k output). Only unmapped/QC-failed reads are dropped;
        secondary alignments (flag 256) are kept so repetitive/TE loci
        contribute signal. Blacklist (te.blacklist=full) and duplicate
        removal (te.remove_duplicates) are applied when configured.
        """
        input:
            bam=lambda wc: aligned_bam_path(wc.run_sample),
            blacklist=te_filter_blacklist_input()
        output:
            bam=te_bam_path("{run_sample}"),
            bai=te_bai_path("{run_sample}"),
            flagstat=te_flagstat_path("{run_sample}")
        params:
            exclude_flags=TE_EXCLUDE_FLAGS,
            blacklist_mode=TE_BLACKLIST_MODE,
            remove_duplicates=TE_REMOVE_DUPLICATES,
            nonamecheck=config.get("bedtools_nonamecheck", True),
            read_type=READ_TYPE,
            sort_mem=config.get("samtools_sort_memory", "2G"),
            pipe_split=pipe_thread_split()
        log:
            f"{LOG_DIR}/bowtie2/{{run_sample}}.te_all.log"
        threads: int(config["threads"]["samtools"])
        resources:
            mem_mb=8000,
            disk_mb=50000
        benchmark:
            f"{OUTDIR}/qc/benchmarks/{ANALYSIS_MODE}/bowtie2/{{run_sample}}.te_all.benchmark.tsv"
        conda:
            "envs/bowtie2.yaml"
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.bam}) $(dirname {output.flagstat})
            THREADS={threads}
            {params.pipe_split}
            STEP_T=$(( THREADS - 1 )); [ "$STEP_T" -lt 0 ] && STEP_T=0
            echo "granted_threads={threads} step_at=$STEP_T serial=$SERIAL" > {log}

            tmp_bam="$(mktemp --suffix=.bam)"
            tmp_ns="$(mktemp --suffix=.bam)"
            tmp_fix="$(mktemp --suffix=.bam)"
            tmp_blbam="$(mktemp --suffix=.bam)"
            tmp_blbed="$(mktemp --suffix=.bed)"
            tmp_blhit="$(mktemp --suffix=.bed)"
            tmp_blnames="$(mktemp)"
            tmp_md="$(mktemp --suffix=.bam)"
            trap 'rm -f "$tmp_bam" "$tmp_ns" "$tmp_fix" "$tmp_blbam" "$tmp_blbed" "$tmp_blhit" "$tmp_blnames" "$tmp_md"' EXIT

            # keep all mapped alignments incl. secondary (multimappers)
            if [ "$SERIAL" = "1" ]; then
                samtools view -b -F {params.exclude_flags} -@ 0 \
                    {input.bam} -o "$tmp_bam" 2>> {log}
            else
                samtools view -b -F {params.exclude_flags} -@ $MAIN_T \
                    {input.bam} | \
                samtools sort -@ $SORT_T -m {params.sort_mem} -o "$tmp_bam" - 2>> {log}
            fi

            if [ "{params.blacklist_mode}" = "full" ]; then
                if [ ! -s "{input.blacklist}" ]; then
                    echo "te.blacklist=full but blacklist BED missing; FAILING." >> {log}
                    exit 1
                fi
                nonamecheck_opt=""
                if [ "{params.nonamecheck}" = "True" ] && \
                    bedtools intersect -h 2>&1 | grep -q -- '--nonamecheck'; then
                    nonamecheck_opt="--nonamecheck"
                fi
                if [ "{params.read_type}" = "PE" ]; then
                    samtools sort -n -@ $STEP_T -m {params.sort_mem} -o "$tmp_ns" "$tmp_bam" 2>> {log}
                    samtools fixmate -m -@ $STEP_T -O bam "$tmp_ns" "$tmp_fix" 2>> {log}
                    samtools view -b -F {params.exclude_flags} -@ $STEP_T \
                        -o "$tmp_blbam" "$tmp_fix" 2>> {log}
                    bedtools bamtobed -i "$tmp_blbam" > "$tmp_blbed" 2>> {log}
                    bedtools intersect -a "$tmp_blbed" -b {input.blacklist} \
                        $nonamecheck_opt -wa > "$tmp_blhit" 2>> {log}
                    awk '{{ seen[$4] = 1 }} END {{ for (q in seen) print q }}' \
                        "$tmp_blhit" > "$tmp_blnames"
                    samtools view -b -N ^"$tmp_blnames" -@ $STEP_T \
                        -o "$tmp_md" "$tmp_fix" 2>> {log}
                    samtools sort -@ $STEP_T -m {params.sort_mem} \
                        -o "$tmp_bam" "$tmp_md" 2>> {log}
                else
                    bedtools intersect -v -abam "$tmp_bam" -b {input.blacklist} \
                        $nonamecheck_opt > "$tmp_md" 2>> {log}
                    mv "$tmp_md" "$tmp_bam"
                fi
            fi

            if [ "{params.remove_duplicates}" = "True" ]; then
                echo "TE BAM: removing duplicates (samtools markdup pipeline)" >> {log}
                if [ {threads} -le 4 ]; then
                    tmp_c="$(mktemp --suffix=.bam)"
                    tmp_f="$(mktemp --suffix=.bam)"
                    trap 'rm -f "$tmp_bam" "$tmp_ns" "$tmp_fix" "$tmp_blnames" "$tmp_md" "$tmp_c" "$tmp_f"' EXIT
                    samtools collate -@ 0 -u -o "$tmp_c" "$tmp_bam" 2>> {log}
                    samtools fixmate -m -@ 0 -O bam -o "$tmp_f" "$tmp_c" 2>> {log}
                    samtools sort -@ 0 -m {params.sort_mem} -o "$tmp_md" "$tmp_f" 2>> {log}
                    samtools markdup -r -@ 0 "$tmp_md" "$tmp_md.out" 2>> {log}
                    mv "$tmp_md.out" "$tmp_bam"
                    rm -f "$tmp_c" "$tmp_f"
                else
                    D_T=$(( ({threads} - 4) / 4 ))
                    [ "$D_T" -lt 1 ] && D_T=1
                    if [ "$SERIAL" = "1" ]; then
                        :
                    fi
                    samtools collate -@ $D_T -O -u "$tmp_bam" 2>> {log} | \
                    samtools fixmate -m -@ $D_T - - 2>> {log} | \
                    samtools sort -@ $D_T -m {params.sort_mem} -o /dev/stdout - 2>> {log} | \
                    samtools markdup -r -@ $D_T - "$tmp_md" 2>> {log}
                    mv "$tmp_md" "$tmp_bam"
                fi
            fi

            mv "$tmp_bam" {output.bam}
            samtools index -@ $SORT_T {output.bam} 2>> {log}

            echo -e "\nTE BAM (all alignments incl. multimappers):" >> {log}
            samtools flagstat {output.bam} >> {log}
            samtools flagstat {output.bam} > {output.flagstat}
            """
