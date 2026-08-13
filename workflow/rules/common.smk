# workflow/rules/common.smk
import os
import re
from snakemake.io import temp

OUTDIR = config.get("output", {}).get("dir", "results")
ANALYSIS_MODE = str(config.get("analysis_mode", "standard")).lower()
LOG_DIR = f"{OUTDIR}/logs"

SAMPLES = list(config.get("samples", {}).keys())
READ_TYPE = config.get("read_type", "PE").upper()


def is_single_end():
    return READ_TYPE == "SE"


FILTER_BLACKLIST = bool(config.get("filter_blacklist", False))
REMOVE_DUPLICATES = bool(config.get("remove_duplicates", True))
KEEP_BAM = bool(config.get("keep_bam", False))
CUTRUN_QC_ENABLED = bool(config.get("cutrun_qc", {}).get("enabled", True))
SPIKEIN_ENABLED = bool(config.get("spikein", {}).get("enabled", False))
SPIKEIN_BIGWIG_ENABLED = (
    SPIKEIN_ENABLED and bool(config.get("spikein", {}).get("generate_bigwig", True))
)
MACS3_CONFIG = config.get("macs3", {})
MACS3_ENABLED = bool(MACS3_CONFIG.get("enabled", True))

GENERATE_RAW_BIGWIG = bool(config.get("generate_raw_bigwig", True))
GENERATE_NORMALIZED_BIGWIG = bool(config.get("generate_normalized_bigwig", True))

# ----------------------------------------------------------------------------
# Analysis mode: "standard" (default) or "te" (transposable-element friendly)
# ----------------------------------------------------------------------------
TE_CONFIG = config.get("te", {})
TE_ENABLED = ANALYSIS_MODE == "te"
TE_ANNOTATION = TE_CONFIG.get("annotation") or None
TE_EXCLUDE_FLAGS = int(TE_CONFIG.get("exclude_flags", 516))
TE_BLACKLIST_MODE = str(TE_CONFIG.get("blacklist", "none")).lower()
TE_REMOVE_DUPLICATES = bool(TE_CONFIG.get("remove_duplicates", False))
TE_GENERATE_BIGWIG = bool(TE_CONFIG.get("generate_bigwig", True))
# TE-specific k takes precedence. Keeping this in the TE section prevents a
# standard-mode bowtie2 tuning value from silently changing TE sensitivity.
TE_K = int(TE_CONFIG.get("k") or config.get("bowtie2", {}).get("k", 1000))
TE_NORMALIZATION = TE_CONFIG.get("normalization", "RPKM")
TE_BIN_SIZE = int(TE_CONFIG.get("bin_size", 10))
TE_FAMILY_COLUMN = int(TE_CONFIG.get("family_column", 9))
TE_CLASS_COLUMN = int(TE_CONFIG.get("class_column", 7))
TE_SUBFAMILY_COLUMN = int(TE_CONFIG.get("subfamily_column", 8))
TE_MIN_OVERLAP_FRAC = float(TE_CONFIG.get("min_overlap_frac", 0.0))

# ----------------------------------------------------------------------------
# Samples / replicates
# ----------------------------------------------------------------------------
# A "condition" (e.g. Factor_IP_Treatment) may contain several biological
# replicates. Each (condition, replicate) pair is processed independently
# through QC/alignment/dedup/peaks; pooled targets merge replicates.
REPLICATES = {
    sample: list(config["samples"][sample].get("replicates", {}).keys())
    for sample in SAMPLES
}
RUN_SAMPLES = [f"{sample}.{rep}" for sample in SAMPLES for rep in REPLICATES[sample]]


def control_replicate_for(sample, chip_rep):
    """Resolve the control replicate for one ChIP replicate.

    Preferred syntax is ``control_replicates: {chip_rep: control_rep}``.
    ``control_replicate`` is retained for one-replicate experiments or for
    mapping every ChIP replicate to ``pooled``.
    """
    sample_cfg = config["samples"][sample]
    mapping = sample_cfg.get("control_replicates", {}) or {}
    if chip_rep in mapping:
        return str(mapping[chip_rep])
    legacy = sample_cfg.get("control_replicate")
    if legacy is not None:
        return str(legacy)
    return chip_rep


def validate_sample_design():
    """Fail before DAG construction on ambiguous sample/control designs."""
    errors = []
    for sample in SAMPLES:
        if "." in sample:
            errors.append(f"sample name must not contain '.': {sample}")
        sample_cfg = config["samples"][sample]
        reps = REPLICATES[sample]
        for rep in reps:
            if "." in rep:
                errors.append(f"replicate name must not contain '.': {sample}/{rep}")
            rep_cfg = sample_cfg["replicates"][rep]
            r1 = rep_cfg.get("R1", [])
            r2 = rep_cfg.get("R2", [])
            if READ_TYPE == "PE" and not r2:
                errors.append(f"PE sample is missing R2: {sample}/{rep}")
            if READ_TYPE == "PE" and len(r1) != len(r2):
                errors.append(
                    f"R1/R2 lane counts differ for {sample}/{rep}: "
                    f"{len(r1)} vs {len(r2)}"
                )

        control = sample_cfg.get("control")
        if not control:
            continue
        if control not in config["samples"]:
            errors.append(f"unknown control '{control}' for sample '{sample}'")
            continue
        control_reps = REPLICATES[control]
        for rep in reps:
            target = control_replicate_for(sample, rep)
            if target != "pooled" and target not in control_reps:
                errors.append(
                    f"control mapping {sample}/{rep} -> {control}/{target} does not exist"
                )
    if errors:
        raise ValueError("Invalid sample/control design:\n- " + "\n- ".join(errors))


validate_sample_design()


def run_sample_parts(run_sample):
    sample, rep = run_sample.rsplit(".", 1)
    return sample, rep


def sample_of(run_sample):
    return run_sample_parts(run_sample)[0]


def replicate_of(run_sample):
    return run_sample_parts(run_sample)[1]


def sample_replicates(sample):
    return [f"{sample}.{rep}" for rep in REPLICATES[sample]]


def pooled_run_sample(sample):
    return f"{sample}.pooled"


def split_threads(total, main_frac=0.66):
    """Split a rule-level thread budget between two piped tools.

    Accounts for the main thread of each tool: for samtools-style tools the
    ``-@ N`` flag starts N *additional* threads on top of the main thread, so
    a pipe using ``-p A`` and ``-@ B`` peaks at A + B + 1 threads.
    """
    if total <= 1:
        return 1, 0
    main = max(1, int(round(total * main_frac)))
    secondary = total - main - 1  # reserve main thread of the second tool
    if secondary < 0:
        main = max(1, total - 1)
        secondary = total - main - 1
    return main, max(0, secondary)


def pipe_thread_split():
    """Shell snippet computing thread budgets from the granted threads
    (shell var THREADS) so that the peak concurrent thread count stays <=
    granted threads.

    Accounting: samtools-style ``-@ N`` means N *additional* threads, so a
    tool given ``-@ N`` really uses N + 1 threads.

    bowtie2 alignment pipe is THREE tools: bowtie2 -p MAIN | samtools view
    | samtools sort -@ SORT. samtools view runs with its single main thread
    (-@ 0), so peak = MAIN + SORT + 2 <= THREADS.

    samtools-only pipes (view | sort) are TWO tools:
    peak = VIEW + SORT + 2 <= THREADS with VIEW=SORT=0 at low cores.

    At THREADS <= 3 the rules fall back to SERIAL mode (one tool at a time,
    each single-threaded).
    """
    return r"""
        SERIAL=0
        if [ "$THREADS" -ge 4 ]; then
            MAIN_T=$(( (THREADS * 2) / 3 ))
            SORT_T=$(( THREADS - MAIN_T - 2 ))
            [ "$SORT_T" -lt 1 ] && SORT_T=1
            # rebalance to keep MAIN + SORT + 2 <= THREADS
            MAIN_T=$(( THREADS - SORT_T - 2 ))
        elif [ "$THREADS" -eq 3 ]; then
            MAIN_T=1; SORT_T=0
        else
            MAIN_T=1; SORT_T=0; SERIAL=1
        fi
    """


def dynamic_threads_snippet(tool_main="ALIGN", tool_sec="SORT"):
    """Return a shell snippet that splits the granted ``threads`` between two
    concurrently running tools inside one rule shell command."""
    return f"""
        {tool_main}=$(( {{threads}} * 2 / 3 ))
        [ "${{{tool_main}}}" -lt 1 ] && {tool_main}=1
        {tool_sec}=$(( {{threads}} - {tool_main} - 1 ))
        [ "${{{tool_sec}}}" -lt 1 ] && {tool_sec}=1
    """


def keep_multimappers():
    """Return True if multimapped alignments must be retained (TE mode)."""
    return TE_ENABLED or bool(config.get("bowtie2", {}).get("keep_multimappers", False))


def bowtie2_k_arg():
    """Return the ``-k`` argument for multimapper reporting, or empty string."""
    if not keep_multimappers():
        return ""
    return f"-k {TE_K}"


def blacklist_config():
    return config.get("blacklist", {})


def blacklist_path():
    blacklist = blacklist_config()
    if blacklist.get("path"):
        return blacklist["path"]
    if blacklist.get("url"):
        cache_dir = blacklist.get("cache_dir", "resources/blacklist")
        filename = os.path.basename(blacklist["url"])
        if filename.endswith(".gz"):
            filename = filename[:-3]
        return os.path.join(cache_dir, filename)
    return None


# --- per-replicate BAM paths --------------------------------------------
def aligned_bam_path(run_sample):
    suffix = "se" if is_single_end() else "pe"
    return f"{OUTDIR}/bowtie2/{run_sample}/{run_sample}_{suffix}.bam"


def unique_bam_path(run_sample):
    return f"{OUTDIR}/bowtie2/{run_sample}/{run_sample}.unique.bam"


def filtered_bam_path(run_sample):
    return f"{OUTDIR}/bowtie2/{run_sample}/{run_sample}.unique.filtered.bam"


def dedup_bam_path(run_sample):
    suffix = "unique.filtered" if FILTER_BLACKLIST else "unique"
    return f"{OUTDIR}/bowtie2/{run_sample}/{run_sample}.{suffix}.dedup.bam"


def duplicate_metrics_path(run_sample):
    return f"{OUTDIR}/bowtie2/{run_sample}/{run_sample}.dedup_metrics.txt"


def final_bam_path(run_sample):
    if REMOVE_DUPLICATES:
        return dedup_bam_path(run_sample)
    if FILTER_BLACKLIST:
        return filtered_bam_path(run_sample)
    return unique_bam_path(run_sample)


def final_bai_path(run_sample):
    return final_bam_path(run_sample) + ".bai"


def pooled_final_bam_path(sample):
    return f"{OUTDIR}/bowtie2/pooled/{sample}/{sample}.pooled.dedup.bam"


def pooled_final_bai_path(sample):
    return pooled_final_bam_path(sample) + ".bai"


# --- TE-mode paths --------------------------------------------------------
def te_bam_path(run_sample):
    return f"{OUTDIR}/bowtie2/{run_sample}/{run_sample}.all.bam"


def te_bai_path(run_sample):
    return te_bam_path(run_sample) + ".bai"


def te_flagstat_path(run_sample):
    return f"{OUTDIR}/qc/bam_stats/{run_sample}/{run_sample}.all.flagstat.txt"


def te_diagnostic_bigwig_path(run_sample):
    return f"{OUTDIR}/bigwig/{run_sample}/{run_sample}.te_diagnostic.bw"


def te_overlap_path(run_sample):
    return f"{OUTDIR}/te/{run_sample}/{run_sample}.te_overlap.tsv"


def te_family_counts_path(run_sample):
    return te_overlap_path(run_sample).replace(".te_overlap.tsv", ".te_family_counts.tsv")


def te_overlap_summary_path():
    return f"{OUTDIR}/te/te_overlap_summary.tsv"


def te_enrichment_path():
    return f"{OUTDIR}/te/te_family_enrichment.tsv"


def te_overlap_mqc_path():
    return f"{OUTDIR}/te/te_overlap_mqc.json"


def te_blacklist_lowconf_path():
    """Blacklist intervals overlapping TE annotation: low-confidence loci that
    must NOT be automatically rescued; reported for manual/sensitivity use."""
    return os.path.join(
        blacklist_config().get("cache_dir", "resources/blacklist"),
        "blacklist.te_overlap_lowconf.bed",
    )


def te_filter_blacklist_input():
    """Blacklist file used for the TE-mode BAM, per te.blacklist mode."""
    if TE_BLACKLIST_MODE == "full":
        return blacklist_path()
    return "/dev/null"


def te_paths(run_sample):
    if not TE_ENABLED:
        return []
    paths = [te_bam_path(run_sample), te_bai_path(run_sample), te_flagstat_path(run_sample)]
    if TE_GENERATE_BIGWIG:
        paths.append(te_diagnostic_bigwig_path(run_sample))
    if TE_ANNOTATION:
        paths.append(te_overlap_path(run_sample))
        paths.append(te_family_counts_path(run_sample))
    return paths


def maybe_temp(path):
    return path if KEEP_BAM else temp(path)


# --- QC paths ------------------------------------------------------------
def aligned_flagstat_path(run_sample):
    suffix = "se" if is_single_end() else "pe"
    return f"{OUTDIR}/qc/bam_stats/{run_sample}/{run_sample}_{suffix}.flagstat.txt"


def unique_flagstat_path(run_sample):
    return f"{OUTDIR}/qc/bam_stats/{run_sample}/{run_sample}.unique.flagstat.txt"


def filtered_flagstat_path(run_sample):
    return f"{OUTDIR}/qc/bam_stats/{run_sample}/{run_sample}.unique.filtered.flagstat.txt"


def dedup_flagstat_path(run_sample):
    suffix = "unique.filtered" if FILTER_BLACKLIST else "unique"
    return f"{OUTDIR}/qc/bam_stats/{run_sample}/{run_sample}.{suffix}.dedup.flagstat.txt"


def bam_flagstat_paths(run_sample):
    paths = [aligned_flagstat_path(run_sample), unique_flagstat_path(run_sample)]
    if FILTER_BLACKLIST:
        paths.append(filtered_flagstat_path(run_sample))
    if REMOVE_DUPLICATES:
        paths.append(dedup_flagstat_path(run_sample))
    return paths


def bowtie2_log_path(run_sample):
    suffix = "se" if is_single_end() else "pe"
    return f"{LOG_DIR}/bowtie2/{run_sample}_{suffix}.log"


def cutrun_usable_fragments_path(run_sample):
    return f"{OUTDIR}/qc/cutrun/{run_sample}/{run_sample}.usable_fragments.tsv"


def fragment_length_path(run_sample):
    return f"{OUTDIR}/qc/cutrun/{run_sample}/{run_sample}.fragment_lengths.tsv"


def cutrun_qc_summary_path():
    return f"{OUTDIR}/qc/cutrun/cutrun_qc_summary.tsv"


def cutrun_fragment_lengths_multiqc_path():
    return f"{OUTDIR}/qc/cutrun/fragment_lengths_mqc.json"


def cutrun_qc_paths(run_sample):
    if not CUTRUN_QC_ENABLED:
        return []
    return [
        cutrun_usable_fragments_path(run_sample),
        fragment_length_path(run_sample),
    ]


# --- spike-in paths -------------------------------------------------------
def spikein_bam_path(run_sample):
    return f"{OUTDIR}/spikein/bowtie2/{run_sample}/{run_sample}.spikein.bam"


def spikein_bai_path(run_sample):
    return spikein_bam_path(run_sample) + ".bai"


def spikein_flagstat_path(run_sample):
    return f"{OUTDIR}/qc/spikein/{run_sample}/{run_sample}.spikein.flagstat.txt"


def spikein_scale_factor_path(run_sample):
    return f"{OUTDIR}/qc/spikein/{run_sample}/{run_sample}.spikein_scale.tsv"


def spikein_scale_factors_path():
    return f"{OUTDIR}/qc/spikein/spikein_scale_factors.tsv"


def spikein_normalized_bigwig_path(run_sample):
    return f"{OUTDIR}/bigwig/{run_sample}/{run_sample}.spikein_normalized.bw"


def spikein_qc_paths(run_sample):
    if not SPIKEIN_ENABLED:
        return []
    return [
        spikein_flagstat_path(run_sample),
        spikein_scale_factor_path(run_sample),
    ]


# --- bigwig paths ----------------------------------------------------------
def raw_bigwig_path(run_sample):
    return f"{OUTDIR}/bigwig/{run_sample}/{run_sample}.raw.bw"


def normalized_bigwig_path(run_sample):
    return f"{OUTDIR}/bigwig/{run_sample}/{run_sample}.normalized.bw"


def pooled_normalized_bigwig_path(sample):
    return f"{OUTDIR}/bigwig/pooled/{sample}/{sample}.pooled.normalized.bw"


def genome_size_path():
    return f"{OUTDIR}/qc/genome_size_{ANALYSIS_MODE}.txt"


# --- macs3 paths ----------------------------------------------------------
def macs3_peak_prefix(run_sample, pooled=False):
    if pooled:
        return f"{OUTDIR}/macs3/pooled/{run_sample}/{run_sample}.pooled"
    return f"{OUTDIR}/macs3/{run_sample}/{run_sample}"


def macs3_xls_path(run_sample, pooled=False):
    return macs3_peak_prefix(run_sample, pooled) + "_peaks.xls"


def macs3_narrow_path(run_sample, pooled=False):
    return macs3_peak_prefix(run_sample, pooled) + "_peaks.narrowPeak"


def macs3_summits_path(run_sample, pooled=False):
    return macs3_peak_prefix(run_sample, pooled) + "_summits.bed"


def macs3_consensus_path(sample):
    return f"{OUTDIR}/macs3/consensus/{sample}.consensus.narrowPeak"


def replicate_concordance_path(sample):
    return f"{OUTDIR}/qc/replicates/{sample}.replicate_concordance.tsv"


def frip_path(run_sample, pooled=False):
    if pooled:
        return f"{OUTDIR}/qc/frip/{run_sample}/{run_sample}.pooled.frip.tsv"
    return f"{OUTDIR}/qc/frip/{run_sample}/{run_sample}.frip.tsv"


def frip_summary_path():
    return f"{OUTDIR}/qc/frip/frip_summary.tsv"


def pooled_control_bam(sample):
    control = config["samples"].get(sample, {}).get("control")
    if control:
        return pooled_final_bam_path(control)
    return None


# Validation
# ----------------------------------------------------------------------------
def validate_workflow_config():
    if READ_TYPE not in {"PE", "SE"}:
        raise ValueError("read_type must be PE or SE.")
    if not SAMPLES:
        raise ValueError("config.samples must contain at least one sample.")

    for sample, sample_cfg in config.get("samples", {}).items():
        reps = sample_cfg.get("replicates")
        if not isinstance(reps, dict) or not reps:
            raise ValueError(f"samples.{sample}.replicates must be a non-empty mapping.")
        for rep, rep_cfg in reps.items():
            r1 = rep_cfg.get("R1", [])
            if not isinstance(r1, list) or not r1:
                raise ValueError(f"samples.{sample}.replicates.{rep}.R1 must be a non-empty list.")
            if is_single_end():
                if rep_cfg.get("R2"):
                    raise ValueError(f"samples.{sample}.replicates.{rep}.R2 is set but read_type is SE.")
            else:
                r2 = rep_cfg.get("R2", [])
                if not isinstance(r2, list) or not r2:
                    raise ValueError(f"samples.{sample}.replicates.{rep}.R2 must be a non-empty list for PE data.")
                if len(r1) != len(r2):
                    raise ValueError(
                        f"samples.{sample}.replicates.{rep}.R1 and R2 must have the same number of FASTQs."
                    )

    if FILTER_BLACKLIST and not blacklist_path():
        raise ValueError(
            "filter_blacklist is enabled but no blacklist.path or blacklist.url is set."
        )

    if SPIKEIN_ENABLED:
        spikein = config.get("spikein", {})
        if not spikein.get("bowtie2_index"):
            raise ValueError("spikein.enabled is true but spikein.bowtie2_index is missing.")
        if float(spikein.get("scale_to", 10000)) <= 0:
            raise ValueError("spikein.scale_to must be greater than 0.")

    if MACS3_ENABLED:
        for key in ["gsize", "qvalue", "extsize"]:
            if key not in MACS3_CONFIG:
                raise ValueError(f"macs3.enabled is true but macs3.{key} is missing.")

    if ANALYSIS_MODE not in {"standard", "te"}:
        raise ValueError("analysis_mode must be 'standard' or 'te'.")

    if TE_ENABLED:
        if TE_BLACKLIST_MODE not in {"none", "full"}:
            raise ValueError("te.blacklist must be one of: none, full.")
        if TE_BLACKLIST_MODE == "full" and not blacklist_path():
            raise ValueError(
                f"te.blacklist is '{TE_BLACKLIST_MODE}' but no blacklist.path or "
                "blacklist.url is set."
            )
        if TE_ANNOTATION and TE_ANNOTATION.startswith("path/to"):
            raise ValueError("te.annotation is still the placeholder path.")

    for sample, sample_cfg in config.get("samples", {}).items():
        control = sample_cfg.get("control")
        if control is not None and control not in config["samples"]:
            raise ValueError(
                f"samples.{sample}.control '{control}' is not in samples."
            )


def check_input_files():
    """Fail fast when referenced files do not exist (checked at DAG time)."""
    ref = config.get("reference", {})
    for key in ("bowtie2_index", "chrom_sizes", "fasta"):
        path = ref.get(key)
        if not path:
            continue
        if key == "bowtie2_index":
            for suffix in (".1.bt2", ".1.bt2l"):
                if os.path.exists(path + suffix):
                    break
            else:
                raise ValueError(f"reference.bowtie2_index does not exist: {path}")
        elif not os.path.exists(path):
            raise ValueError(f"reference.{key} does not exist: {path}")
    explicit_blacklist = blacklist_config().get("path")
    if FILTER_BLACKLIST and explicit_blacklist and not explicit_blacklist.startswith("/dev/null"):
        if not os.path.exists(explicit_blacklist):
            raise ValueError(f"blacklist.path does not exist: {explicit_blacklist}")
    # The unified config documents TE resources even in standard mode. Do not
    # require those files until the TE branch is explicitly enabled.
    if TE_ENABLED:
        if not TE_ANNOTATION:
            raise ValueError("analysis_mode is 'te' but te.annotation is missing.")
        if not os.path.exists(TE_ANNOTATION):
            raise ValueError(f"te.annotation does not exist: {TE_ANNOTATION}")
        dfam_metadata = TE_CONFIG.get("dfam_metadata")
        if dfam_metadata and not os.path.exists(dfam_metadata):
            raise ValueError(f"te.dfam_metadata does not exist: {dfam_metadata}")
    for sample in SAMPLES:
        for rep, rep_cfg in config["samples"][sample]["replicates"].items():
            for key in ("R1", "R2"):
                for fastq in rep_cfg.get(key, []):
                    if not os.path.exists(fastq):
                        raise ValueError(f"samples.{sample}.{rep}.{key} does not exist: {fastq}")


validate_workflow_config()


def genome_size_exclude_chrom_re():
    """Regex used to exclude chrM/scaffolds from bigwig/effective size."""
    return r"(?i)(^|\|)(chrM|MT|M)(\||$)|\b(scaffold|random|un_|_un|alt|fix|hap)\b|_mapped_Scaffold"


def effective_genome_size_cmd():
    """Shell command writing the effective genome size to a file.

    Uses only the contigs that are kept in the bigWig tracks so that
    RPKM-style normalization is not inflated by excluded contigs.
    """
    remove_chrm = config.get("bigwig", {}).get("remove_chrM_and_scaffolds", True)
    remove_opt = "--remove-chrm" if remove_chrm else ""
    return (
        f"python workflow/rules/scripts/compute_effective_genome_size.py "
        f"--chrom-sizes {config['reference']['chrom_sizes']} {remove_opt} "
        f"--out {genome_size_path()}"
    )


def exclude_bed_commands():
    """Shell fragment producing $EXCLUDE_BED (contigs to omit from bigWig)."""
    chrom_sizes = config["reference"]["chrom_sizes"]
    if not config.get("bigwig", {}).get("remove_chrM_and_scaffolds", True):
        return 'EXCLUDE_BED=""'
    return (
        'EXCLUDE_BED=$(mktemp)\n'
        "awk 'BEGIN{IGNORECASE=1} {\n"
        "    chrom=$1; size=$2;\n"
        "    if (chrom ~ /_/ || chrom ~ /scaffold|random|un|alt|fix|hap/ || "
        'chrom == "chrM" || chrom == "MT" || chrom == "M") {\n'
        '        print chrom"\\t0\\t"size\n'
        "    }\n"
        "}' " + chrom_sizes + ' > "$EXCLUDE_BED"\n'
        'if [ ! -s "$EXCLUDE_BED" ]; then EXCLUDE_BED=""; fi'
    )


rule faidx_reference:
    input:
        fa=config["reference"]["fasta"]
    output:
        fai=config["reference"]["fasta"] + ".fai"
    threads: 2
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools faidx {input.fa}"

check_input_files()

if FILTER_BLACKLIST and blacklist_config().get("url"):
    rule get_blacklist:
        output:
            bed=blacklist_path()
        conda:
            "envs/samtools.yaml"
        script:
            "scripts/get_blacklist.py"
