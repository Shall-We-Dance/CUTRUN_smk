# CUT&RUN / CUT&Tag / ChIP-seq Snakemake pipeline

A replicate-aware Snakemake workflow for paired-end or single-end chromatin
profiling. It supports a conventional unique-mapper workflow and an optional
TE-aware mode that retains multimappers and reports family-level 1/N-weighted
enrichment against matched input controls.

## What it does

- merges sequencing lanes within each biological replicate;
- runs fastp QC and Bowtie2 alignment;
- creates a standard MAPQ-filtered branch;
- optionally removes blacklist-overlapping templates and PCR duplicates;
- reports usable fragments, fragment lengths, true template-level FRiP,
  BigWigs, MACS3 peaks, replicate concordance, and MultiQC;
- optionally aligns spike-in reads and creates spike-in-scaled tracks;
- in TE mode, retains secondary alignments, assigns each alignment to its
  best-overlapping TE record, weights every template by 1/N, and summarizes
  family/subclass enrichment using explicitly matched controls.

The repository contains no biological data, real sample identifiers,
accessions, local filesystem paths, results, or Snakemake runtime metadata.

## Requirements

- Linux
- Snakemake >= 9.3.3
- Conda or Mamba for per-rule environments

All bioinformatics dependencies are declared under `workflow/rules/envs/`.

## Quick start

1. Edit the thoroughly commented `config.yaml`. Set `analysis_mode` to
   `standard` (unique-mapper analysis) or `te` (adds the TE-aware branch).

2. Validate the DAG before starting a large run:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config.yaml \
  --use-conda --cores 16 --dry-run
```

3. Run the workflow:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config.yaml \
  --use-conda --cores 16 --rerun-incomplete --printshellcmds
```

For parameter-sensitivity runs, use a different `output.dir` for every run so
results and their associated logs cannot overwrite one another.

`assay_type` is optional and defaults to `cutrun`. Set it to `chipseq` for a
ChIP-seq design. This does not force an Input/control or replicate labels: both
remain optional for incomplete public datasets and simple experiments.

## Assay type, samples, lanes, and optional controls

Assay type and analysis mode are independent:

| Setting | Purpose | Default |
| --- | --- | --- |
| `assay_type: cutrun \| chipseq` | Describes the assay/sample design | `cutrun` |
| `analysis_mode: standard \| te` | Selects unique-read or TE-aware mapping | `standard` |

For a single biological sample, use the simple flat layout. A control/input and
a replicate label are both optional. Multiple FASTQs are sequencing lanes and
are concatenated in list order. For PE data, R1 and R2 lists must have equal
lengths and matching lane order.

```yaml
samples:
  Example_H3K4me3:
    R1:
      - data/example_L001_R1.fastq.gz
      - data/example_L002_R1.fastq.gz
    R2:
      - data/example_L001_R2.fastq.gz
      - data/example_L002_R2.fastq.gz
```

This produces outputs named with `Example_H3K4me3`, without an artificial
`.rep1` suffix. It also does not produce redundant `.pooled` outputs. MACS3
runs without `-c` when no control is configured.

Use the explicit layout only when biological replicates and/or matched controls
need to be represented:

```yaml
samples:
  Example_IP:
    control: Example_INPUT
    control_replicates:
      rep1: input_a
      rep2: input_b
    replicates:
      rep1:
        R1:
          - data/ip_rep1_L001_R1.fastq.gz
          - data/ip_rep1_L002_R1.fastq.gz
        R2:
          - data/ip_rep1_L001_R2.fastq.gz
          - data/ip_rep1_L002_R2.fastq.gz
      rep2:
        R1: [data/ip_rep2_L001_R1.fastq.gz]
        R2: [data/ip_rep2_L001_R2.fastq.gz]
  Example_INPUT:
    call_peaks: false
    replicates:
      input_a:
        R1:
          - data/input_a_L001_R1.fastq.gz
          - data/input_a_L002_R1.fastq.gz
        R2:
          - data/input_a_L001_R2.fastq.gz
          - data/input_a_L002_R2.fastq.gz
      input_b:
        R1: [data/input_b_R1.fastq.gz]
        R2: [data/input_b_R2.fastq.gz]
```

When IP and input replicate names are identical, omit `control_replicates`.
Set a mapping value to `pooled` only when pooling controls is justified by the
experimental design. The same mapping is used by MACS3 and TE enrichment, and
invalid mappings fail before DAG construction. Sample and explicit replicate
names must not contain a period.

Pooled fragment tables, normalized BigWigs, peaks, and FRiP rows are generated only for a
sample with at least two explicitly configured biological replicates. Multiple
FASTQ lanes within one replicate do not trigger pooling. Pooled BAMs are also
retained when `keep_bam: true`. If a replicated IP has a single-replicate
control, that control's fragment table is used directly without making a
redundant pooled control.

Without a control, TE overlap and family-count tables are still produced, but
ChIP/input enrichment tables contain headers only because enrichment cannot be
estimated without a comparator.

## Standard versus TE mode

Standard mode uses the MAPQ-filtered, optionally blacklist-filtered and
deduplicated branch as its primary result. MACS3, FRiP, quantitative BigWigs and
usable-fragment QC consume one shared fragment table. `macs3.keep_dup` must be
`all`; use upstream `remove_duplicates` for duplicate sensitivity runs.

MACS3 output follows the configured peak mode: `broad: false` produces
`narrowPeak` plus `summits.bed`, while `broad: true` produces `broadPeak` and
`gappedPeak`. Broad mode does not request summits because MACS3 forbids
combining `--broad` with `--call-summits`. FRiP, pooled peak calling, and
replicate concordance automatically use the matching primary peak file.

TE mode adds a separate all-alignment branch:

- Bowtie2 reports up to `te.k` alignments;
- secondary alignments are retained;
- the recommended primary TE configuration uses `te.blacklist: none` and
  `te.remove_duplicates: false`;
- `te.blacklist: full`, duplicate removal, overlap thresholds, and alternative
  `k` values should be separate sensitivity runs;
- the TE diagnostic BigWig counts alignments and is not a quantitative family
  abundance track;
- family tables use 1/N weights and matched-input RPM ratios.

`te_family_enrichment.curated.tsv` contains only names verified through the
configured Dfam metadata. Unverified names are placed in `unknown`, never
silently promoted to curated.

## Important interpretation limits

See [consensus, IDR and fragment-counting documentation](docs/consensus-and-fragments.md)
for exact matching rules, settings, outputs, and migration notes. Consensus
uses fixed pooled candidates and configurable k-of-n replicate support; broad
peaks also have an independently evaluated, basewise-supported core.

- A single biological replicate cannot establish reproducibility; consensus
  and direction-consistency are reported as `NA`.
- Symmetric peak overlap/Jaccard is QC, not a replacement for TF IDR.
- 1/N weighting is a transparent family-level estimator, not an EM/CSEM/T3E
  locus-resolution model. Use an orthogonal method and parameter sensitivity
  for publication-level TE conclusions.
- Pooled peaks are supplementary. Prefer replicate-supported or IDR peaks.
- Low FRiP or low replicate concordance is a data-quality warning, even when
  every workflow job succeeds.

## Main outputs

```text
<output.dir>/
├── bowtie2/                         per-replicate and pooled BAMs
├── fragments/                       shared BED3 molecule tables
├── bigwig/                          replicate and pooled tracks
├── macs3/                           replicate, pooled, and consensus peaks
├── qc/
│   ├── fastp/
│   ├── bam_stats/
│   ├── cutrun/
│   ├── frip/
│   ├── replicates/
│   └── multiqc/multiqc_report.html
└── te/                              TE mode only
    ├── <sample>.<rep>/*.te_overlap.tsv
    ├── <sample>.<rep>/*.te_family_counts.tsv
    ├── te_overlap_summary.tsv
    └── te_family_enrichment.{tsv,curated.tsv,novel.tsv,unknown.tsv,subclass.tsv}
```

Logs are written to `<output.dir>/logs/`.

The MultiQC report includes fastp, Bowtie2, processing-stage flagstat metrics,
duplicate metrics, CUT&RUN alignment/usable-fragment counts and span-length distributions,
and FRiP. Spike-in columns and TE-overlap QC are added when those analyses are
enabled. Custom CUT&RUN/FRiP/TE sections are checked after report generation so
a silently dropped MultiQC custom-content module fails the workflow.

## Tests

The tiny integration fixtures contain synthetic sequences only.

```bash
# Unit tests require samtools on PATH (or SAMTOOLS=/path/to/samtools).
python tests/test_te_quant.py
python tests/test_te_quant_pe.py
python tests/test_frip.py
python tests/test_multiqc_summaries.py
python tests/test_consensus.py
# Run with the deeptools rule environment (pyBigWig) and samtools on PATH.
python tests/test_fragment_coverage.py
# Run with the IDR environment; IDR=/path/to/idr also selects the executable.
# The real-fit test is skipped if IDR is not installed; the other tests still run.
python tests/test_idr.py
# DAG-only positive/negative configuration checks (Snakemake + PyYAML).
python tests/test_workflow_config.py

# End-to-end PE standard, PE TE, and SE smoke tests.
snakemake --snakefile workflow/Snakefile \
  --configfile tests/config.tiny.yaml --use-conda --cores 4
snakemake --snakefile workflow/Snakefile \
  --configfile tests/config.tiny_te.yaml --use-conda --cores 4
snakemake --snakefile workflow/Snakefile \
  --configfile tests/config.tiny_se.yaml --use-conda --cores 4

# Full regressions: pooling, mixed modes, three replicates/optional IDR,
# and spike-in counting plus MultiQC. These can also be checked with --dry-run.
snakemake --snakefile workflow/Snakefile \
  --configfile tests/config.tiny_replicates.yaml --use-conda --cores 4
snakemake --snakefile workflow/Snakefile \
  --configfile tests/config.tiny_broad.yaml --use-conda --cores 4
snakemake --snakefile workflow/Snakefile \
  --configfile tests/config.tiny_consensus.yaml --use-conda --cores 4
snakemake --snakefile workflow/Snakefile \
  --configfile tests/config.tiny_spikein.yaml --use-conda --cores 4
# After all seven integration runs, verify counts, every BigWig base value,
# pooled molecule multiplicity, candidate/core support, and IDR/MultiQC status.
# Use the deeptools rule environment (pyBigWig).
python tests/check_tiny_outputs.py
```

The unit tests cover SE and real paired-end 1/N weighting, nested TE records,
secondary alignments, fragment-level FRiP, orphan mates, and exact BED/SAM
coordinate boundaries. New tests additionally cover fixed candidate boundaries,
reciprocal/summit matching, per-replicate support, non-bridging broad cores,
shared BAM-to-fragment counts, BigWig values/normalizations, and real IDR fitting
on synthetic ranked peaks. The tiny sequencing fixture has too few peaks for a
meaningful IDR fit and must report `insufficient_peaks`, not success.

## Reproducibility recommendations

- Commit the exact config used for each analysis after replacing private paths
  with portable relative paths or a site-specific profile.
- Record the Git commit, config checksum, reference/annotation checksums,
  Snakemake version, and Conda environment exports with every result set.
- Keep sensitivity runs in separate output and log directories.
- Do not commit `.snakemake/`, results, logs, local configs, or private paths;
  the included `.gitignore` covers the usual cases.

## License

MIT. See `LICENSE`.
