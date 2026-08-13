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

## Sample and control model

Each sample contains one or more biological replicates. Multiple FASTQ paths
inside a replicate are sequencing lanes and are concatenated in list order.
For PE data, R1 and R2 lists must have equal lengths and matching lane order.
Sample and replicate names must not contain a period.

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
invalid mappings fail before DAG construction.

## Standard versus TE mode

Standard mode uses the MAPQ-filtered, blacklist-filtered, deduplicated branch
as its primary result. If that final BAM has already been coordinate-
deduplicated, keep `macs3.keep_dup: all`; otherwise MACS3 applies a second
coordinate-level duplicate filter.

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

## Tests

The tiny integration fixtures contain synthetic sequences only.

```bash
# Unit tests require samtools on PATH (or SAMTOOLS=/path/to/samtools).
python tests/test_te_quant.py
python tests/test_te_quant_pe.py
python tests/test_frip.py

# End-to-end PE standard, PE TE, and SE smoke tests.
snakemake --snakefile workflow/Snakefile \
  --configfile tests/config.tiny.yaml --use-conda --cores 4
snakemake --snakefile workflow/Snakefile \
  --configfile tests/config.tiny_te.yaml --use-conda --cores 4
snakemake --snakefile workflow/Snakefile \
  --configfile tests/config.tiny_se.yaml --use-conda --cores 4
```

The unit tests cover SE and real paired-end 1/N weighting, nested TE records,
secondary alignments, fragment-level FRiP, orphan mates, and exact BED/SAM
coordinate boundaries.

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
