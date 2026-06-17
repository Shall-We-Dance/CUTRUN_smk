# CUT&RUN/CUT&Tag/ChIP-seq Snakemake Pipeline

This repository provides a modular Snakemake workflow for paired-end and single-end CUT&RUN/CUT&Tag/ChIP-seq upstream processing. The pipeline focuses on QC, alignment, optional blacklist filtering, optional duplicate removal, and bigWig generation.

## Features

- **Per-lane fastp QC + sample-level reports**
- **Bowtie2 alignment with unique-mapper filtering**
- **Optional ENCODE/Boyle-Lab (or custom) blacklist filtering**
- **Optional PCR duplicate removal (Picard or samtools)**
- **CUT&RUN/CUT&Tag QC: usable fragments and fragment lengths**
- **Optional spike-in normalization and spike-in-scaled bigWigs**
- **Raw + normalized bigWig outputs**
- **MACS3 peak calling (optional)**
- **MultiQC summary report**

## Quick Start

```bash
git clone https://github.com/Shall-We-Dance/CUTRUN_smk.git
cd CUTRUN_smk

# edit config.yaml
snakemake --use-conda --cores 16
```

## Configuration (`config.yaml`)

Key sections:

- **`reference`**: bowtie2 index, chromosome sizes, fasta.
- **`read_type`**: set to `PE` (paired-end, default) or `SE` (single-end).
- **`samples`**: FASTQ lists per sample. For `SE`, provide only `R1` entries.
- **`filter_blacklist`**: enable/disable blacklist filtering.
- **`blacklist`**: manual input for blacklist files.
- **`remove_duplicates`**: enable/disable PCR duplicate removal.
- **`fastp.dedup_adapter.dedup`**: enable/disable fastp dedup for the main output workflow.
- **`bigwig`**: bin size + normalization + optional removal of chrM/scaffolds in bigWig output.
- **`cutrun_qc`**: usable fragment/read depth and fragment length outputs.
- **`spikein`**: optional spike-in genome alignment, scale factors, and spike-in-scaled bigWigs.
- **`output.dir`**: output directory (default `results`).
- **`macs3`**: MACS3 peak calling toggle and parameters.

BigWig-specific option:
- `bigwig.remove_chrM_and_scaffolds`: when `true` (default), excludes `chrM`/`MT`/`M` and scaffold-like contigs (`_`, `scaffold`, `random`, `Un`, `alt`, `fix`, `hap`) during bigWig generation.

### Blacklist (Manual Input)

Blacklist files are **manual**. You can either:

1. **Use a local file**
   ```yaml
   blacklist:
     path: "/path/to/custom.blacklist.bed"
   ```

2. **Download from a URL**
   ```yaml
   blacklist:
     url: "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.v2.bed.gz"
     cache_dir: "resources/blacklist"
   ```

Reference repositories:
- https://github.com/Boyle-Lab/Blacklist
- https://github.com/Boyle-Lab/Blacklist/tree/master/lists

## Workflow Steps

1. **fastp QC per lane** → unit-level HTML/JSON reports
2. **Merge lanes** → per-sample FASTQs
3. **Sample-level fastp report-only** → merged HTML/JSON
4. **Bowtie2 alignment**
5. **Unique mapper filtering**
6. **Optional blacklist filtering**
7. **Optional duplicate removal**
8. **CUT&RUN/CUT&Tag QC** (usable fragments, fragment length)
9. **Optional spike-in alignment and scaling**
10. **BigWig generation** (raw + normalized + optional spike-in scaled)
11. **MACS3 peak calling** (optional)
12. **MultiQC report**

## Outputs (default `results/`)

```
results/
├── qc/
│   ├── fastp/
│   ├── cutrun/
│   ├── spikein/
│   └── multiqc/
├── bowtie2/
│   └── <sample>/<sample>.<suffix>.bam
├── spikein/
│   └── bowtie2/<sample>/<sample>.spikein.bam
├── macs3/
│   └── <sample>/<sample>_peaks.narrowPeak
└── bigwig/
    └── <sample>/<sample>.raw.bw
    └── <sample>/<sample>.normalized.bw
    └── <sample>/<sample>.spikein_normalized.bw
```

## Notes

- The pipeline supports both **paired-end** (`read_type: PE`) and **single-end** (`read_type: SE`) samples.
- When `filter_blacklist: true`, you must provide **either** `blacklist.path` **or** `blacklist.url`.
- When `spikein.enabled: true`, `spikein.bowtie2_index` must point to a spike-in Bowtie2 index.
- MultiQC scans `results/qc`, `results/bowtie2`, and `logs/` by default.

## License

MIT License
