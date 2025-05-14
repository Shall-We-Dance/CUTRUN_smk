# Snakemake Workflow for CUT&RUN Upstream Analysis

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

This repository contains a modular and scalable [Snakemake](https://github.com/snakemake/snakemake) workflow for analyzing CUT&RUN (or ChIP-seq) data.

## ✨ Features

- 🧠 **Automatic Sample Detection**  
  Supports various naming conventions including `_R1.fastq.gz`, `_1.fastq.gz`, `.fastq.gz`, `_R1_001.fastq.gz`, etc.

- 🔁 **SE/PE Mode Auto-Detection**  
  Automatically routes samples through the correct pipeline depending on whether data is single-end or paired-end.

- ⚙️ **Flexible and Configurable**  
  Centralized `config.yaml` to set input paths, number of threads, STAR index, genome size, bin size, and more.

- 🧬 **Multimapping Handling**  
  Retains multi-mapping reads during STAR alignment, and includes a post-mapping `multimap_weight` function to adjust for `NH` tag weights (for accurate peak calling).

- 🚫 **Blacklist Filtering (Optional)**  
  When `filter_blacklist: true` is set in `config.yaml`, ENCODE blacklist regions will be automatically downloaded (based on genome) and applied to `bamCoverage` using `--blackListFileName`. This step replaces the older repeat masking logic.

- 📊 **BigWig Generation with Normalization**  
  Converts BAM to bigWig using `deeptools` with and without normalization (e.g., RPKM), while excluding PCR duplicates (`--samFlagExclude 1024`).

---

## 🧬 Workflow Overview

1. **Sample Detection**  
   Automatically detects sample names based on filenames.

2. **Quality Trimming**  
   Uses `fastp` to trim adapters and remove low-quality reads.

3. **Alignment**  
   Aligns reads to the reference genome using `STAR`, retaining up to 100 multi-mapped hits.

4. **Multimap Weighting**  
   Applies fractional weighting to multi-mapped reads based on their `NH` tag values.

5. **Blacklist Filtering (Optional)**  
   Filters signal from known artefact regions via ENCODE blacklist when `filter_blacklist` is enabled.

6. **BigWig Conversion**  
   Generates normalized (`RPKM`) and unnormalized bigWig files for visualization.

7. **Peak Calling**  
   Uses `MACS3` to call peaks from the aligned BAM files. 
   
---

## 🚀 Quick Start

1. Clone the repository:

```bash
git clone https://github.com/Shall-We-Dance/CUTRUN_smk.git
cd CUTRUN_smk
```

2. Edit `config.yaml` to specify your paths and parameters.

3. Activate Snakemake and run the pipeline:

```bash
snakemake --use-conda --cores 16
```

---

## 📁 Project Structure

```
CUTRUN_smk/
├── config/
│   └── config.yaml              # Main configuration file
│
├── workflow/
│   ├── Snakefile                # Entry point Snakefile
│   ├── rules/                   # Modular rule files
│   │   ├── fastp.smk
│   │   ├── star.smk
│   │   ├── macs3.smk
│   │   ├── bam_to_bigwig.smk
│   │   └── detect_samples.smk
│   └── envs/                    # Conda environments
│       ├── fastp.yaml
│       ├── star.yaml
│       ├── macs3.yaml
│       ├── bedtools.yaml
│       └── deeptools.yaml
│
├── results/                     # Final and intermediate output files
│   ├── fastp/
│   ├── star/
│   └── bigwig/
│
├── logs/                        # Log files for each step
│   ├── fastp/
│   ├── star/
│   └── bigwig/
│
├── .gitignore
├── LICENSE
└── README.md
```

---

## 📝 Notes

* STAR genome index must be prebuilt.
* For blacklist functionality, genome name must match those recognized by ENCODE (e.g., `hg38`, `mm10`).
* `samtools`, `deeptools`, and other tools will auto-scale to the number of available threads (default `max/4`).

---

## License

MIT License

---

## Contact

For questions, issues, or contributions, please open an issue or pull request on GitHub.
