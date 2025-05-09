# Snakemake workflow for CUT&RUN Upstream Analysis

This repository contains a modular and scalable [Snakemake](https://github.com/snakemake/snakemake) workflow for analyzing CUT&RUN data.

## Features

- 🧠 **Automatic Sample Detection**  
  Supports various naming conventions including `_R1.fastq.gz`, `_1.fastq.gz`, `.fastq.gz`, etc.

- 🔁 **SE/PE Mode Auto-Detection**  
  Automatically routes samples through the correct pipeline depending on whether data is single-end or paired-end.

- 🎛️ **Customizable Configuration**  
  Centralized `config.yaml` for specifying input paths, threading, STAR index, repeat BED file, and more.

- 🔍 **Repeat Masking (Optional)**  
  If enabled in the config, `bedtools` will be used to filter out alignments overlapping repeat regions.

- 📊 **BigWig Generation**  
  Supports BAM-to-bigWig conversion with and without normalization (e.g., RPKM).

## Workflow Overview

1. **Sample Detection**  
   Automatically detects sample names based on filename patterns.

2. **Quality Trimming**  
   Uses `fastp` to trim adapters and filter reads.

3. **Alignment**  
   Maps reads to the reference genome using `STAR`.

4. **Repeat Masking (Optional)**  
   Removes reads overlapping with specified repeat regions.

5. **BigWig Conversion**  
   Converts BAM files to bigWig tracks for visualization.

## Usage

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
│   │   ├── bam_to_bigwig.smk
│   │   └── detect_samples.smk
│   └── envs/                    # Conda environments
│       ├── fastp.yaml
│       ├── star.yaml
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
│   └── star/
│
├── .gitignore
├── LICENSE
└── README.md
```

## License

MIT License

## Contact

For questions, issues, or contributions, please open an issue or pull request on GitHub.
