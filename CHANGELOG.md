# Changelog

## [v1.1.0] - 2025-05-14

### 🚀 Performance Improvements
- Default CPU threads changed from `max/8` → `max/4` for more efficient usage on multi-core systems.
- Enabled multi-threading in all `samtools` steps via `-@ {threads}`.

### 🎯 Multi-Mapping Support
- `--outSAMmultNmax` in STAR increased from `1` to `100`, improving retention of multi-mapped reads.
- Introduced `multimap_weight` function to post-process alignments with `NH` tags, preparing for downstream MACS3 integration.

### 📊 BigWig Generation
- `--samFlagExclude 1024` added to `bamCoverage` to exclude PCR duplicates from signal tracks.
- New `filter_blacklist` option in `config.yaml`:  
  When set to `true`, automatically downloads and applies ENCODE blacklist BED files (based on genome) during bigWig generation with `--blackListFileName`.

### 🧼 Cleanup & Removal
- Deprecated `mask_repeats` step. Use `filter_blacklist` instead for repeat/artefact exclusion.

### 🧠 Usability Enhancements
- `detect_samples` now supports filenames like `*_R1_001.fq.gz`, improving compatibility with common Illumina naming formats.

### 🧬 Peak Calling
- Integrated **MACS3** for peak calling from aligned BAM files.

---

## [v1.0.0] - 2025-05-09

### ✨ Initial Release
- End-to-end CUT&RUN/ChIP-seq pipeline based on Snakemake.
- Modules included:
  - Read trimming (fastp)
  - Alignment (STAR with genome index)
  - BAM filtering and sorting
  - bigWig generation (deeptools)
  - Quality control reports
- Configurable via `config.yaml` (samples, genome, bin size, etc.)
- Designed for paired-end and single-end data.
- Initial environment support via Conda (deeptools, samtools, STAR).

---

