# rules/detect_samples.smk

import os
import re
import glob
from snakemake.io import glob_wildcards

def detect_samples(raw_dir, recursive=False):
    patterns = [
        "*_R1.fastq.gz", "*_1.fastq.gz", "*_R1_001.fastq.gz",
        "*_R1.fq.gz", "*_1.fq.gz", "*_R1_001.fq.gz",
        "*.fastq.gz", "*.fq.gz"
    ]
    sample_set = set()
    for pattern in patterns:
        for filepath in glob.glob(os.path.join(raw_dir, "**" if recursive else "", pattern), recursive=recursive):
            filename = os.path.basename(filepath)
            if re.search(r"(_R2|_2|_R2_001)\.f(ast)?q\.gz$", filename):
                continue
            name = re.sub(r"(_R1|_1|_R1_001)?\.f(ast)?q\.gz$", "", filename)
            sample_set.add(name)
    return sorted(sample_set)

def get_sample_list(config):
    recursive = config.get("auto_detect_subdirs", False)
    if config["samples"] is None or config["samples"] == "all":
        samples = detect_samples(config["raw_dir"], recursive=recursive)
        print(f"Detected samples: {samples}")
    elif isinstance(config["samples"], list):
        samples = config["samples"]
        print(f"Using provided samples: {samples}")
    elif isinstance(config["samples"], str):
        samples = [config["samples"]]
        print(f"Using provided sample: {samples}")
    else:
        raise ValueError("Invalid value for config['samples']")
    return samples

def get_sample_type_map(raw_dir, samples=None):
    patterns_r1 = [
        "*_R1.fastq.gz", "*_1.fastq.gz", "*_R1_001.fastq.gz",
        "*_R1.fq.gz", "*_1.fq.gz", "*_R1_001.fq.gz",
        "*.fastq.gz", "*.fq.gz"
    ]
    sample_types_detected = {}
    recursive = config.get("auto_detect_subdirs", False)

    if samples is None:
        samples = detect_samples(raw_dir)

    for sample in samples:
        # Check both SE and PE patterns
        r1_patterns = [
            f"{sample}_R1.fastq.gz", f"{sample}_1.fastq.gz", f"{sample}_R1_001.fastq.gz",
            f"{sample}_R1.fq.gz", f"{sample}_1.fq.gz", f"{sample}_R1_001.fq.gz",
            f"{sample}.fastq.gz", f"{sample}.fq.gz"
        ]
        r2_patterns = [
            f"{sample}_R2.fastq.gz", f"{sample}_2.fastq.gz", f"{sample}_R2_001.fastq.gz",
            f"{sample}_R2.fq.gz", f"{sample}_2.fq.gz", f"{sample}_R2_001.fq.gz"
        ]

        found_r1 = False
        found_r2 = False
        for pat in r1_patterns:
            if glob.glob(os.path.join(raw_dir, "**" if recursive else "", pat), recursive=recursive):
                found_r1 = True
                break
        for pat in r2_patterns:
            if glob.glob(os.path.join(raw_dir, "**" if recursive else "", pat), recursive=recursive):
                found_r2 = True
                break

        if found_r1:
            sample_types_detected[sample] = "PE" if found_r2 else "SE"
        else:
            print(f"⚠️ Warning: Sample {sample} not found in {raw_dir}")

    return sample_types_detected
