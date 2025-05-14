# rules/detect_samples.smk

import os
import re
import glob
from snakemake.io import glob_wildcards

def detect_samples(raw_dir):
    patterns = [
        "*_R1.fastq.gz", "*_1.fastq.gz", "*_R1_001.fastq.gz",
        "*_R1.fq.gz", "*_1.fq.gz", "*_R1_001.fq.gz",
        "*.fastq.gz", "*.fq.gz"
    ]
    recursive = config.get("auto_detect_subdirs", False)
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

def get_sample_type_map(raw_dir):
    patterns_r1 = [
        "*_R1.fastq.gz", "*_1.fastq.gz", "*_R1_001.fastq.gz",
        "*_R1.fq.gz", "*_1.fq.gz", "*_R1_001.fq.gz",
        "*.fastq.gz", "*.fq.gz"
    ]
    sample_types_detected = {}
    recursive = config.get("auto_detect_subdirs", False)
    for pattern in patterns_r1:
        for filepath in glob.glob(os.path.join(raw_dir, "**" if recursive else "", pattern), recursive=recursive):
            filename = os.path.basename(filepath)
            if re.search(r"(_R2|_2|_R2_001)\.f(ast)?q\.gz$", filename):
                continue
            sample = re.sub(r"(_R1|_1|_R1_001)?\.f(ast)?q\.gz$", "", filename)
            dir_path = os.path.dirname(filepath)
            r2_patterns = [
                f"{sample}_R2.fastq.gz", f"{sample}_2.fastq.gz", f"{sample}_R2_001.fastq.gz",
                f"{sample}_R2.fq.gz", f"{sample}_2.fq.gz", f"{sample}_R2_001.fq.gz"
            ]
            has_r2 = any(
                glob.glob(os.path.join(dir_path, r2_pat))
                for r2_pat in r2_patterns
            )
            sample_types_detected[sample] = "PE" if has_r2 else "SE"
    return sample_types_detected
