# rules/bowtie2.smk

# This rule is for paired-end reads
rule bowtie2_pe:
    input:
        r1 = "results/fastp/{sample}_pe_fastp_R1.fastq.gz",
        r2 = "results/fastp/{sample}_pe_fastp_R2.fastq.gz"
    params:
        genome = config["bowtie2_index"]
    output:
        bam = "results/bowtie2/{sample}/{sample}_pe.bam",
    log:
        "logs/bowtie2/{sample}.log"
    threads: config["threads"]
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        bowtie2 \
            -x {params.genome} \
            -1 {input.r1} \
            -2 {input.r2} \
            --end-to-end \
            --very-sensitive \
            -X 500 \
            --no-mixed \
            --no-discordant \
            -p {threads} | \
            samtools view -bS - | \
            samtools sort -@ {threads} -o {output.bam} -
        echo "Running samtools index for {wildcards.sample}..." >> {log} 2>&1
        samtools index -@ {threads} {output.bam} >> {log} 2>&1
        """

# This rule is for single-end reads
rule biwtie2_se:
    input:
        r1 = "results/fastp/{sample}_se_fastp.fastq.gz"
    params:
        genome = config["bowtie2_index"]
    output:
        bam = "results/bowtie2/{sample}/{sample}_se.bam",
    log:
        "logs/bowtie2/{sample}.log"
    threads: config["threads"]
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2 \
            -x {params.genome} \
            -U {input.r1} \
            --end-to-end \
            --very-sensitive \
            -p {threads} | \
            samtools view -bS - | \
            samtools sort -@ {threads} -o {output.bam} -
        echo "Running samtools index for {wildcards.sample}..." >> {log} 2>&1
        samtools index -@ {threads} {output.bam} >> {log} 2>&1
        """

def get_raw_bam(sample):
    if sample in PE_SAMPLES:
        return f"results/bowtie2/{sample}/{sample}_pe.bam"
    else:
        return f"results/bowtie2/{sample}/{sample}_se.bam"


def get_filtered_bam(sample):
        return f"results/bowtie2/{sample}/{sample}.filtered.bam"


#def get_analysis_bam(sample):
#    if config.get("filter_blacklist", False):
#        return get_filtered_bam(sample)
#    else:
#        return get_raw_bam(sample)

def get_blacklist_path(genome):
    """
    Download and cache ENCODE blacklist BED file for a given genome.
    """
    url_dict = {
        "mm10": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz",
        "hg19": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz",
        "hg38": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz",
        "ce11": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce11-blacklist.v2.bed.gz",
        "ce10": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce10-blacklist.v2.bed.gz",
        "dm3": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm3-blacklist.v2.bed.gz",
        "dm6": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.v2.bed.gz",
    }
    if genome not in url_dict:
        raise ValueError(f"Unsupported genome: {genome}")

    cache_dir = "resources/blacklist"
    os.makedirs(cache_dir, exist_ok=True)

    url = url_dict[genome]
    gz_path = os.path.join(cache_dir, os.path.basename(url))
    bed_path = gz_path.replace(".gz", "")

    if not os.path.exists(bed_path):
        if not os.path.exists(gz_path):
            print(f"Downloading blacklist for {genome}...")
            urllib.request.urlretrieve(url, gz_path)

        print(f"Unzipping blacklist for {genome}...")
        with gzip.open(gz_path, "rb") as f_in, open(bed_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    return bed_path

rule filter_blacklist_bam:
    input:
        bam = lambda wildcards: get_raw_bam(wildcards.sample),
    output:
        bam = "results/bowtie2/{sample}/{sample}.filtered.bam",
    params:
        blacklist = lambda wildcards: get_blacklist_path(config["genome"])
    threads: config["threads"]
    log:
        "logs/bowtie2/{sample}.filter.log"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        r"""
        set -euo pipefail
        bedtools intersect -v \
            -abam {input.bam} \
            -b {params.blacklist} \
            -nonamecheck | \
        samtools sort -@ {threads} -o {output.bam} >> {log} 2>&1

        samtools index -@ {threads} {output.bam} >> {log} 2>&1
        """