# workflow/scripts/get_blacklist.py
import os
import gzip
import shutil
import urllib.request
import urllib.error

genome = snakemake.config["genome"]
cache_dir = snakemake.config.get("blacklist", {}).get("cache_dir", "resources/blacklist")
os.makedirs(cache_dir, exist_ok=True)

# Boyle-Lab blacklist repository (raw GitHub)
url_dict = {
    "mm10": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz",
    "mm39": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm39-blacklist.v2.bed.gz",
    "hg19": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz",
    "hg38": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz",
    "ce11": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce11-blacklist.bed.gz",
    "ce10": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce10-blacklist.bed.gz",
    "dm3": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm3-blacklist.bed.gz",
    "dm6": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.bed.gz",
}

if genome not in url_dict:
    raise ValueError(
        f"Unsupported genome: {genome}. "
        f"Available: {', '.join(sorted(url_dict.keys()))}"
    )

url = url_dict[genome]
gz_path = os.path.join(cache_dir, os.path.basename(url))
out_bed = snakemake.output.bed

# Download if needed
if not os.path.exists(gz_path):
    print(f"[get_blacklist] Downloading {genome} blacklist:")
    print(f"  URL: {url}")
    print(f"  -> {gz_path}")
    try:
        urllib.request.urlretrieve(url, gz_path)
    except (urllib.error.URLError, urllib.error.HTTPError) as e:
        raise RuntimeError(f"Failed to download blacklist from {url}: {e}")

# Unzip to requested output path
print(f"[get_blacklist] Unzipping:")
print(f"  {gz_path} -> {out_bed}")
os.makedirs(os.path.dirname(out_bed), exist_ok=True)

with gzip.open(gz_path, "rb") as f_in, open(out_bed, "wb") as f_out:
    shutil.copyfileobj(f_in, f_out)

print("[get_blacklist] Done.")
