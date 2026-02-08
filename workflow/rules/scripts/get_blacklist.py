# workflow/scripts/get_blacklist.py
import gzip
import os
import shutil
import urllib.error
import urllib.request

blacklist_config = snakemake.config.get("blacklist", {})
url = blacklist_config.get("url")
out_bed = snakemake.output.bed

if not url:
    raise ValueError("No blacklist.url provided; cannot download blacklist.")

os.makedirs(os.path.dirname(out_bed), exist_ok=True)

print("[get_blacklist] Downloading blacklist:")
print(f"  URL: {url}")
print(f"  Output: {out_bed}")

tmp_path = out_bed + ".gz" if url.endswith(".gz") else out_bed + ".tmp"
try:
    urllib.request.urlretrieve(url, tmp_path)
except (urllib.error.URLError, urllib.error.HTTPError) as exc:
    raise RuntimeError(f"Failed to download blacklist from {url}: {exc}")

if tmp_path.endswith(".gz"):
    print("[get_blacklist] Unzipping downloaded blacklist.")
    with gzip.open(tmp_path, "rb") as f_in, open(out_bed, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(tmp_path)
else:
    shutil.move(tmp_path, out_bed)

print("[get_blacklist] Done.")
