# workflow/scripts/get_blacklist.py
import gzip
import os
import shutil
import urllib.error
import urllib.request
from contextlib import suppress

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


def download_blacklist(download_url, destination, retries=3, timeout=30):
    last_error = None
    headers = {"User-Agent": "CUTRUN_smk/blacklist-fetch"}
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(download_url, headers=headers)
            with urllib.request.urlopen(request, timeout=timeout) as response, open(
                destination, "wb"
            ) as f_out:
                shutil.copyfileobj(response, f_out)
            if os.path.getsize(destination) == 0:
                raise RuntimeError("downloaded file is empty")
            return
        except (urllib.error.URLError, urllib.error.HTTPError, RuntimeError) as exc:
            last_error = exc
            with suppress(FileNotFoundError):
                os.remove(destination)
            if attempt < retries:
                print(f"[get_blacklist] Download attempt {attempt} failed: {exc}")
                print("[get_blacklist] Retrying...")
    raise RuntimeError(f"Failed to download blacklist from {download_url}: {last_error}")


download_blacklist(url, tmp_path)

if tmp_path.endswith(".gz"):
    print("[get_blacklist] Unzipping downloaded blacklist.")
    try:
        with gzip.open(tmp_path, "rb") as f_in, open(out_bed, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    except OSError as exc:
        raise RuntimeError(
            "Downloaded blacklist is not a valid gzip file. "
            f"Check the URL or network/proxy settings. Original error: {exc}"
        )
    finally:
        with suppress(FileNotFoundError):
            os.remove(tmp_path)
else:
    shutil.move(tmp_path, out_bed)

print("[get_blacklist] Done.")
