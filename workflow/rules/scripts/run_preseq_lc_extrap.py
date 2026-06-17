from pathlib import Path
import shlex
import shutil
import subprocess
import sys


PRESEQ_HEADER = "TOTAL_READS\tEXPECTED_DISTINCT\n"
NONFATAL_PRESEQ_ERRORS = {
    "max count before zero is less than min required count": (
        "SKIPPED_INSUFFICIENT_COUNTS",
        "preseq could not fit a library complexity curve because the duplicate count distribution is insufficient",
    ),
}


def classify_preseq_failure(log_text):
    lower_log = log_text.lower()
    for needle, result in NONFATAL_PRESEQ_ERRORS.items():
        if needle in lower_log:
            return result
    return None


def clean_message(message):
    return " ".join(str(message).split())


def write_status(path, sample, status, message):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("sample\tpreseq_status\tpreseq_message\n")
        handle.write(f"{sample}\t{status}\t{clean_message(message)}\n")


def log_text(path):
    try:
        return Path(path).read_text(encoding="utf-8", errors="replace")
    except FileNotFoundError:
        return ""


def run_preseq_lc_extrap(smk):
    sample = smk.wildcards.sample
    input_bam = Path(smk.input.bam)
    output_txt = Path(smk.output.txt)
    status_tsv = Path(smk.output.status)
    log_path = Path(smk.log[0])
    tmp_output = Path(f"{output_txt}.tmp")

    output_txt.parent.mkdir(parents=True, exist_ok=True)
    status_tsv.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_output.unlink(missing_ok=True)

    cmd = ["preseq", "lc_extrap"]
    pe_arg = str(smk.params.pe_arg)
    if pe_arg:
        cmd.append(pe_arg)
    cmd.extend(["-output", str(tmp_output), str(input_bam)])

    with open(log_path, "w", encoding="utf-8") as log_handle:
        log_handle.write("$ " + " ".join(shlex.quote(arg) for arg in cmd) + "\n")
        log_handle.flush()
        proc = subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT)

    if proc.returncode == 0:
        if not tmp_output.exists() or tmp_output.stat().st_size == 0:
            message = "preseq exited successfully but did not write a non-empty lc_extrap output"
            write_status(status_tsv, sample, "FAILED_EMPTY_OUTPUT", message)
            with open(log_path, "a", encoding="utf-8") as log_handle:
                log_handle.write(f"\nERROR: {message}\n")
            return 1
        shutil.move(str(tmp_output), output_txt)
        write_status(status_tsv, sample, "OK", "preseq lc_extrap completed")
        return 0

    failure = classify_preseq_failure(log_text(log_path))
    tmp_output.unlink(missing_ok=True)
    if failure:
        status, message = failure
        output_txt.write_text(PRESEQ_HEADER, encoding="utf-8")
        write_status(status_tsv, sample, status, message)
        with open(log_path, "a", encoding="utf-8") as log_handle:
            log_handle.write(
                "\nWARNING: preseq lc_extrap could not estimate this sample; "
                "wrote a header-only output so downstream QC can continue.\n"
            )
        return 0

    write_status(status_tsv, sample, "FAILED", f"preseq exited with code {proc.returncode}")
    return proc.returncode


sys.exit(run_preseq_lc_extrap(snakemake))
