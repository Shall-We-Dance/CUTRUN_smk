#!/usr/bin/env python3
"""Regression tests for CUT&RUN and FRiP MultiQC custom-content payloads."""

import json
import runpy
import shutil
import subprocess
import tempfile
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow" / "rules" / "scripts"


def run_cutrun_summary(tmpdir):
    usable = tmpdir / "sample.usable_fragments.tsv"
    usable.write_text(
        "metric\tvalue\n"
        "sample\tSampleA\n"
        "read_type\tPE\n"
        "total_reads_after_processing\t100\n"
        "mapped_reads_after_processing\t90\n"
        "usable_reads_after_processing\t80\n"
        "proper_pair_reads_after_processing\t60\n"
        "usable_fragments_after_processing\t30\n"
    )
    fragments = tmpdir / "SampleA.fragment_lengths.tsv"
    fragments.write_text("fragment_length\tcount\n100\t3\n150\t2\n")
    outputs = SimpleNamespace(
        tsv=str(tmpdir / "cutrun_qc_summary.tsv"),
        fragments_mqc=str(tmpdir / "fragment_lengths_mqc.json"),
        summary_mqc=str(tmpdir / "cutrun_qc_summary_mqc.json"),
    )
    snakemake = SimpleNamespace(
        input=SimpleNamespace(usable=[str(usable)], fragments=[str(fragments)], spikein=[]),
        output=outputs,
    )
    runpy.run_path(
        str(SCRIPTS / "summarize_cutrun_qc.py"),
        init_globals={"snakemake": snakemake},
    )
    summary = json.loads(Path(outputs.summary_mqc).read_text())
    fragments_payload = json.loads(Path(outputs.fragments_mqc).read_text())
    assert summary["pconfig"]["id"] == "cutrun_qc_table"
    assert summary["data"]["SampleA"]["mapped_fraction"] == 0.9
    assert summary["data"]["SampleA"]["usable_fragments"] == 30
    assert fragments_payload["pconfig"]["id"] == "cutrun_fragment_lengths_plot"
    assert fragments_payload["data"]["SampleA"] == {"100": 3, "150": 2}
    return outputs


def write_frip(path, usable, in_peaks):
    path.write_text(
        "metric\tvalue\n"
        "total_reads\t200\n"
        "mapped_reads\t180\n"
        f"usable_templates\t{usable}\n"
        f"templates_in_peaks\t{in_peaks}\n"
        f"fragment_frip\t{in_peaks / usable:.6f}\n"
        "reads_in_peaks\t80\n"
        "read_frip\t0.444444\n"
    )


def run_frip_summary(tmpdir):
    replicate = tmpdir / "SampleA.rep1.frip.tsv"
    pooled = tmpdir / "SampleA.pooled.frip.tsv"
    write_frip(replicate, 100, 50)
    write_frip(pooled, 200, 120)
    outputs = SimpleNamespace(
        tsv=str(tmpdir / "frip_summary.tsv"),
        mqc=str(tmpdir / "frip_summary_mqc.json"),
    )
    snakemake = SimpleNamespace(
        input=SimpleNamespace(per_rep=[str(replicate)], pooled=[str(pooled)]),
        output=outputs,
    )
    runpy.run_path(
        str(SCRIPTS / "summarize_frip.py"),
        init_globals={"snakemake": snakemake},
    )
    payload = json.loads(Path(outputs.mqc).read_text())
    assert payload["pconfig"]["id"] == "frip_table"
    assert payload["data"]["SampleA.rep1"]["fragment_frip"] == 0.5
    assert payload["data"]["SampleA (pooled)"]["fragment_frip"] == 0.6
    return outputs


def run_multiqc_if_available(tmpdir, cutrun_outputs, frip_outputs):
    multiqc = shutil.which("multiqc")
    if not multiqc:
        return False
    report_dir = tmpdir / "multiqc"
    subprocess.run(
        [
            multiqc,
            "--force",
            "-o",
            str(report_dir),
            cutrun_outputs.summary_mqc,
            cutrun_outputs.fragments_mqc,
            frip_outputs.mqc,
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    report_data = "\n".join(
        path.read_text(errors="replace")
        for path in (report_dir / "multiqc_data").iterdir()
        if path.is_file()
    )
    for section_id in ("cutrun_qc_table", "cutrun_fragment_lengths_plot", "frip_table"):
        assert section_id in report_data, f"MultiQC dropped custom section {section_id}"
    return True


def main():
    with tempfile.TemporaryDirectory(prefix="multiqc_summary_test_") as tmp:
        tmpdir = Path(tmp)
        cutrun_outputs = run_cutrun_summary(tmpdir)
        frip_outputs = run_frip_summary(tmpdir)
        rendered = run_multiqc_if_available(tmpdir, cutrun_outputs, frip_outputs)
    suffix = " and render in MultiQC" if rendered else ""
    print(f"PASS: MultiQC custom-content summaries have valid keyed data{suffix}")


if __name__ == "__main__":
    main()
