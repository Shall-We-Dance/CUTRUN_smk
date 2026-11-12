#!/usr/bin/env bash
# Shared local/Actions entrypoint. Missing dependencies must fail, not skip.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
case "${1:-}" in
  core)
    command -v samtools
    command -v snakemake
    python -c 'import pyBigWig, yaml'
    python -m compileall -q workflow/rules/scripts tests
    python tests/test_consensus.py
    python tests/test_fragment_coverage.py
    python tests/test_frip.py
    python tests/test_te_quant.py
    python tests/test_te_quant_pe.py
    python tests/test_workflow_config.py
    ;;
  multiqc)
    command -v multiqc
    python -c 'import multiqc'
    python tests/test_multiqc_summaries.py
    ;;
  idr)
    command -v idr
    python tests/test_idr.py
    ;;
  *)
    echo 'Usage: bash tests/run_ci.sh {core|multiqc|idr}' >&2
    exit 2
    ;;
esac
