#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
export MPLCONFIGDIR="$PWD/.mpl" NUMBA_CACHE_DIR="$PWD/.numba" CELLTYPIST_FOLDER="$PWD/data/celltypist" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
case "${1:-}" in
 checkpoint) .venv/bin/python scripts/02_checkpoint.py 2>&1 | tee logs/checkpoint.log;;
 metadata) .venv/bin/python scripts/01_metadata.py 2>&1 | tee logs/metadata.log;;
 audit) .venv/bin/python scripts/04_input_audit.py 2>&1 | tee logs/input_audit.log;;
 verify) .venv/bin/python scripts/05_verify_gate.py 2>&1 | tee logs/gate_verification.log;;
 report) .venv/bin/python scripts/06_gate_report.py 2>&1 | tee logs/gate_report.log;;
 test) .venv/bin/python -m pytest tests -q 2>&1 | tee logs/tests.log;;
 pilot) .venv/bin/python scripts/03_pilot.py 2>&1 | tee logs/pilot.log;;
 main) .venv/bin/python scripts/11_main.py;;
 main-summary) .venv/bin/python scripts/13_summarize.py;;
 secondary) .venv/bin/python scripts/14_scib_secondary.py;;
 *) echo 'Usage: bash run_all.sh {checkpoint|metadata|test|pilot|audit|verify|report|main}'; exit 2;;
esac
