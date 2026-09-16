#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript is not on PATH." >&2
  exit 2
fi
if [[ -x .venv_ag3/bin/python ]]; then
  PY=.venv_ag3/bin/python
elif [[ -x .venv_ag3/Scripts/python.exe ]]; then
  PY=.venv_ag3/Scripts/python.exe
else
  echo "ERROR: Python environment missing. Run bash scripts/setup_ag3_python.sh first." >&2
  exit 2
fi

echo "[1/3] Running frozen-panel preflight..."
Rscript --vanilla R/00_preflight_and_freeze.R

echo "[2/3] Acquiring ONLY targeted Ag3 metadata/haplotypes/accessibility..."
"$PY" scripts/acquire_ag3_targeted.py

echo "[3/3] Running local-only R analysis from Step 01 onward..."
MOSQEDIT_START_STEP=01 Rscript --vanilla run_all.R
