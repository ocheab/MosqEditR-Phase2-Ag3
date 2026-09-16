#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${MOSQEDIT_PYTHON:-python3}
"$PY" scripts/acquire_phase2_discovery.py "$@"
if [[ " $* " == *" --dry-run "* ]]; then
  exit 0
fi
Rscript --vanilla R/00b_acquire_phase2_phased_remote.R
