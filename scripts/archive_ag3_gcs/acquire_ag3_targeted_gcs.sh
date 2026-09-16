#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ -x .venv_ag3/bin/python ]]; then
  PY=.venv_ag3/bin/python
elif [[ -x .venv_ag3/Scripts/python.exe ]]; then
  PY=.venv_ag3/Scripts/python.exe
else
  echo "ERROR: .venv_ag3 not found. First run: bash scripts/setup_ag3_python.sh" >&2
  exit 2
fi
exec "$PY" scripts/acquire_ag3_targeted.py "$@"
