#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

choose_python() {
  for exe in python3.12 python3.11 python3.10 python3 python; do
    if command -v "$exe" >/dev/null 2>&1; then
      if "$exe" -c 'import sys; raise SystemExit(0 if (3,10) <= sys.version_info[:2] < (3,13) else 1)' >/dev/null 2>&1; then
        echo "$exe"
        return 0
      fi
    fi
  done
  echo "ERROR: Python 3.10-3.12 is required. Python 3.12 is recommended." >&2
  exit 2
}

PY="$(choose_python)"
echo "Using $PY: $($PY --version 2>&1)"
"$PY" -m venv .venv_ag3
if [[ -x .venv_ag3/bin/python ]]; then
  VPY=.venv_ag3/bin/python
elif [[ -x .venv_ag3/Scripts/python.exe ]]; then
  VPY=.venv_ag3/Scripts/python.exe
else
  echo "ERROR: virtual environment was created but its Python executable was not found." >&2
  exit 2
fi
"$VPY" -m pip install --upgrade pip setuptools wheel
"$VPY" -m pip install -r scripts/python_requirements_ag3.txt
"$VPY" - <<'PY'
import sys, malariagen_data, pandas, numpy, yaml
print("Python:", sys.version.split()[0])
print("malariagen_data:", getattr(malariagen_data, "__version__", "unknown"))
print("pandas:", pandas.__version__)
print("numpy:", numpy.__version__)
print("Targeted Ag3 environment ready.")
PY
mkdir -p metadata
"$VPY" -m pip freeze > metadata/python_environment_freeze.txt
echo "Environment freeze written to metadata/python_environment_freeze.txt"
