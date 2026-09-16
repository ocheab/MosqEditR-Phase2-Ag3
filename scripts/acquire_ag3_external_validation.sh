#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
exec Rscript --vanilla R/05c_acquire_ag3_external_public_sanger.R "$@"
