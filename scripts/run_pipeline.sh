#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${SCRIPT_DIR}/run_rnaseq_core.R" "$@"
Rscript "${SCRIPT_DIR}/run_rnaseq_extension.R" "$@"
