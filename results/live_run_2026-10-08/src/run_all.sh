#!/bin/bash
# Run the whole pipeline in order, logging each step to logs/.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

STAMP="$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOGS_DIR"
START=$SECONDS
for step in 01_download_public 02_merge 03_qc 04_nearest_subset 05_align 06_trim 07_tree 08_plot 09_summary; do
    t0=$SECONDS
    "$PROJECT_DIR/src/$step.sh" 2>&1 | tee "$LOGS_DIR/${step}_${STAMP}.log"
    echo "[run_all] $step: $((SECONDS - t0)) s"
done
echo "[run_all] total: $((SECONDS - START)) s"
