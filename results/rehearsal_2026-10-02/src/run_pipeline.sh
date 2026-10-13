#!/bin/bash
# Run the whole pipeline, or a range of steps: ./src/run_pipeline.sh [FIRST [LAST]]
# Each step logs to logs/{script}_{YYYYMMDD_HHMMSS}.log
set -euo pipefail
source "$(dirname "$0")/../config.sh"

FIRST="${1:-1}"
LAST="${2:-9}"
# 08_table runs before 07_plot: the plot labels groups using the table
STEPS=(01_download 02_merge 03_qc 04_align 05_trim 06_tree 08_table 07_plot 09_summary)

mkdir -p "$LOGS_DIR"
for i in $(seq "$FIRST" "$LAST"); do
    step="${STEPS[$((i - 1))]}"
    log="$LOGS_DIR/${step}_$(date +%Y%m%d_%H%M%S).log"
    echo "=== step $step (log: $log)"
    "$SRC_DIR/$step.sh" 2>&1 | tee "$log"
done
