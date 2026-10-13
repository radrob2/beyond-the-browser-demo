#!/bin/bash
# Step 05: trim gappy alignment ends (> MAX_GAP_FRACTION of taxa gapped).
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda

OUTDIR="$WORK_DIR/05_trim"
mkdir -p "$OUTDIR"
touch "$OUTDIR/.claude-sealed"
python "$SRC_DIR/05_trim.py"
echo "[05] done -> $OUTDIR/05_trim_cox1.fasta"
