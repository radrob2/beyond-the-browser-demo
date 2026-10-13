#!/bin/bash
# Step 7: approximate maximum-likelihood tree with FastTree (GTR+Gamma,
# SH-like local supports). Chosen for speed: no model selection, no bootstrap.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda

OUTDIR="$WORK_DIR/07_tree"
mkdir -p "$OUTDIR"
IN="$WORK_DIR/06_trim/06_trim_cox1.fasta"
OUT="$OUTDIR/07_tree_cox1.treefile"

echo "[07] FastTreeMP -nt -gtr -gamma, $THREADS threads"
OMP_NUM_THREADS="$THREADS" FastTreeMP -nt -gtr -gamma -quiet \
    -log "$OUTDIR/07_tree_fasttree.log" "$IN" > "$OUT"
echo "[07] tips: $(tr ',' '\n' < "$OUT" | wc -l | tr -d ' ')"
echo "[07] wrote $OUT"
