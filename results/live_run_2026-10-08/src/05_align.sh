#!/bin/bash
# Step 5: multiple sequence alignment of the tree subset with MAFFT.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda

OUTDIR="$WORK_DIR/05_align"
mkdir -p "$OUTDIR"
touch "$OUTDIR/.claude-sealed"
IN="$WORK_DIR/04_subset/04_subset_cox1.fasta"
OUT="$OUTDIR/05_align_cox1.fasta"

echo "[05] mafft --auto on $(grep -c '>' "$IN") sequences, $THREADS threads"
mafft --auto --preservecase --thread "$THREADS" "$IN" > "$OUT" 2> "$OUTDIR/05_align_mafft.log"
"$UTILS_DIR/normalize_fasta.py" "$OUT"
grep -E '^Strategy|^ [A-Z]-INS|FFT-NS' "$OUTDIR/05_align_mafft.log" | tail -2 || true
echo "[05] aligned: $(grep -c '>' "$OUT") sequences x $(sed -n 2p "$OUT" | tr -d '\n' | wc -c | tr -d ' ') columns"
echo "[05] wrote $OUT"
