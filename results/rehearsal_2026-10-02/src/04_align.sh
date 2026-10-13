#!/bin/bash
# Step 04: align the QC'd tree set with MAFFT (L-INS-i).
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda

OUTDIR="$WORK_DIR/04_align"
mkdir -p "$OUTDIR"
touch "$OUTDIR/.claude-sealed"
IN="$WORK_DIR/03_qc/03_qc_cox1.fasta"
OUT="$OUTDIR/04_align_cox1.fasta"

echo "[04] aligning $(grep -c '^>' "$IN") sequences with mafft L-INS-i, $THREADS threads"
mafft --localpair --maxiterate 1000 --thread "$THREADS" --quiet "$IN" > "$OUT"
"$UTILS_DIR/normalize_fasta.py" "$OUT"
echo "[04] alignment: $(grep -c '^>' "$OUT") sequences, $(awk 'NR==2{print length($0)}' "$OUT") columns -> $OUT"
