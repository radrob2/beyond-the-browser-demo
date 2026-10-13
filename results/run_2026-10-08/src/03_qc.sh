#!/bin/bash
# Step 3: QC. length filter -> identity to the Hyalella cox1 reference (also
# fixes strand) -> deduplicate identical sequences -> unique headers.
# The reference is used for the comparison only and never enters the dataset.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda

OUTDIR="$WORK_DIR/03_qc"
mkdir -p "$OUTDIR"
touch "$OUTDIR/.claude-sealed"
IN="$WORK_DIR/02_merge/02_merge_cox1.fasta"
SRC="$(dirname "$0")"

echo "[03] length filter (>= $MIN_LENGTH bp)"
python3 "$SRC/03_qc.py" length "$IN" "$OUTDIR" "$MIN_LENGTH"

echo "[03] identity to reference (>= $MIN_REF_IDENTITY), both strands"
vsearch --usearch_global "$OUTDIR/03_qc_length_cox1.fasta" \
    --db "$REFERENCE_FASTA" \
    --id "$MIN_REF_IDENTITY" --strand both \
    --maxaccepts 1 --maxrejects 0 \
    --threads "$THREADS" \
    --userout "$OUTDIR/03_qc_reference_hits.tsv" \
    --userfields query+id+qstrand+alnlen \
    --quiet

echo "[03] filter, orient, deduplicate, check headers"
python3 "$SRC/03_qc.py" finish "$IN" "$OUTDIR" "$REFERENCE_FASTA"
echo "[03] wrote $OUTDIR/03_qc_cox1.fasta"
