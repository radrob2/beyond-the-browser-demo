#!/bin/bash
# Step 4: (a) nearest public relative of every sample of ours, searched
# against ALL public sequences that passed QC -> share/table/;
# (b) pick the subset of sequences that goes into the alignment and tree.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda

OUTDIR="$WORK_DIR/04_subset"
mkdir -p "$OUTDIR" "$SHARE_DIR/table"
touch "$OUTDIR/.claude-sealed"
SRC="$(dirname "$0")"

python3 "$SRC/04_nearest_subset.py" split "$WORK_DIR" "$SHARE_DIR"

echo "[04] vsearch: every sample of ours vs every QC-passed public sequence"
vsearch --usearch_global "$OUTDIR/04_own_all.fasta" \
    --db "$OUTDIR/04_public_all.fasta" \
    --id 0.5 --iddef 2 --strand plus \
    --maxaccepts 0 --maxrejects 0 --maxhits 200 \
    --threads "$THREADS" \
    --userout "$OUTDIR/04_nearest_hits.tsv" \
    --userfields query+target+id+alnlen+mism+gaps \
    --quiet

python3 "$SRC/04_nearest_subset.py" finish "$WORK_DIR" "$SHARE_DIR" \
    "$SAME_LINEAGE_IDENTITY" "$CLOSE_RELATIVE_IDENTITY" "$MIN_HIT_OVERLAP" \
    "$NEIGHBOURS_PER_SAMPLE" "$MAX_PER_SPECIES" "$OUTGROUP_N"
echo "[04] wrote $SHARE_DIR/table/04_nearest_public_relatives.tsv"
echo "[04] wrote $OUTDIR/04_subset_cox1.fasta"
