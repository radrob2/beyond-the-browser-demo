#!/bin/bash
# Step 02: merge public cox1 (Hyalella + outgroups) with our PCR and
# genome-extracted cox1. The QC reference is deliberately NOT merged.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda

OUTDIR="$WORK_DIR/02_merge"
mkdir -p "$OUTDIR"
touch "$OUTDIR/.claude-sealed"   # holds our sequences
OUT="$OUTDIR/02_merge_cox1.fasta"

cat "$WORK_DIR/01_download/01_download_hyalella_cox1_ncbi.fasta" \
    "$WORK_DIR/01_download/01_download_outgroup_cox1_ncbi.fasta" \
    "$WORK_DIR/01_download/01_download_hyalella_coi_bold.fasta" \
    "$PCR_DIR"/*.fa \
    "$GENOMES_DIR"/*.fa > "$OUT.tmp"
"$UTILS_DIR/normalize_fasta.py" "$OUT.tmp" "$OUT"
rm "$OUT.tmp"

n_pub=$(grep -c '^>' "$WORK_DIR/01_download/01_download_hyalella_cox1_ncbi.fasta")
n_bold=$(grep -c "^>" "$WORK_DIR/01_download/01_download_hyalella_coi_bold.fasta")
n_out=$(grep -c '^>' "$WORK_DIR/01_download/01_download_outgroup_cox1_ncbi.fasta")
n_pcr=$(cat "$PCR_DIR"/*.fa | grep -c '^>')
n_gen=$(cat "$GENOMES_DIR"/*.fa | grep -c '^>')
echo "[02] public Hyalella GenBank: $n_pub  BOLD-only: $n_bold  outgroups: $n_out  pcr: $n_pcr  genome-extracted: $n_gen"
echo "[02] merged: $(grep -c '^>' "$OUT") records -> $OUT"
