#!/bin/bash
# Step 01: download public Hyalella + outgroup cox1 from NCBI (with metadata).
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda

OUTDIR="$WORK_DIR/01_download"
mkdir -p "$OUTDIR"
echo "[01] downloading public cox1 -> $OUTDIR"
python "$SRC_DIR/01_download.py"
"$UTILS_DIR/normalize_fasta.py" "$OUTDIR/01_download_hyalella_cox1_ncbi.fasta"
"$UTILS_DIR/normalize_fasta.py" "$OUTDIR/01_download_outgroup_cox1_ncbi.fasta"
"$UTILS_DIR/normalize_fasta.py" "$OUTDIR/01_download_hyalella_coi_bold.fasta"
echo "[01] done"
