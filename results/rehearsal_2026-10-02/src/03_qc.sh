#!/bin/bash
# Step 03: QC (length, identity to reference, header uniqueness, dedup),
# nearest-public-relative search, and selection of the tree set.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda

OUTDIR="$WORK_DIR/03_qc"
mkdir -p "$OUTDIR"
python "$SRC_DIR/03_qc.py"
"$UTILS_DIR/normalize_fasta.py" "$OUTDIR/03_qc_cox1_all.fasta"
"$UTILS_DIR/normalize_fasta.py" "$OUTDIR/03_qc_cox1.fasta"
echo "[03] done"
