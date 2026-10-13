#!/bin/bash
# Step 2: merge public cox1 with our PCR and genome-extracted cox1.
# The QC reference is deliberately NOT merged.
set -euo pipefail
source "$(dirname "$0")/../config.sh"

OUTDIR="$WORK_DIR/02_merge"
mkdir -p "$OUTDIR"
OUT="$OUTDIR/02_merge_cox1.fasta"

echo "[02] merging public + pcr + genome-extracted cox1"
# awk adds a newline after each file so records never run together
awk 'FNR==1 && NR>1 {print ""} {print}' \
    "$WORK_DIR/01_download/01_download_hyalella_cox1.fasta" \
    "$WORK_DIR/01_download/01_download_outgroup_cox1.fasta" \
    "$PCR_DIR"/*.fa \
    "$GENOME_DIR"/*.fa \
    | grep -v '^[[:space:]]*$' > "$OUT"
"$UTILS_DIR/normalize_fasta.py" "$OUT"
touch "$OUTDIR/.claude-sealed"

python3 - "$OUT" "$(dirname "$0")" <<'PY'
import sys, collections
sys.path.insert(0, sys.argv[2])
from cox1_common import read_fasta, classify
recs = read_fasta(sys.argv[1])
counts = collections.Counter(classify(h) for h, _ in recs)
print(f"[02] merged records: {len(recs)}")
for k in ('public', 'pcr', 'genome', 'unknown'):
    print(f"[02]   {k}: {counts.get(k, 0)}")
if counts.get('unknown'):
    sys.exit("[02] ERROR: headers that match no source pattern")
PY
echo "[02] wrote $OUT"
