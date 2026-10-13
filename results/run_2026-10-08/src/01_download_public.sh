#!/bin/bash
# Step 1: download public cox1 sequences (all Hyalella + outgroup candidates)
# and their GenBank source metadata from NCBI via direct E-utilities HTTP.
#
# Usage: ./src/01_download_public.sh 2>&1 | tee logs/01_download_public_$(date +%Y%m%d_%H%M%S).log
set -euo pipefail
source "$(dirname "$0")/../config.sh"

OUTDIR="$WORK_DIR/01_download"
mkdir -p "$OUTDIR"

START=$SECONDS
echo "=== 01_download_public: $(date '+%Y-%m-%d %H:%M:%S') ==="
echo "ingroup=$INGROUP_TAXON  outgroup_candidates=$OUTGROUP_CANDIDATES  length=${SLEN_MIN}..${SLEN_MAX} bp"

# stdlib-only python; no conda env needed for this step
python3 "$(dirname "$0")/01_download_public.py" \
    --outdir "$OUTDIR" \
    --ingroup "$INGROUP_TAXON" \
    --outgroup-n "$OUTGROUP_CANDIDATES" \
    --slen-min "$SLEN_MIN" \
    --slen-max "$SLEN_MAX"

echo ""
echo "normalizing fasta (single-line, uppercase)"
for f in "$OUTDIR/01_download_hyalella_cox1.fasta" "$OUTDIR/01_download_outgroup_cox1.fasta"; do
    python3 "$UTILS_DIR/normalize_fasta.py" "$f"
    echo "  $(grep -c '^>' "$f") records  $f"
done
echo "  $(($(wc -l < "$OUTDIR/01_download_public_metadata.tsv") - 1)) rows     $OUTDIR/01_download_public_metadata.tsv"
echo "=== done in $((SECONDS - START))s ==="
