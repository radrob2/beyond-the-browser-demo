#!/bin/bash
# Step 6: drop alignment columns where more than MAX_GAP_FRACTION of taxa
# have a gap. With mixed-length cox1 this is what trims the ragged ends; the
# same rule is applied to interior columns. Sequences left with too little
# overlap to place reliably are dropped and reported.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda

OUTDIR="$WORK_DIR/06_trim"
mkdir -p "$OUTDIR"
touch "$OUTDIR/.claude-sealed"
IN="$WORK_DIR/05_align/05_align_cox1.fasta"
OUT="$OUTDIR/06_trim_cox1.fasta"
MIN_OVERLAP=200

# trimal -gt is the minimum fraction of sequences WITHOUT a gap
GT="$(python3 -c "print(1 - $MAX_GAP_FRACTION)")"
echo "[06] trimal -gt $GT (drop columns with > $MAX_GAP_FRACTION gaps)"
trimal -in "$IN" -out "$OUTDIR/06_trim_raw.fasta" -fasta -gt "$GT"

python3 - "$OUTDIR/06_trim_raw.fasta" "$OUT" "$MIN_OVERLAP" "$(dirname "$0")" <<'PY'
import sys
sys.path.insert(0, sys.argv[4])
from cox1_common import read_fasta, write_fasta, seq_id, classify
min_overlap = int(sys.argv[3])
recs = [(seq_id(h), s) for h, s in read_fasta(sys.argv[1])]
keep = [(h, s) for h, s in recs if len(s) - s.count('-') >= min_overlap]
dropped = [h for h, s in recs if len(s) - s.count('-') < min_overlap]
own_dropped = [h for h in dropped if classify(h) != 'public']
write_fasta(sys.argv[2], keep)
print(f"[06] trimmed alignment: {len(keep)} sequences x {len(keep[0][1])} columns")
print(f"[06] dropped for < {min_overlap} bp left after trimming: {len(dropped)} {dropped if dropped else ''}")
if own_dropped:
    print(f"[06] WARNING: own samples dropped: {own_dropped}")
PY
rm "$OUTDIR/06_trim_raw.fasta"
echo "[06] wrote $OUT"
