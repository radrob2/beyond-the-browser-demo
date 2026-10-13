#!/bin/bash
# Step 06: maximum-likelihood tree with IQ-TREE 3
# (ModelFinder, 1000 ultrafast bootstraps + 1000 SH-aLRT replicates).
# Topological constraint: Hyalella is monophyletic with respect to the
# outgroup genera (relationships inside Hyalella and among outgroups are
# free). Without it, saturated cox1 pulls long-branch outgroups into Hyalella.
# Rooting on the outgroups happens at plot time (step 07).
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda

OUTDIR="$WORK_DIR/06_tree"
mkdir -p "$OUTDIR"
IN="$WORK_DIR/05_trim/05_trim_cox1.fasta"
PREFIX="$OUTDIR/06_tree_cox1"
CONSTRAINT="$OUTDIR/06_tree_constraint.nwk"

# constraint: ((all ingroup), outgroup1, outgroup2, ...);
awk -F'\t' 'NR > 1 { if ($2 == "outgroup") out = out "," $1; else ing = ing (ing ? "," : "") $1 }
    END { print "((" ing ")" out ");" }' "$WORK_DIR/03_qc/03_qc_tree_set.tsv" > "$CONSTRAINT"
echo "[06] constraint: Hyalella + our samples monophyletic vs $(awk -F'\t' '$2 == "outgroup"' "$WORK_DIR/03_qc/03_qc_tree_set.tsv" | wc -l) outgroups"

echo "[06] IQ-TREE on $(grep -c '^>' "$IN") sequences, $THREADS threads"
iqtree3 -s "$IN" --prefix "$PREFIX" -m MFP -g "$CONSTRAINT" -B 1000 -alrt 1000 \
    -T "$THREADS" --seed 12345 -redo --quiet
# IQ-TREE writes a sequence copy when it collapses identical sequences
rm -f "$PREFIX.uniqueseq.phy"
echo "[06] best model: $(grep 'Best-fit model' "$PREFIX.iqtree" | sed 's/.*: //')"
echo "[06] log-likelihood: $(grep 'Log-likelihood of the tree' "$PREFIX.iqtree" | awk '{print $5}')"
echo "[06] done -> $PREFIX.treefile"
