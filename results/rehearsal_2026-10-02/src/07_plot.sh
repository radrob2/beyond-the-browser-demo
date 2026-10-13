#!/bin/bash
# Step 07: plot the rooted tree (PDF + PNG) into share/figure/.
# Needs the step 08 table (group brackets); run_pipeline.sh runs 08 first.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda
[ -s "$SHARE_DIR/08_nearest_relatives.tsv" ] || { echo "[07] run step 08 first" >&2; exit 1; }
python "$SRC_DIR/07_plot.py"
