#!/bin/bash
# Step 9: build the self-contained projector page share/summary.html.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda

python3 "$(dirname "$0")/09_summary.py" \
    "$SHARE_DIR/table/08_lineage_groups.tsv" \
    "$SHARE_DIR/table/04_nearest_public_relatives.tsv" \
    "$SHARE_DIR/figure/08_tree_cox1.png" \
    "$WORK_DIR/01_download/01_download_public_metadata.tsv" \
    "$SHARE_DIR/summary.html"
