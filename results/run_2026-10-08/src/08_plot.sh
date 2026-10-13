#!/bin/bash
# Step 8: plot the tree (PDF + PNG) and write the grouped lineage table.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda

mkdir -p "$SHARE_DIR/figure" "$SHARE_DIR/table"
python3 "$(dirname "$0")/08_plot.py" \
    "$WORK_DIR/07_tree/07_tree_cox1.treefile" \
    "$WORK_DIR/01_download/01_download_public_metadata.tsv" \
    "$SHARE_DIR/table/04_nearest_public_relatives.tsv" \
    "$SHARE_DIR/figure/08_tree_cox1" \
    "$SHARE_DIR/table/08_lineage_groups.tsv" \
    "$CLUSTER_MAX_DEPTH" "$SUPPORT_SHOWN"
