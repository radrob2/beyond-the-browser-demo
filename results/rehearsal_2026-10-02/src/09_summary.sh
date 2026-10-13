#!/bin/bash
# Step 09: build the projector-friendly share/summary.html.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda
python "$SRC_DIR/09_summary.py"
