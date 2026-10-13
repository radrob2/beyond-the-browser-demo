#!/bin/bash
# Step 08: table of nearest public relatives of our samples (share/).
set -euo pipefail
source "$(dirname "$0")/../config.sh"
activate_conda
python "$SRC_DIR/08_table.py"
