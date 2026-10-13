#!/bin/bash
# Central configuration: paths, parameters, conda activation.
# Source from scripts in src/:  source "$(dirname "$0")/../config.sh"

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$PROJECT_DIR/data"
WORK_DIR="$PROJECT_DIR/work"
SHARE_DIR="$PROJECT_DIR/share"
LOGS_DIR="$PROJECT_DIR/logs"
UTILS_DIR="$PROJECT_DIR/utils"

# Inputs (sealed; processed by tools, never modified)
PCR_DIR="$DATA_DIR/pcr"
GENOME_DIR="$DATA_DIR/genomes_extracted"
REFERENCE_FASTA="$DATA_DIR/reference/hyalella_cox1_reference.fa"  # QC only

# Compute
THREADS="$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)"
CONDA_ENV="hyalella-cox1-v1"

# Step 1: public download
INGROUP_TAXON="Hyalella"
OUTGROUP_CANDIDATES=40   # over-fetch; QC and step 4 cut this to OUTGROUP_N
SLEN_MIN=300
SLEN_MAX=2000            # excludes whole mitogenomes (traded for speed)

# Step 3: QC
MIN_LENGTH=300
MIN_REF_IDENTITY=0.70

# Step 4: tree subset (speed + readability)
OUTGROUP_N=10
NEIGHBOURS_PER_SAMPLE=1  # extra top hits kept per sample, besides the reported nearest relative
MAX_PER_SPECIES=1        # backbone: public seqs kept per formally named species

# Step 6: trimming
MAX_GAP_FRACTION=0.5

activate_conda() {
    # conda's activate scripts reference unset variables
    set +u
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
    set -u
}

# Step 4: reading the % identity to the nearest public relative.
# Rule-of-thumb cox1 barcode cut-offs, not formal species delimitation.
SAME_LINEAGE_IDENTITY=97   # >= this: matches a lineage already in the databases
CLOSE_RELATIVE_IDENTITY=90 # between: related to, but distinct from, known lineages
MIN_HIT_OVERLAP=250        # bp aligned for a public hit to count

# Step 8: plot
CLUSTER_MAX_DEPTH=0.05   # subs/site; tight clades holding our samples get a labelled bracket
SUPPORT_SHOWN=0.9        # SH-like support at or above this is drawn as a dot
