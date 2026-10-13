#!/bin/bash
# Central configuration for the Hyalella cox1 phylogeny pipeline.
# Sourced by every script in src/.

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$PROJECT_DIR/data"
WORK_DIR="$PROJECT_DIR/work"
SHARE_DIR="$PROJECT_DIR/share"
LOGS_DIR="$PROJECT_DIR/logs"
SRC_DIR="$PROJECT_DIR/src"
UTILS_DIR="$PROJECT_DIR/utils"

THREADS=24
CONDA_ENV="hyalella_dryrun"

# input data (sealed)
PCR_DIR="$DATA_DIR/pcr"
GENOMES_DIR="$DATA_DIR/genomes_extracted"
REFERENCE_FA="$DATA_DIR/reference/hyalella_cox1_reference.fa"

# QC parameters
MIN_LEN=300
MIN_IDENTITY=70          # % identity to the reference
MAX_GAP_FRACTION=0.5     # trim alignment-end columns with more gaps than this

# public-data subsampling for the tree
# (the nearest-relative table always uses ALL QC-passed public sequences)
TOP_HITS_PER_SAMPLE=2    # public nearest relatives of each of our samples, always kept
CLUSTER_ID=0.95          # public sequences clustered at this identity; one representative each

# NCBI
NCBI_EMAIL="${NCBI_EMAIL:-}"
NCBI_API_KEY="${NCBI_API_KEY:-}"

activate_conda() {
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
}

export PROJECT_DIR DATA_DIR WORK_DIR SHARE_DIR LOGS_DIR SRC_DIR UTILS_DIR THREADS CONDA_ENV
export PCR_DIR GENOMES_DIR REFERENCE_FA MIN_LEN MIN_IDENTITY MAX_GAP_FRACTION
export CLUSTER_ID TOP_HITS_PER_SAMPLE NCBI_EMAIL NCBI_API_KEY
