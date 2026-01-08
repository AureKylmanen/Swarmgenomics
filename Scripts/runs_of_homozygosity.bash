#!/bin/bash
# Script Name: runs_of_homozygosity.bash
# Description: Analysing and plotting runs of homozygosity
# Author: Aure Kylmänen
# Created: 08.01.2026
# Current working version
# ====================================
set -euo pipefail

# ============================
# USAGE
# ============================
# ./roh.sh <species> <params.txt>

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <species> <params.txt>"
    exit 1
fi

SPECIES="$1"
PARAMS_FILE="$2"

# ============================
# LOAD PARAMETERS
# ============================

if [[ ! -f "$PARAMS_FILE" ]]; then
    echo "ERROR: Parameter file not found: $PARAMS_FILE"
    exit 1
fi

source "$PARAMS_FILE"

# ============================
# REQUIRED PARAMETERS CHECK
# ============================

: "${WORKING_DIR:?Missing WORKING_DIR in params file}"
: "${ROH_MIN_GQ:?Missing ROH_MIN_GQ in params file}"
: "${ROH_DEFAULT_AF:?Missing ROH_DEFAULT_AF in params file}"
: "${RoH_PLOT:?Missing RoH_PLOT in params file}"

# ============================
# DIRECTORIES
# ============================


# Create results directory if needed
mkdir -p "$RESULTS_DIR"

# ============================
# VALIDATION
# ============================

if [[ ! -d "$VCF_DIR" ]]; then
    echo "ERROR: VCF directory not found: $VCF_DIR"
    echo "Have you run variant calling?"
    exit 1
fi

shopt -s nullglob
VCF_FILES=("${VCF_DIR}"/*.vcf.gz)
shopt -u nullglob

if [[ ${#VCF_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No VCF files found in $VCF_DIR"
    exit 1
fi

# ============================
# RUN ROH
# ============================

echo "======================================"
echo "Running Runs of Homozygosity (RoH)"
echo "Species: $SPECIES"
echo "Min genotype quality: $ROH_MIN_GQ"
echo "Default allele frequency: $ROH_DEFAULT_AF"
echo "======================================"

echo "Merging VCF files..."
bcftools merge --force-samples -O z \
    -o "${VCF_DIR}/merged.vcf.gz" \
    "${VCF_FILES[@]}"

echo "Running bcftools roh..."
bcftools roh \
    -G "$ROH_MIN_GQ" \
    --AF-dflt "$ROH_DEFAULT_AF" \
    "${VCF_DIR}/merged.vcf.gz" \
    > "${VCF_DIR}/roh_results.txt"

# ============================
# PROCESS OUTPUT
# ============================

echo "Extracting RG block from RoH output..."
awk '/^# RG/ {print; in_block=1} in_block && /^RG/' \
    "${VCF_DIR}/roh_results.txt" \
    > "${VCF_DIR}/RG.txt"

# ============================
# PLOTTING
# ============================

echo "Generating RoH plots..."
cd "$VCF_DIR"

set +e
Rscript "$RoH_PLOT"
plot_status=$?
set -e

if [[ $plot_status -ne 0 ]]; then
    echo "⚠️  RoH plotting failed, continuing without plots."
fi

# ============================
# COPY RESULTS
# ============================

echo "Copying results to ${RESULTS_DIR}..."
cp roh_results.txt RG.txt roh_bar_plots.png \
    "$RESULTS_DIR/" 2>/dev/null || \
    echo "⚠️  Some result files could not be copied."

echo "RoH analysis completed successfully."

