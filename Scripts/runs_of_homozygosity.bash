#!/bin/bash
# Script Name: runs_of_homozygosity.bash
# Description: Running and plotting runs of homozygosity analysis
# Author: Aure Kylmänen
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
    echo "❌ ERROR: Parameter file not found: $PARAMS_FILE"
    exit 1
fi

source "$PARAMS_FILE"


# ============================
# REQUIRED PARAMETERS CHECK
# ============================
: "${WORKING_DIR:?Missing WORKING_DIR}"
: "${VCF_DIR:?Missing VCF_DIR}"
: "${RESULTS_DIR:?Missing RESULTS_DIR}"
: "${BCFTOOLS:?Missing BCFTOOLS}"
: "${ROH_MIN_GQ:?Missing ROH_MIN_GQ}"
: "${ROH_DEFAULT_AF:?Missing ROH_DEFAULT_AF}"
: "${ROH_PLOT_SCRIPT:?Missing ROH_PLOT_SCRIPT}"
: "${ROH_COUNT_COLOR:?Missing ROH_COUNT_COLOR}"
: "${ROH_LENGTH_COLOR:?Missing ROH_LENGTH_COLOR}"
: "${ROH_BINS:?Missing ROH_BINS}"
: "${ROH_BIN_LABELS:?Missing ROH_BIN_LABELS}"

mkdir -p "$RESULTS_DIR"


# ============================
# VALIDATION
# ============================

if [[ ! -d "$VCF_DIR" ]]; then
    echo "❌ ERROR: VCF directory not found: $VCF_DIR"
    echo "Have you run variant calling?"
    exit 1
fi

shopt -s nullglob
VCF_FILES=("${VCF_DIR}"/*.vcf.gz)
shopt -u nullglob

if [[ ${#VCF_FILES[@]} -eq 0 ]]; then
    echo "❌ ERROR: No VCF files found in $VCF_DIR"
    exit 1
fi

# ============================
# RUN ROH
# ============================
echo "Merging vcf files"
"${BCFTOOLS}/bcftools" merge --force-samples -O z \
  -o "${VCF_DIR}/merged.vcf.gz" \
  "${VCF_FILES[@]}"

echo "Starting RoH analysis"
"${BCFTOOLS}/bcftools" roh \
  -G "$ROH_MIN_GQ" \
  --AF-dflt "$ROH_DEFAULT_AF" \
  "${VCF_DIR}/merged.vcf.gz" \
  > "${VCF_DIR}/roh_results.txt"

# ============================
# EXTRACT RG BLOCK
# ============================
awk '/^# RG/ {in_block=1} in_block && /^RG/' \
  "${VCF_DIR}/roh_results.txt" \
  > "${VCF_DIR}/RG.txt"

# ============================
# PLOTTING
# ============================
echo "Generating RoH plots..."

Rscript "$ROH_PLOT_SCRIPT" \
  "${VCF_DIR}/RG.txt" \
  "$ROH_BINS" \
  "$ROH_BIN_LABELS" \
  "$ROH_COUNT_COLOR" \
  "$ROH_LENGTH_COLOR" \
  "${RESULTS_DIR}/roh_bar_plots.png" || \
  echo "⚠️ RoH plotting failed."

# ============================
# COPY RESULTS
# ============================
cp -v "${VCF_DIR}/roh_results.txt" "${VCF_DIR}/RG.txt" "${RESULTS_DIR}/" 2>/dev/null || true

echo "✅ RoH analysis completed successfully."
