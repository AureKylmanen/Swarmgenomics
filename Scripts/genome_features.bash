#!/bin/bash
# Script Name: genome_features.bash
# Description: Generate genome feature stats from BAM (idxstats), estimate read length,
#              and produce plots with R.
# Author: Aure Kylmänen
# ====================================

set -euo pipefail

# ============================
# USAGE
# ============================
# ./genome_features.bash <species> <params.txt>

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
# REQUIRED PARAMETERS
# ============================
: "${WORKING_DIR:?Missing WORKING_DIR}"
: "${BAM_FILE:?Missing BAM_FILE}"
: "${IDXSTATS_INPUT:?Missing IDXSTATS_INPUT}"
: "${IDXSTATS_CSV:?Missing IDXSTATS_CSV}"
: "${idxstats_PLOT:?Missing idxstats_PLOT}"
: "${RESULTS_DIR:?Missing RESULTS_DIR}"
: "${IDXSTATS_TOP_SCAFFOLDS:?Missing IDXSTATS_TOP_SCAFFOLDS}"
: "${IDXSTATS_COLOR_PALETTE:?Missing IDXSTATS_COLOR_PALETTE}"

mkdir -p "$RESULTS_DIR"

# ============================
# GENERATE IDXSTATS
# ============================
echo "Running samtools idxstats on BAM..."
samtools idxstats "$BAM_FILE" > "$IDXSTATS_INPUT"

if [[ ! -s "$IDXSTATS_INPUT" ]]; then
    echo "❌ ERROR: idxstats output file is empty or missing: $IDXSTATS_INPUT"
    exit 1
fi

# ============================
# ESTIMATE READ LENGTH
# ============================
echo "Estimating read length from BAM..."
READ_LENGTH=$(samtools view "$BAM_FILE" | head -n 1000 | \
    awk '{if($10 != "*") print length($10)}' | sort | uniq -c | sort -nr | head -n1 | awk '{print $2}' || true)

if [[ -z "$READ_LENGTH" ]]; then
    echo "❌ ERROR: Failed to detect read length from BAM file"
    exit 1
fi

echo "✅ Detected read length: $READ_LENGTH"

# ============================
# CREATE CLEANED CSV
# ============================
echo "Creating cleaned CSV file..."
awk 'BEGIN {OFS=","; print "chrom,length,mapped,unmapped"} {print $1,$2,$3,$4}' "$IDXSTATS_INPUT" > "$IDXSTATS_CSV"

if [[ ! -s "$IDXSTATS_CSV" ]]; then
    echo "❌ ERROR: Failed to create idxstats CSV: $IDXSTATS_CSV"
    exit 1
fi

# ============================
# RUN IDXSTATS PLOTTING
# ============================
echo "Running R script for idxstats plotting..."
Rscript "$idxstats_PLOT" \
  "$IDXSTATS_CSV" \
  "$READ_LENGTH" \
  "$IDXSTATS_TOP_SCAFFOLDS" \
  "$IDXSTATS_COLOR_PALETTE" \
  "${IDXSTATS_CUSTOM_COLORS:-}"

# ============================
# COPY GENERATED FILES TO RESULTS
# ============================
echo "Copying generated PNG plots to results directory..."

# Check for PNGs and move/rename them by species
if compgen -G "*.png" > /dev/null; then
    for f in idxstats_scaffold_lengths.png idxstats_unmapped_reads_per_mbp.png idxstats_coverage.png; do
        if [[ -f "$f" ]]; then
            new_name="${SPECIES}_$f"
            mv -v "$f" "$RESULTS_DIR/$new_name"
        fi
    done
else
    echo "⚠️ No PNG files found in current directory after plotting."
fi

# Copy CSV, input, and R script to results directory
echo "Copying CSV, input, and R script to results directory..."
cp -v "$IDXSTATS_CSV" "$IDXSTATS_INPUT" "$idxstats_PLOT" "$RESULTS_DIR/"

echo "✅ idxstats analysis completed. All results saved in $RESULTS_DIR."
