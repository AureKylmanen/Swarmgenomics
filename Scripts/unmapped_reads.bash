#!/bin/bash
# Script Name: unmapped_reads.bash
# Description: Extract unmapped reads, classify with Kraken2, and generate Krona plots
# Author: Nikolas Vellnow
# Edited by: Aure Kylmänen
# ====================================

set -euo pipefail

# ============================
# USAGE
# ============================
# ./unmapped_reads.bash <species> <params.txt>
# Example: ./unmapped_reads.bash golden_eagle params.txt

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
# CHECK REQUIRED VARIABLES
# ============================
: "${WORKING_DIR:?Missing WORKING_DIR}"
: "${BAM_FILE:?Missing BAM_FILE}"
: "${THREADS:?Missing THREADS}"
: "${KRAKEN_DB:?Missing KRAKEN_DB}"

: "${UNMAPPED_BAM:?Missing UNMAPPED_BAM}"
: "${UNMAPPED_FASTQ:?Missing UNMAPPED_FASTQ}"
: "${KRAKEN_OUTPUT:?Missing KRAKEN_OUTPUT}"
: "${KRAKEN_REPORT:?Missing KRAKEN_REPORT}"
: "${CLASSIFIED_OUT:?Missing CLASSIFIED_OUT}"
: "${UNCLASSIFIED_OUT:?Missing UNCLASSIFIED_OUT}"
: "${KREPORT_TO_KRONA_SCRIPT:?Missing KREPORT_TO_KRONA_SCRIPT}"
: "${KRONA_HTML:?Missing KRONA_HTML}"
: "${RESULTS_DIR:?Missing RESULTS_DIR}"

# ============================
# SETUP DIRECTORIES
# ============================
cd "${WORKING_DIR}/${SPECIES}"
mkdir -p "${RESULTS_DIR}"

# ============================
# EXTRACT UNMAPPED READS
# ============================
echo "Extracting unmapped reads from BAM..."
samtools view -b -f 4 "$BAM_FILE" > "$UNMAPPED_BAM"

# ============================
# CONVERT BAM TO FASTQ
# ============================
echo "Converting BAM to FASTQ..."
bedtools bamtofastq -i "$UNMAPPED_BAM" -fq "$UNMAPPED_FASTQ"

# ============================
# RUN KRAKEN2
# ============================
echo "Classifying unmapped reads with Kraken2..."

source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate "$KRAKEN_CONDA_ENV"

kraken2 \
  --db "$KRAKEN_DB" \
  --output "$KRAKEN_OUTPUT" \
  --report "$KRAKEN_REPORT" \
  --use-names \
  --classified-out "$CLASSIFIED_OUT" \
  --unclassified-out "$UNCLASSIFIED_OUT" \
  --confidence 0.1 \
  --threads "$THREADS" \
  "$UNMAPPED_FASTQ"

# ============================
# GENERATE KRONA PLOT
# ============================
echo "Converting Kraken2 report to Krona format..."
python3 "$KREPORT_TO_KRONA_SCRIPT" -r "$KRAKEN_REPORT" -o "${SPECIES}_SAMPLE.krona"

echo "Generating Krona visualization..."
ktImportText "${SPECIES}_SAMPLE.krona" -o "$KRONA_HTML"

# Copy results to central results directory
cp "$KRONA_HTML" "$RESULTS_DIR/" 2>/dev/null
cp "$KRAKEN_OUTPUT" "$RESULTS_DIR/" 2>/dev/null
cp "$KRAKEN_REPORT" "$RESULTS_DIR/" 2>/dev/null

conda deactivate
set -u

echo "✅ Unmapped reads analysis complete. Results are available in $RESULTS_DIR"
