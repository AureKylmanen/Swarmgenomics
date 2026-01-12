#!/bin/bash
# Script Name: numt_detection.bash
# Description: Detect and visualize NUMTs (Nuclear Mitochondrial DNA segments)
# Author: Yu-Chi Chen
# Adapted to bash by: Aure Kylmänen
# ====================================

set -euo pipefail

# ============================
# USAGE
# ============================
# ./numt_detection.bash <species> <params.txt>

# Usage:
# ./numt_detection.bash <species> <genome_file> <params.txt>
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <species> <genome_file> <params.txt>"
    exit 1
fi

SPECIES="$1"
GENOME_FILE="$2"
PARAMS_FILE="$3"


# ============================
# LOAD PARAMETERS
# ============================
if [[ ! -f "$PARAMS_FILE" ]]; then
    echo "❌ ERROR: Parameter file not found: $PARAMS_FILE"
    exit 1
fi

source "$PARAMS_FILE"

# ============================
# INITIALIZE CONDA
# ============================
source "$(conda info --base)/etc/profile.d/conda.sh"

# ============================
# REQUIRED PARAMETERS
# ============================
: "${WORKING_DIR:?Missing WORKING_DIR}"
: "${THREADS:?Missing THREADS}"
: "${RESULTS_DIR:?Missing RESULTS_DIR}"
: "${MT_OUT_DIR:?Missing MT_OUT_DIR}"
: "${GENOME_FILE:?Missing GENOME_FILE}"

: "${NUMT_BLAST_OUT:?Missing NUMT_BLAST_OUT}"
: "${NUMT_BLAST_CLEAN:?Missing NUMT_BLAST_CLEAN}"
: "${NUMT_SEP_FA:?Missing NUMT_SEP_FA}"

: "${NUMT_WORD_SIZE:?Missing NUMT_WORD_SIZE}"
: "${NUMT_EVALUE:?Missing NUMT_EVALUE}"

: "${NUMT_PLOT_SCRIPT:?Missing NUMT_PLOT_SCRIPT}"

# ============================
# PREPARE GENOME REFERENCE
# ============================
REF_BASENAME=$(basename "$GENOME_FILE" | sed 's/.gz$//')
GENOME_UNCOMPRESSED="${WORKING_DIR}/${SPECIES}/${REF_BASENAME}"

# Only unzip if the uncompressed file doesn't already exist
if [[ ! -f "$GENOME_UNCOMPRESSED" ]]; then
    if [[ "$GENOME_FILE" == *.gz ]]; then
        echo "Unzipping reference genome..."
        gunzip -c "$GENOME_FILE" > "$GENOME_UNCOMPRESSED"
    else
        echo "Copying reference genome..."
        cp "$GENOME_FILE" "$GENOME_UNCOMPRESSED"
    fi
else
    echo "Reference genome already prepared: $GENOME_UNCOMPRESSED"
fi

# ============================
# LOCATE MITOGENOME
# ============================
MITO_OUT="$MT_OUT_DIR/${SPECIES}_mitogenome"
FIXED_MITO="$MITO_OUT/mitogenome_cleaned.fasta"

if [[ ! -f "$FIXED_MITO" ]]; then
    echo "❌ ERROR: Mitogenome not found at $FIXED_MITO"
    echo "Please run mitogenome assembly first."
    exit 1
fi

# ============================
# CREATE BLAST DATABASES
# ============================
echo "Creating BLAST databases..."
makeblastdb -in "$FIXED_MITO" -parse_seqids -dbtype nucl
makeblastdb -in "$GENOME_UNCOMPRESSED" -parse_seqids -dbtype nucl

# ============================
# BLAST SEARCH FOR NUMTs
# ============================
echo "Running BLAST to detect NUMTs..."
blastn \
    -db "$GENOME_UNCOMPRESSED" \
    -query "$FIXED_MITO" \
    -outfmt 7 \
    -word_size "$NUMT_WORD_SIZE" \
    -evalue "$NUMT_EVALUE" \
    -num_threads "$THREADS" \
    -out "$NUMT_BLAST_OUT"

# ============================
# CLEAN BLAST OUTPUT
# ============================
echo "Processing BLAST results..."
grep -v "#" "$NUMT_BLAST_OUT" | cut -f2,3,4,9,10 > "$NUMT_BLAST_CLEAN"

# ============================
# GENERATE BED FILE
# ============================
# Activate R environment
set +u
conda activate r_env
set -u

cd "${WORKING_DIR}/${SPECIES}"   # <- run R from species folder
ls -1 "$NUMT_BLAST_CLEAN" > "files_list.txt"
Rscript "${SCRIPTS}/get_bed2.R"

set +u
conda deactivate
set -u


# ============================
# EXTRACT NUMT SEQUENCES
# ============================
BED_FILE="${NUMT_BLAST_CLEAN}_AllFrags.bed"

echo "Extracting NUMT sequences..."
bedtools getfasta \
    -fi "$GENOME_UNCOMPRESSED" \
    -bed "$BED_FILE" \
    -fo "$NUMT_SEP_FA"

# ============================
# COUNT NUMTs
# ============================
NUMT_COUNT=$(grep -c ">" "$NUMT_SEP_FA" || echo "0")

echo "Number of NUMTs detected: $NUMT_COUNT"
echo "$NUMT_COUNT" > "${RESULTS_DIR}/numt_count.txt"

# ============================
# VISUALIZE NUMTs
# ============================
if [[ "$NUMT_COUNT" -gt 0 ]]; then
    echo "Generating NUMT visualization..."
    
    set +u
    conda activate r_env
    set -u
    
    # Run R from species folder (where .sep.fa and .blast.out exist)
    cd "${WORKING_DIR}/${SPECIES}"  
    Rscript "$NUMT_PLOT_SCRIPT"
    
    set +u
    conda deactivate
    set -u
    
    # Copy visualization to results
    if [[ -f "${RESULTS_DIR}/numt.png" ]]; then
        echo "✅ NUMT visualization saved: ${RESULTS_DIR}/numt.png"
    fi
else
    echo "⚠️  No NUMTs detected - skipping visualization"
fi

# ============================
# COPY RESULTS
# ============================
echo "Copying results to $RESULTS_DIR..."

cp "$NUMT_BLAST_CLEAN" "$RESULTS_DIR/" 2>/dev/null || true
cp "$NUMT_SEP_FA" "$RESULTS_DIR/" 2>/dev/null || true
cp "$BED_FILE" "$RESULTS_DIR/" 2>/dev/null || true

# Copy visualization to RESULTS_DIR
if [[ -f "${WORKING_DIR}/${SPECIES}/results/numt.png" ]]; then
    echo "✅ NUMT visualization saved: ${WORKING_DIR}/${SPECIES}/results/numt.png"
fi
