#!/bin/bash
# Script Name: repeat_annotation.bash
# Description: De novo repeat annotation using RepeatModeler2 and RepeatMasker
# Author: Justin Wilcox
# ====================================

set -euo pipefail

# ============================
# USAGE
# ============================
# ./repeat_annotation.bash <species> <reference.fasta> <params.txt>
# Example: ./repeat_annotation.bash "golden_eagle" params.txt genome.fna.gz

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <species> <reference.fasta> <params.txt>"
    exit 1
fi

SPECIES="$1"
REFERENCE="$2"
PARAMS_FILE="$3"

# ============================
# LOAD PARAMETERS
# ============================
if [[ ! -f "$PARAMS_FILE" ]]; then
    echo "ERROR: Parameter file not found: $PARAMS_FILE"
    exit 1
fi

source "$PARAMS_FILE"

# ============================
# SETUP DIRECTORIES
# ============================

# Workdir for this analysis
WORKDIR="${WORKING_DIR}/${SPECIES}"
RM2_WORKDIR="${RM2_DIR}"
RM_WORKDIR="${RM_DIR}"

# Make sure directories exist
mkdir -p "$WORKDIR" "$RM2_WORKDIR" "$RM_WORKDIR"


# ============================
# PREPARE REFERENCE
# ============================
REF_BASENAME=$(basename "$REFERENCE" | sed 's/.gz$//')
REF_UNCOMPRESSED="${RM2_WORKDIR}/${REF_BASENAME}"

# Only unzip if the uncompressed file doesn't already exist
if [[ ! -f "$REF_UNCOMPRESSED" ]]; then
    if [[ "$REFERENCE" == *.gz ]]; then
        echo "Unzipping reference genome..."
        gunzip -c "$REFERENCE" > "$REF_UNCOMPRESSED"
    else
        echo "Copying reference genome..."
        cp "$REFERENCE" "$REF_UNCOMPRESSED"
    fi
else
    echo "Reference genome already prepared: $REF_UNCOMPRESSED"
fi
# ============================
# BUILD REPEATMODEL DATABASE
# ============================
DB_NAME="${SPECIES}_RM"

echo "Building RepeatModeler2 database..."
cd "$RM2_WORKDIR"
$REPEATMODELER2/BuildDatabase -engine ncbi -name "$DB_NAME" "$REF_UNCOMPRESSED"

# ============================
# RUN REPEATMODELER2
# ============================
echo "Running RepeatModeler2..."
$REPEATMODELER2/RepeatModeler \
    -threads "$RM_THREADS" \
    -database "$DB_NAME" \
    -engine ncbi \
    -genomeSampleSizeMax "$RM_GENOME_SAMPLE" \
    -LTRStruct

# Copy families.fa for RepeatMasker
FAMILIES_FA="${RM2_WORKDIR}/${DB_NAME}-families.fa"
cp "$FAMILIES_FA" "$RM_WORKDIR/"

# ============================
# RUN REPEATMASKER
# ============================
echo "Running RepeatMasker..."
cd "$RM_WORKDIR"
cp "$REF_UNCOMPRESSED" "$RM_WORKDIR/"

$REPEATMASKER/RepeatMasker \
    -e rmblast \
    -gccalc \
    -s \
    -a \
    -pa "$THREADS" \
    -lib "$RM_WORKDIR/${DB_NAME}-families.fa" \
    "$RM_WORKDIR/${REF_BASENAME}"

# ============================
# CREATE REPEAT LANDSCAPE
# ============================
GENOME_SIZE=$(grep -v ">" "$RM_WORKDIR/$REF_BASENAME" | wc -m)

perl "$REPEATMASKER/util/calcDivergenceFromAlign.pl" \
    -s "$RM_WORKDIR/${REF_BASENAME}.divsum" \
    -a "$RM_WORKDIR/${REF_BASENAME}.align" \
    "$RM_WORKDIR/${REF_BASENAME}.align"

perl "$REPEATMASKER/util/createRepeatLandscape.pl" \
    -div "$RM_WORKDIR/${REF_BASENAME}.divsum" \
    -t "$SPECIES Repeat Landscape" \
    -g "$GENOME_SIZE" \
    > "$RM_WORKDIR/${SPECIES}_repeat_landscape.html"

# ============================
# COPY RESULTS TO MAIN RESULTS DIR
# ============================
mkdir -p "${RESULTS_DIR}/repeat_annotation"
cp -r "$RM_WORKDIR"/* "${RESULTS_DIR}/repeat_annotation/"

echo "Repeat annotation completed successfully."
