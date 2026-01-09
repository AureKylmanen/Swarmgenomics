#!/bin/bash
# Script Name: build_mt_database.bash
# Description: Download and build a shared mitochondrial BLAST database
# Author: Aure Kylmänen
# ====================================
set -euo pipefail

# ============================
# USAGE
# ============================
# ./build_mt_database.bash <species> <params.txt>

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <species> <params.txt>"
    exit 1
fi

SPECIES="$1"
PARAMS_FILE="$2"

if [[ ! -f "$PARAMS_FILE" ]]; then
    echo "ERROR: Parameter file not found: $PARAMS_FILE"
    exit 1
fi

source "$PARAMS_FILE"

# ============================
# REQUIRED PARAMETERS
# ============================

: "${WORKING_DIR:?Missing WORKING_DIR in params file}"

MITO_DB_DIR="${WORKING_DIR}/mt_database"
MT_FASTA="${MITO_DB_DIR}/genome.fna"
MT_FASTA_GZ="${MITO_DB_DIR}/genome.fna.gz"
MITO_DB="${MITO_DB_DIR}/mito_db"


# ============================
# SETUP
# ============================

mkdir -p "$MITO_DB_DIR"
cd "$MITO_DB_DIR"

# ============================
# DOWNLOAD DATABASE IF NEEDED
# ============================

if [[ ! -f "$MT_FASTA" ]]; then
    echo "Mitochondrial FASTA not found."

    if [[ ! -f "$MT_FASTA_GZ" ]]; then
        echo "Downloading mitochondrial genome database..."
        wget -O "$MT_FASTA_GZ" "$MT_DB_URL"
    else
        echo "Compressed FASTA already exists. Skipping download."
    fi

    echo "Unzipping genome.fna.gz..."
    gunzip -f "$MT_FASTA_GZ"
else
    echo "Mitochondrial FASTA already exists. Skipping download."
fi

# ============================
# BUILD BLAST DATABASE
# ============================

if [[ -f "${MITO_DB}.nhr" ]]; then
    echo "BLAST database already exists. Skipping makeblastdb."
else
    echo "Building mitochondrial BLAST database..."
    makeblastdb \
        -in "$MT_FASTA" \
        -dbtype nucl \
        -out "$MITO_DB"
fi

echo "======================================"
echo "Mitochondrial database setup complete"
echo "FASTA: $MT_FASTA"
echo "BLAST DB prefix: $MITO_DB"
echo "======================================"
