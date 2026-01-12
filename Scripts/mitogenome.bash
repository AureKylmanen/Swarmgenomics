#!/bin/bash
# Script Name: mitogenome.bash
# Description: Assemble mitochondrial genome, BLAST against shared database,
#              and infer phylogenetic placement
# Author: Aure Kylmänen
# ====================================

set -euo pipefail

# ============================
# USAGE
# ============================
# ./mitogenome.bash <species> <params.txt>

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
# REQUIRED PARAMETERS
# ============================
: "${WORKING_DIR:?Missing WORKING_DIR}"
: "${FASTQ_DIR:?Missing FASTQ_DIR}"
: "${MT_OUT_DIR:?Missing MT_OUT_DIR}"
: "${THREADS:?Missing THREADS}"

: "${MITO_DB:?Missing MITO_DB}"
: "${MT_FASTA:?Missing MT_FASTA}"

: "${GETORGANELLE_TYPE:?Missing GETORGANELLE_TYPE}"
: "${GETORGANELLE_READ_ROUNDS:?Missing GETORGANELLE_READ_ROUNDS}"
: "${N_BLAST_HITS:?Missing N_BLAST_HITS}"

: "${BLAST_MAX_TARGET_SEQS:?Missing BLAST_MAX_TARGET_SEQS}"
: "${BLAST_WORD_SIZE:?Missing BLAST_WORD_SIZE}"

: "${MAFFT_OPTS:?Missing MAFFT_OPTS}"
: "${FASTTREE_OPTS:?Missing FASTTREE_OPTS}"

# ============================
# DIRECTORIES
# ============================
mkdir -p "$MT_OUT_DIR" "$RESULTS_DIR"

# ============================
# LOCATE INPUT READS
# ============================
FASTQ1=$(ls "$FASTQ_DIR"/*_1_paired.fastq.gz | head -n 1)
FASTQ2=$(ls "$FASTQ_DIR"/*_2_paired.fastq.gz | head -n 1)

if [[ -z "$FASTQ1" || -z "$FASTQ2" ]]; then
    echo "❌ ERROR: Paired FASTQ files not found in $FASTQ_DIR"
    exit 1
fi

MITO_OUT="$MT_OUT_DIR/${SPECIES}_mitogenome"

echo "FASTQ1: $FASTQ1"
echo "FASTQ2: $FASTQ2"
echo "Output directory: $MITO_OUT"

# ============================
# RUN GETORGANELLE
# ============================
echo "Running GetOrganelle..."
source "$(conda info --base)/etc/profile.d/conda.sh"
set +u
conda activate getorganelle

get_organelle_from_reads.py \
    -1 "$FASTQ1" \
    -2 "$FASTQ2" \
    -F "$GETORGANELLE_TYPE" \
    -o "$MITO_OUT" \
    -R "$GETORGANELLE_READ_ROUNDS" \
    -t "$THREADS"

conda deactivate
set -u

# ============================
# IDENTIFY LONGEST SCAFFOLD
# ============================
FILE_MITO=$(find "$MITO_OUT" -name "*.path_sequence.fasta" | head -n 1)

if [[ -z "$FILE_MITO" ]]; then
    echo "❌ ERROR: No mitogenome FASTA produced by GetOrganelle"
    exit 1
fi

echo "✅ Selecting longest mitochondrial scaffold..."
FIXED_MITO="$MITO_OUT/mitogenome_cleaned.fasta"
awk '/^>/{print ">seq" NR; next} {print}' "$FILE_MITO" > "$FIXED_MITO"

LONGEST_MITO="$MITO_OUT/mito_longest.fasta"
seqkit sort -l -r "$FIXED_MITO" | seqkit head -n 1 > "$LONGEST_MITO"

QUERY_FASTA="$MITO_OUT/mito_query.fasta"
sed "s/^>.*/>${SPECIES}_query/" "$LONGEST_MITO" > "$QUERY_FASTA"

# ============================
# BLAST AGAINST SHARED DATABASE
# ============================
echo "Running BLAST against mitochondrial database..."
BLAST_OUT="$MITO_OUT/blast_hits.tsv"

blastn \
    -query "$QUERY_FASTA" \
    -db "$MITO_DB" \
    -outfmt "6 sseqid bitscore" \
    -max_target_seqs "$BLAST_MAX_TARGET_SEQS" \
    -word_size "$BLAST_WORD_SIZE" \
    -num_threads "$THREADS" \
    > "$BLAST_OUT"

TOP_IDS="$MITO_OUT/top_hits.ids"
sort -k2,2nr "$BLAST_OUT" | \
    awk '!seen[$1]++' | \
    head -n "$N_BLAST_HITS" | \
    cut -f1 > "$TOP_IDS"

TOP_FASTA="$MITO_OUT/top_hits.fasta"
seqkit grep -f "$TOP_IDS" "$MT_FASTA" > "$TOP_FASTA"
cat "$QUERY_FASTA" >> "$TOP_FASTA"

# ============================
# ALIGNMENT AND TREE
# ============================
echo "Aligning sequences..."
mafft $MAFFT_OPTS "$TOP_FASTA" > "$MITO_OUT/aligned.fasta"

echo "Building phylogenetic tree..."
FastTree $FASTTREE_OPTS "$MITO_OUT/aligned.fasta" > "$MITO_OUT/tree.nwk"

# ============================
# PLOTTING
# ============================

if [[ -n "${PHYLOTREE_SCRIPT:-}" ]]; then
    echo "Plotting phylogenetic tree..."
    
    # Make sure output PNG path is in RESULTS_DIR
    TREE_PNG="${RESULTS_DIR}/mito_tree.png"
    
    python "$PHYLOTREE_SCRIPT" \
        "$MT_FASTA" \
        "$QUERY_FASTA" \
        "$MITO_OUT/tree.nwk" \
        "$SPECIES" \
        "$TREE_PNG"
    
    echo "✅ Phylogenetic tree plotted: $TREE_PNG"
fi


# ============================
# COPY RESULTS
# ============================
echo "Copying results to $RESULTS_DIR..."

# Copy key output files
cp "$QUERY_FASTA" "$RESULTS_DIR/" 2>/dev/null || true
cp "$MITO_OUT/aligned.fasta" "$RESULTS_DIR/" 2>/dev/null || true
cp "$MITO_OUT/tree.nwk" "$RESULTS_DIR/" 2>/dev/null || true
cp "$BLAST_OUT" "$RESULTS_DIR/" 2>/dev/null || true
cp "$TOP_FASTA" "$RESULTS_DIR/" 2>/dev/null || true

# Copy tree plot if it exists (already saved to RESULTS_DIR by Python script)
if [[ -f "$TREE_PNG" ]]; then
    echo "✅ Tree plot saved: $TREE_PNG"
else
    echo "⚠️  Warning: Tree plot not found at $TREE_PNG"
fi

echo "✅ Mitogenome analysis completed successfully. Results saved in $RESULTS_DIR."
