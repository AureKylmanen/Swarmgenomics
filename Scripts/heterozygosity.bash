#!/bin/bash
# Script Name: heterozygosity.bash
# Description: Calculate and plot heterozygosity per scaffold
# Author: Aure Kylmänen
# ====================================
set -euo pipefail

# ============================
# USAGE
# ============================
# ./heterozygosity.bash <species> <params.txt>

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
: "${VCF_DIR:?Missing VCF_DIR in params file}"
: "${HET_RESULTS:?Missing HET_RESULTS in params file}"
: "${HET_PLOT_SCRIPT:?Missing HET_PLOT_SCRIPT in params file}"

# ============================
# DIRECTORIES
# ============================
RESULTS_DIR="${WORKING_DIR}/${SPECIES}/results"
mkdir -p "$RESULTS_DIR"

# ============================
# CALCULATE HETEROZYGOSITY
# ============================
echo "Calculating heterozygosity..."
rm -f "$HET_RESULTS"
echo -e "Chromosome\tHeterozygosity\tLength" >> "$HET_RESULTS"

shopt -s nullglob
VCF_FILES=("$VCF_DIR"/*.vcf)
shopt -u nullglob

if [[ ${#VCF_FILES[@]} -eq 0 ]]; then
    echo "❌ ERROR: No VCF files found in $VCF_DIR"
    exit 1
fi

for vcf_file in "${VCF_FILES[@]}"; do
    chrom=$(basename "$vcf_file" .vcf)
    scaffold_length=$(grep -m 1 "^##contig=<ID=${chrom}," "$vcf_file" | grep -oP "length=\K\d+")
    
    if [[ -z "$scaffold_length" ]]; then
        echo "Skipping $vcf_file: Scaffold length not found"
        continue
    fi

    heterozygosity=$(${BCFTOOLS}/bcftools query -f '[%GT\n]' "$vcf_file" | \
        awk 'BEGIN {het=0; hom=0}
             {if ($0 == "0/1" || $0 == "1/0") het+=1; else if ($0 == "0/0" || $0 == "1/1") hom+=1}
             END {if (het+hom > 0) printf "%.6f\n", het/(het+hom); else print 0}')

    echo -e "$chrom\t$heterozygosity\t$scaffold_length" >> "$HET_RESULTS"
done

# ============================
# PLOT RESULTS
# ============================
# Plot heterozygosity results
if [[ -f "$HET_RESULTS" ]]; then
    echo "Plotting heterozygosity results..."
    Rscript "$HET_PLOT_SCRIPT" "$HET_RESULTS" "$TOP_SCAFFOLDS" "${RESULTS_DIR}/heterozygosity_plot.png"
else
    echo "❌ ERROR: Heterozygosity results file not found!"
    exit 1
fi

echo "✅ Heterozygosity analysis completed successfully. Results saved in $RESULTS_DIR."
