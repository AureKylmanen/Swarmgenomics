#!/bin/bash
# Script Name: circos.bash
# Description: Prepare Circos input tracks and plot (with caching & replot-only mode)
# Author: Aure Kylmänen
# ====================================

set -euo pipefail

# ============================
# USAGE
# ============================
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <species> <reference.fna(.gz)> <params.txt>"
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
# REQUIRED PARAMETERS
# ============================
: "${WORKING_DIR:?Missing WORKING_DIR}"
: "${RESULTS_DIR:?Missing RESULTS_DIR}"
: "${BCFTOOLS:?Missing BCFTOOLS}"
: "${TOP_SCAFFOLDS:?Missing TOP_SCAFFOLDS}"
: "${BIN_SIZE:?Missing BIN_SIZE}"
: "${THREADS:?Missing THREADS}"
: "${CHR_NAMING_MODE:?Missing CHR_NAMING_MODE}"
: "${CIRCOS_BIN:?Missing CIRCOS_BIN}"
: "${CIRCOS_CONF_TEMPLATE:?Missing CIRCOS_CONF_TEMPLATE}"
: "${BAM_FILE:?Missing BAM_FILE}"
: "${VCF_DIR:?Missing VCF_DIR}"
: "${REPEATS_BED:?Missing REPEATS_BED}"

# ============================
# COLORS SHIM
# ============================
KARYO_COL="${__KARYO_COL__:-$KARYO_COL:-}"
GC_COL="${__GC_COL__:-$GC_COL:-}"
COV_COL="${__COV_COL__:-$COV_COL:-}"
HET_COL="${__HET_COL__:-$HET_COL:-}"
REPEAT_COL="${__REPEAT_COL__:-$REPEAT_COL:-}"

: "${KARYO_COL:?Missing KARYO_COL}"
: "${GC_COL:?Missing GC_COL}"
: "${COV_COL:?Missing COV_COL}"
: "${HET_COL:?Missing HET_COL}"
: "${REPEAT_COL:?Missing REPEAT_COL}"

# ============================
# SETUP WORKING DIR
# ============================
cd "${WORKING_DIR}/${SPECIES}"
mkdir -p "$RESULTS_DIR"

# ============================
# INPUT CACHE CHECK
# ============================
REQUIRED_INPUTS=(bins.bed coverage.txt gc_content.txt repeat_density.txt heterozygosity_density.txt)
inputs_exist_and_nonempty() {
    for f in "${REQUIRED_INPUTS[@]}"; do
        [[ -s "$f" ]] || return 1
    done
    return 0
}

# ============================
# SAFE MIN/MAX FUNCTION
# ============================
get_min_max() {
    local file="$1"
    if [[ ! -s "$file" ]]; then
        echo "0 1"
        return
    fi
    awk '{print $4}' "$file" | sort -n | awk 'NR==1{min=$1} {max=$1} END{if(min>max){tmp=min;min=max;max=tmp}; print min,max}'
}

# ============================
# PREPARE INPUT FILES IF MISSING
# ============================
if inputs_exist_and_nonempty; then
    echo "✅ Existing Circos inputs detected — skipping preprocessing (except karyotype)"
else
    echo "🛠 Preparing Circos input files"

    # ----- Genome -----
    [[ -f "$GENOME_FILE" ]] || { echo "❌ Genome file not found: $GENOME_FILE"; exit 1; }
    if [[ "$GENOME_FILE" == *.gz ]]; then
        REF_FA=$(basename "$GENOME_FILE" .gz)
        [[ -f "$REF_FA" ]] || bgzip -cd "$GENOME_FILE" > "$REF_FA"
    else
        REF_FA="$GENOME_FILE"
    fi
    [[ -f "${REF_FA}.fai" ]] || samtools faidx "$REF_FA"

    # ----- Top scaffolds -----
    sort -k2,2nr -k1,1 "${REF_FA}.fai" | head -n "$TOP_SCAFFOLDS" > top_scaffolds.fai

    # ----- Bins -----
    bedtools makewindows -g top_scaffolds.fai -w "$BIN_SIZE" | sort -k1,1 -k2,2n > bins.bed

    # ----- Coverage -----
    [[ -f "${BAM_FILE}.bai" ]] || samtools index "$BAM_FILE"
    mosdepth --by "$BIN_SIZE" -t "$THREADS" "$SPECIES" "$BAM_FILE"
    REGIONS="${SPECIES}.regions.bed.gz"
    if [[ ! -f "$REGIONS" ]]; then
        echo "❌ mosdepth output missing: $REGIONS"
        exit 1
    fi
    zcat "$REGIONS" > coverage.txt

    # ----- GC content -----
    bedtools nuc -fi "$REF_FA" -bed bins.bed | awk 'NR>1{print $1,$2,$3,$5}' > gc_content.txt

    # ----- Repeat density -----
    bedtools intersect -a bins.bed -b "$REPEATS_BED" -c | awk '{d=$4/($3-$2);print $1,$2,$3,d}' > repeat_density.txt

    # ----- Heterozygosity -----
    HET_VCF="output.vcf.gz"
    [[ -f "$HET_VCF" ]] || { echo "❌ Heterozygosity VCF not found: $HET_VCF"; exit 1; }
    bcftools view -g het -Ov "$HET_VCF" > heterozygous.vcf
    bedtools intersect -a bins.bed -b heterozygous.vcf -c | awk '{d=$4/($3-$2);print $1,$2,$3,d}' > heterozygosity_density.txt
fi

# ============================
# COMPUTE SAFE MIN/MAX VALUES
# ============================
read COV_MIN COV_MAX < <(get_min_max coverage.txt)
read GC_MIN GC_MAX < <(get_min_max gc_content.txt)
read HET_MIN HET_MAX < <(get_min_max heterozygosity_density.txt)
read REPEAT_MIN REPEAT_MAX < <(get_min_max repeat_density.txt)

# ----- ensure max > min -----
for var in COV GC HET REPEAT; do
    eval min=\$${var}_MIN
    eval max=\$${var}_MAX
    if (( $(echo "$min == $max" | bc -l) )); then
        max=$(echo "$min*1.01" | bc -l)
        eval ${var}_MAX=$max
    fi
done

# ============================
# ALWAYS GENERATE KARYOTYPE
# ============================
echo "Generating karyotype using mode: $CHR_NAMING_MODE"
case "$CHR_NAMING_MODE" in
    LG)
        awk -v color="$KARYO_COL" '{print "chr -",$1,"LG" NR,0,$2,color}' top_scaffolds.fai > karyotype.txt
        ;;
    SCAFFOLD)
        awk -v color="$KARYO_COL" '{print "chr -",$1,$1,0,$2,color}' top_scaffolds.fai > karyotype.txt
        ;;
    PREFIX)
        : "${CHR_PREFIX:?CHR_PREFIX must be set for PREFIX mode}"
        awk -v p="$CHR_PREFIX" -v color="$KARYO_COL" '{print "chr -",$1,p NR,0,$2,color}' top_scaffolds.fai > karyotype.txt
        ;;
    CUSTOM)
        : "${CHR_NAME_MAP:?CHR_NAME_MAP must be set for CUSTOM mode}"
        awk -v color="$KARYO_COL" 'NR==FNR{map[$1]=$2;next}{print "chr -",$1,map[$1],0,$2,color}' "$CHR_NAME_MAP" top_scaffolds.fai > karyotype.txt
        ;;
    *)
        echo "❌ ERROR: Unknown CHR_NAMING_MODE: $CHR_NAMING_MODE"
        exit 1
        ;;
esac

# ============================
# PREPARE circos.conf
# ============================
CIRCOS_CONF="${WORKING_DIR}/${SPECIES}/circos.conf"
cp "$CIRCOS_CONF_TEMPLATE" "$CIRCOS_CONF"

sed -i \
    -e "s|__COV_MIN__|$COV_MIN|g" \
    -e "s|__COV_MAX__|$COV_MAX|g" \
    -e "s|__GC_MIN__|$GC_MIN|g" \
    -e "s|__GC_MAX__|$GC_MAX|g" \
    -e "s|__HET_MIN__|$HET_MIN|g" \
    -e "s|__HET_MAX__|$HET_MAX|g" \
    -e "s|__REPEAT_MIN__|$REPEAT_MIN|g" \
    -e "s|__REPEAT_MAX__|$REPEAT_MAX|g" \
    -e "s|__GC_COL__|$GC_COL|g" \
    -e "s|__COV_COL__|$COV_COL|g" \
    -e "s|__HET_COL__|$HET_COL|g" \
    -e "s|__REPEAT_COL__|$REPEAT_COL|g" \
    -e "s|__COLORS_CONF__|$CIRCOS_COLORS_CONF|g" \
    -e "s|__BREWER_CONF__|$CIRCOS_BREWER_CONF|g" \
    -e "s|__FONTS_CONF__|$CIRCOS_FONTS_CONF|g" \
    -e "s|__IMAGE_CONF__|$CIRCOS_IMAGE_CONF|g" \
    -e "s|__HOUSEKEEPING_CONF__|$CIRCOS_HOUSEKEEPING_CONF|g" \
    "$CIRCOS_CONF"

# ============================
# RUN CIRCOS
# ============================
echo "🎨 Running Circos..."
"$CIRCOS_BIN" -conf "$CIRCOS_CONF"

# ============================
# LEGEND OVERLAY
# ============================
if [[ -n "${CIRCOS_OVERLAY_SCRIPT:-}" && -f "${CIRCOS_OVERLAY_SCRIPT}" ]]; then
    echo "Adding legend overlay..."
    if [[ -f circos.png ]]; then
        cp "${CIRCOS_LEGEND_IMAGE}" legend.png
        python "${CIRCOS_OVERLAY_SCRIPT}"
        [[ -f circos_with_legend.png ]] && echo "✅ Legend overlay added"
    else
        echo "⚠️ circos.png not found — skipping overlay"
    fi
fi

# ============================
# COPY RESULTS
# ============================
echo "Copying Circos input files and outputs..."
for f in karyotype.txt bins.bed coverage.txt gc_content.txt repeat_density.txt heterozygosity_density.txt circos.conf; do
    [[ -f "$f" ]] && cp -v "$f" "${RESULTS_DIR}/"
done

# Copy Circos images
[[ -f circos.png ]] && mv circos.png "${RESULTS_DIR}/${SPECIES}_circos.png"
[[ -f circos_with_legend.png ]] && mv circos_with_legend.png "${RESULTS_DIR}/${SPECIES}_circos_with_legend.png"
[[ -f circos.svg ]] && mv circos.svg "${RESULTS_DIR}/${SPECIES}_circos.svg"

echo "✅ Circos preparation complete."
echo "Results saved to: ${RESULTS_DIR}"
