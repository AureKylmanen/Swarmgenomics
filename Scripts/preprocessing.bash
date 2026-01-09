#!/bin/bash
# Script Name: preprocessing.bash
# Description: Read preprocessing, alignment, variant calling, VCF preparation and diploid.fq creation
# Author: Aure Kylmänen
# ====================================
set -euo pipefail

# ============================
# USAGE
# ============================
# ./preprocessing.bash <species> <reference.fna.gz> <SRA_ID> <params.txt>

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <species> <reference.fna.gz> <SRA_ID> <params.txt>"
    exit 1
fi

SPECIES="$1"
GENOME_FILE="$2"
SRA_FILE="$3"
PARAMS_FILE="$4"

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

: "${WORKING_DIR:?Missing WORKING_DIR}"
: "${THREADS:?Missing THREADS}"
: "${SRATOOLS:?Missing SRATOOLS}"
: "${BCFTOOLS:?Missing BCFTOOLS}"
: "${ADAPTERS:?Missing ADAPTERS}"
: "${TRIM_SLIDINGWINDOW:?Missing TRIM_SLIDINGWINDOW}"
: "${TRIM_MINLEN:?Missing TRIM_MINLEN}"
: "${MPILEUP_C:?Missing MPILEUP_C}"
: "${MIN_SCAFFOLD_LENGTH:?Missing MIN_SCAFFOLD_LENGTH}"
: "${VCF_MIN_DP:?Missing VCF_MIN_DP}"

# ============================
# DIRECTORIES
# ============================

SPECIES_DIR="${WORKING_DIR}/${SPECIES}"
FASTQ_DIR="${SPECIES_DIR}/fastq"
FASTQC_DIR="${SPECIES_DIR}/fastqc"
VCF_DIR="${SPECIES_DIR}/vcf"

mkdir -p "$FASTQ_DIR" "$FASTQC_DIR" "$VCF_DIR"

# ============================
# INPUT VALIDATION
# ============================

if [[ ! -f "$GENOME_FILE" ]]; then
    echo "ERROR: Reference genome not found: $GENOME_FILE"
    exit 1
fi

# ============================
# STEP 1: SRA → FASTQ
# ============================

echo "Downloading and converting SRA to FASTQ..."

SRA_PATH="${SPECIES_DIR}/${SRA_FILE}/${SRA_FILE}.sra"

"${SRATOOLS}/fastq-dump" \
    --split-files \
    --outdir "$FASTQ_DIR" \
    "$SRA_PATH"

# ============================
# STEP 2: FASTQC (RAW)
# ============================

echo "Running FastQC on raw reads..."

fastqc -t "$THREADS" -o "$FASTQC_DIR" \
    "${FASTQ_DIR}/${SRA_FILE}_1.fastq" \
    "${FASTQ_DIR}/${SRA_FILE}_2.fastq"

# ============================
# STEP 3: TRIMMOMATIC
# ============================

echo "Trimming reads with Trimmomatic..."

trimmomatic PE -threads "$THREADS" \
    "${FASTQ_DIR}/${SRA_FILE}_1.fastq" \
    "${FASTQ_DIR}/${SRA_FILE}_2.fastq" \
    "${FASTQ_DIR}/${SRA_FILE}_1_paired.fastq.gz" \
    "${FASTQ_DIR}/${SRA_FILE}_1_unpaired.fastq.gz" \
    "${FASTQ_DIR}/${SRA_FILE}_2_paired.fastq.gz" \
    "${FASTQ_DIR}/${SRA_FILE}_2_unpaired.fastq.gz" \
    SLIDINGWINDOW:${TRIM_SLIDINGWINDOW} \
    MINLEN:${TRIM_MINLEN} \
    ILLUMINACLIP:${ADAPTERS}:2:30:10:2:keepBothReads

# ============================
# STEP 4: FASTQC (TRIMMED)
# ============================

fastqc -t "$THREADS" -o "$FASTQC_DIR" \
    "${FASTQ_DIR}/${SRA_FILE}_1_paired.fastq.gz" \
    "${FASTQ_DIR}/${SRA_FILE}_2_paired.fastq.gz"

# ============================
# STEP 5: ALIGNMENT (BWA)
# ============================

echo "Indexing reference genome..."
bwa index "$GENOME_FILE"

echo "Aligning reads..."
bwa mem -t "$THREADS" "$GENOME_FILE" \
    "${FASTQ_DIR}/${SRA_FILE}_1_paired.fastq.gz" \
    "${FASTQ_DIR}/${SRA_FILE}_2_paired.fastq.gz" \
    > "${SPECIES_DIR}/${SRA_FILE}.sam"

# ============================
# STEP 6: BAM PROCESSING
# ============================

samtools sort -@ "$THREADS" \
    "${SPECIES_DIR}/${SRA_FILE}.sam" \
    -o "${SPECIES_DIR}/bwa.sorted.bam"

samtools index "${SPECIES_DIR}/bwa.sorted.bam"

# ============================
# STEP 7: PREPARE REFERENCE
# ============================

REF_UNZIPPED="${GENOME_FILE%.gz}"

if [[ ! -f "$REF_UNZIPPED" ]]; then
    gunzip -c "$GENOME_FILE" > "$REF_UNZIPPED"
fi

samtools faidx "$REF_UNZIPPED"

# ============================
# STEP 8: VARIANT CALLING
# ============================

echo "Calling variants..."

"${BCFTOOLS}/bcftools" mpileup \
    -C "$MPILEUP_C" \
    -f "$REF_UNZIPPED" \
    "${SPECIES_DIR}/bwa.sorted.bam" | \
"${BCFTOOLS}/bcftools" call \
    -c \
    -o "${SPECIES_DIR}/output.vcf"

# ============================
# STEP 9: REPEAT MASKING
# ============================

echo "Generating repeat BED file..."

perl -lne '
if(/^(>.*)/){ $head=$1 }
else { $fa{$head}.=$_ }
END{
for $s (keys %fa){
print "$s\n$fa{$s}"
}}' "$REF_UNZIPPED" | \
perl -lne '
if(/^>(\S+)/){ $n=$1 }
else{
while(/([a-z]+)/g){
printf("%s\t%d\t%d\n",$n,pos($_)-length($1),pos($_))
}}' > "${SPECIES_DIR}/repeats.bed"

bedtools subtract -header \
    -a "${SPECIES_DIR}/output.vcf" \
    -b "${SPECIES_DIR}/repeats.bed" \
    > "${SPECIES_DIR}/output_norepeats.vcf"

mv "${SPECIES_DIR}/output_norepeats.vcf" "${SPECIES_DIR}/output.vcf"

# ============================
# STEP 10: FILTER LARGE SCAFFOLDS
# ============================

TARGET_LEN="${SPECIES_DIR}/target_len.lst"

awk '$0 ~ ">" {
    if (NR > 1) print c;
    c=0;
    printf substr($0,2) "\t";
}
$0 !~ ">" { c+=length($0) }
END { print c }' "$REF_UNZIPPED" | \
awk "\$2 >= ${MIN_SCAFFOLD_LENGTH} {print \$1}" \
> "$TARGET_LEN"

# ============================
# STEP 11: SPLIT VCF BY SCAFFOLD
# ============================

"${BCFTOOLS}/bcftools" view \
    "${SPECIES_DIR}/output.vcf" -Oz \
    -o "${SPECIES_DIR}/output.vcf.gz"

tabix -p vcf "${SPECIES_DIR}/output.vcf.gz"

cd "$VCF_DIR"

while read chr; do
    tabix -h "${SPECIES_DIR}/output.vcf.gz" "$chr" > "${chr}.vcf"
done < "$TARGET_LEN"

# ============================
# STEP 12: VCF FILTERING
# ============================

for file in *.vcf; do
    vcftools \
        --vcf "$file" \
        --minDP "$VCF_MIN_DP" \
        --min-alleles "$VCF_MIN_ALLELES" \
        --max-alleles "$VCF_MAX_ALLELES" \
        --recode --recode-INFO-all \
        --out "${file%.vcf}_filtered"
done

for file in *_filtered.recode.vcf; do
    bgzip -@ "$THREADS" "$file"
    "${BCFTOOLS}/bcftools" index "${file}.gz"
done

# ============================
# STEP 13: GENERATE CONSENSUS FASTQ FILES PER SCAFFOLD
# ============================

cd "$VCF_DIR"
echo "Generating consensus FASTQ for each scaffold..."

while read chr; do
    echo "Processing $chr..."
    samtools faidx "$REF_UNZIPPED" "$chr" > "${chr}.fa"

    $BCFTOOLS/bcftools consensus \
        -f "${chr}.fa" \
        "${chr}_filtered.recode.vcf.gz" | \
    seqtk seq -F 'I' - > "${chr}.fq" || echo "⚠️ FAILED: $chr"
done < "$TARGET_LEN"

# ============================
# STEP 14: CONCATENATE AND COMPRESS DIPLOID FASTQ
# ============================

echo "Concatenating all scaffold FASTQs into ${DIPLOID_FASTQ}..."
cat *.fq | gzip -c > "$DIPLOID_FASTQ"

echo "✅ Diploid FASTQ creation completed successfully: $DIPLOID_FASTQ"

echo "Preprocessing completed successfully."
