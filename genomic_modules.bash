#!/bin/bash
# Version 1.0 
# 21.05.2025

# Input variables
SPECIES=$1

# Directories
WORKING_DIR="/vol/storage/swarmgenomics" # Change this
VCF_DIR="${WORKING_DIR}/${SPECIES}/vcf"
VCF_HET_DIR="${WORKING_DIR}/${SPECIES}/vcf_het"
ROH_DIR="${WORKING_DIR}/${SPECIES}/roh"
RESULTS_DIR="${WORKING_DIR}/${SPECIES}/results"

# Tools
BCFTOOLS="/vol/storage/software/bcftools-1.17/bcftools" # Change this
PSMC="/vol/storage/software/psmc/psmc" # Change this
PSMC_UTILS="/vol/storage/software/psmc/utils" # Change this

# File paths
BAM_FILE="${WORKING_DIR}/${SPECIES}/bwa.sorted.bam" # Change this if different bam file name
DIPLOID_FASTQ="${WORKING_DIR}/${SPECIES}/diploid.fq.gz"
DIPLOID_PSMCFA="${WORKING_DIR}/${SPECIES}/diploid.psmcfa"
DIPLOID_PSMC="${WORKING_DIR}/${SPECIES}/diploid.psmc"
IDXSTATS_INPUT="${WORKING_DIR}/${SPECIES}/${SPECIES}_idxstats_output.txt"
IDXSTATS_CSV="${WORKING_DIR}/${SPECIES}/idxstats_clean.csv"
HET_RESULTS="heterozygosity_results.txt"
RoH_PLOT="${WORKING_DIR}/${SPECIES}/all_roh_bar_plots.R"
idxstats_PLOT="${WORKING_DIR}/${SPECIES}/idxstats_plot.R"

# Create necessary directories
mkdir -p "${VCF_HET_DIR}" "${ROH_DIR}" "${RESULTS_DIR}"

# Step 0: Generate idxstats and auto-detect read length
samtools idxstats "${BAM_FILE}" > "${IDXSTATS_INPUT}"

READ_LENGTH=$(samtools view "${BAM_FILE}" | head -n 1000 | \
    awk '{if($10 != "*") print length($10)}' | \
    sort | uniq -c | sort -nr | head -n1 | awk '{print $2}')

echo "Detected read length: $READ_LENGTH"

awk 'BEGIN {OFS=","; print "chrom,length,mapped,unmapped"} {print $1,$2,$3,$4}' "$IDXSTATS_INPUT" > "$IDXSTATS_CSV"

cp "${IDXSTATS_CSV}" "${IDXSTATS_INPUT}" "${idxstats_PLOT}" "${RESULTS_DIR}/"

# Step 2: Generate RoH
cd "${VCF_DIR}"

${BCFTOOLS} merge --force-samples -O z -o merged.vcf.gz *.vcf.gz
${BCFTOOLS} roh -G30 --AF-dflt 0.4 merged.vcf.gz > roh_results.txt

grep -E "^# RG" roh_results.txt > RG.txt
grep -E "^RG" roh_results.txt >> RG.txt
awk '/^# RG/ {print; seen=1} seen && /^RG/' roh_results.txt > RG.txt

# Plot RoHs
Rscript "${RoH_PLOT}"
cp roh_results.txt RG.txt roh_bar_plots.png "${RESULTS_DIR}/" 2>/dev/null

# Step 3: Copy VCF files for heterozygosity calculations
cp *.vcf "${VCF_HET_DIR}"

# Step 4: Calculate heterozygosity
echo "Calculating heterozygosity..."
cd "${VCF_HET_DIR}"
rm -f ${HET_RESULTS}
echo -e "Chromosome\tHeterozygosity\tLength" >> ${HET_RESULTS}

for vcf_file in "${VCF_HET_DIR}"/*.vcf; do
    chrom=$(basename "$vcf_file" .vcf)
    scaffold_length=$(grep -m 1 -oP "(?<=##contig=<ID=$chrom,length=)\d+" "$vcf_file")

    if [ -z "$scaffold_length" ]; then
        echo "Skipping $vcf_file: Scaffold length not found"
        continue
    fi

    heterozygosity=$(${BCFTOOLS} query -f '[%GT\n]' "$vcf_file" | \
        awk 'BEGIN {het=0; hom=0}
             {if ($0 == "0/1" || $0 == "1/0") het+=1; else if ($0 == "0/0" || $0 == "1/1") hom+=1}
             END {if (het+hom > 0) print het/(het+hom); else print 0}')

    echo -e "$chrom\t$heterozygosity\t$scaffold_length" >> ${HET_RESULTS}
done

echo "Heterozygosity calculation complete."

# Step 5: Plot heterozygosity results
if [ -f "${HET_RESULTS}" ]; then
    echo "Plotting heterozygosity results..."
    Rscript "${WORKING_DIR}/${SPECIES}/plot_heterozygosity.R"
    cp "${HET_RESULTS}" heterozygosity_top_20_longest_plot.png "${RESULTS_DIR}/" 2>/dev/null
else
    echo "Heterozygosity results file not found!"
    exit 1
fi


# Step 1: Run idxstats R plotting
Rscript "${idxstats_PLOT}" "${RESULTS_DIR}/idxstats_clean.csv" "${RESULTS_DIR}/heterozygosity_results.txt" $READ_LENGTH


# Step 6: Prepare input for PSMC
cd "${WORKING_DIR}/${SPECIES}"
${PSMC_UTILS}/fq2psmcfa -q20 "${DIPLOID_FASTQ}" > "${DIPLOID_PSMCFA}"

# Step 7: Run PSMC
"${PSMC}" -N25 -t15 -r5 -p "2+2+25*2+4+6" -o "${DIPLOID_PSMC}" "${DIPLOID_PSMCFA}"

# Step 8: Generate PSMC historical scripts and plots
"${PSMC_UTILS}/psmc2history.pl" "${DIPLOID_PSMC}" | \
"${PSMC_UTILS}/history2ms.pl" > "${WORKING_DIR}/${SPECIES}_ms-cmd.sh"

"${PSMC_UTILS}/psmc_plot.pl" -p "${WORKING_DIR}/${SPECIES}/diploid" "${DIPLOID_PSMC}"

# Copy all final PSMC outputs
cp "${DIPLOID_PSMCFA}" "${DIPLOID_PSMC}" "${WORKING_DIR}/${SPECIES}_ms-cmd.sh" \
   diploid.* "${RESULTS_DIR}/" 2>/dev/null
  

echo "All analysis complete. Results are saved in ${RESULTS_DIR}/"
