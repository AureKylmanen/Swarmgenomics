#!/bin/bash

# Directory containing VCF files
VCF_DIR="/vol/storage/swarmGenomics/giant_panda/vcf"
OUTPUT_FILE="heterozygosity_results.txt"
BCFTOOLS="/vol/storage/software/bcftools-1.19"

# Remove the output file if it exists
rm -f $OUTPUT_FILE

# Write header to the output file
echo -e "Chromosome\tHeterozygosity" >> $OUTPUT_FILE

# Loop over each VCF file in the directory
for vcf_file in "$VCF_DIR"/*.vcf; do
    # Extract chromosome name (assuming filename is like chr1.vcf or scaffold1.vcf)
    chrom=$(basename "$vcf_file" .vcf)
    
    # Calculate heterozygosity using FORMAT/GT
    heterozygosity=$($BCFTOOLS/bcftools query -f '[%GT\n]' "$vcf_file" | \
                     awk 'BEGIN {het=0; hom=0} 
                          {if ($0 == "0/1" || $0 == "1/0") het+=1; else if ($0 == "0/0" || $0 == "1/1") hom+=1} 
                          END {if (het+hom > 0) print het/(het+hom); else print 0}')
    
    # Print the result to the output file
    echo -e "$chrom\t$heterozygosity" >> $OUTPUT_FILE
done

echo "Heterozygosity calculation complete. Results saved in $OUTPUT_FILE."


