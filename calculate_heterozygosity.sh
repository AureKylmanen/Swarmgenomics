#!/bin/bash

# Directory containing VCF files
VCF_DIR="/vol/storage/swarmGenomics/giant_panda/vcf" # change this
OUTPUT_FILE="/vol/storage/swarmGenomics/giant_panda/heterozygosity_results.txt" # change this
BCFTOOLS="/vol/storage/software/bcftools-1.19" # change this

# Remove the output file if it exists
rm -f $OUTPUT_FILE

# Write header to the output file
echo -e "Chromosome\tHeterozygosity\tLength" >> $OUTPUT_FILE

# Loop over each VCF file in the directory
for vcf_file in "$VCF_DIR"/*.vcf; do
    # Extract chromosome name (assuming filename is like chr1.vcf or scaffold1.vcf)
    chrom=$(basename "$vcf_file" .vcf)
    
    # Extract scaffold length from VCF header
    length=$(grep -m 1 -oP "(?<=##contig=<ID=$chrom,length=)\d+" "$vcf_file")

    # Calculate heterozygosity using FORMAT/GT
    heterozygosity=$($BCFTOOLS/bcftools query -f '[%GT\n]' "$vcf_file" | \
                     awk 'BEGIN {het=0; hom=0} 
                          {if ($0 == "0/1" || $0 == "1/0") het+=1; else if ($0 == "0/0" || $0 == "1/1") hom+=1} 
                          END {if (het+hom > 0) print het/(het+hom); else print 0}')
    
    # Print the result to the output file
    echo -e "$chrom\t$heterozygosity\t$length" >> $OUTPUT_FILE
done

echo "Heterozygosity calculation complete. Results saved in $OUTPUT_FILE."
