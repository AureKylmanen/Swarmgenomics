# Version 1.0 
# 21.11.2024

# Input variables
# 1 = species name (e.g., golden_eagle)
# 2 = path to the reference genome file (e.g., /vol/storage/swarmgenomics/reference.fna.gz)
# 3 = SRA file identifier (e.g., ERR3316068)

# To run use nohup and & (e.g., nohup ./preprocessing.bash "golden_eagle" "/vol/storage/swarmgenomics/reference.fna.gz" "ERR3316068" &)

# Adjust number of threads according to your machine
# To check how many threads use command "nproc"
THREADS=26  # Change this to match available resources

SPECIES=$1
GENOME_FILE=$2
SRA_FILE=$3

# Directories
WORKING_DIR="/vol/storage/etc" # Change this
FASTQ_DIR="${WORKING_DIR}/${SPECIES}/fastq"
FASTQC_DIR="${WORKING_DIR}/${SPECIES}/fastqc"
VCF_DIR="${WORKING_DIR}/${SPECIES}/vcf"
TARGET_LEN="${WORKING_DIR}/${SPECIES}/target_len.lst"
ADAPTERS="/vol/storage/software/trimmomatics/adapters.fa" # Change this if installations in different location
BCFTOOLS="/vol/storage/software/bcftools-1.19" # Change this if installations in different location
SRATOOLS="/vol/storage/software/sratoolkit.3.1.1-centos_linux64/bin" # Change this if installations in different location
SRA_PATH="${WORKING_DIR}/${SPECIES}/${SRA_FILE}/${SRA_FILE}.sra"

# Step 1: Dump SRA files into FASTQ
$SRATOOLS/fastq-dump --outdir $FASTQ_DIR --split-files $SRA_PATH

# Step 2: Run FastQC on raw FASTQ files
mkdir -p $FASTQC_DIR
fastqc -t $THREADS -o $FASTQC_DIR -f fastq ${FASTQ_DIR}/${SRA_FILE}_1.fastq ${FASTQ_DIR}/${SRA_FILE}_2.fastq

# Step 3: Trimming the reads using Trimmomatic
trimmomatic PE -threads $THREADS ${FASTQ_DIR}/${SRA_FILE}_1.fastq ${FASTQ_DIR}/${SRA_FILE}_2.fastq \
    ${FASTQ_DIR}/${SRA_FILE}_1_paired.fastq.gz ${FASTQ_DIR}/${SRA_FILE}_1_unpaired.fastq.gz \
    ${FASTQ_DIR}/${SRA_FILE}_2_paired.fastq.gz ${FASTQ_DIR}/${SRA_FILE}_2_unpaired.fastq.gz \
    SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:${ADAPTERS}:2:30:10:2:keepBothReads

# Step 4: Run FastQC on the trimmed FASTQ files
fastqc -t $THREADS -o $FASTQC_DIR -f fastq \
    ${FASTQ_DIR}/${SRA_FILE}_1_paired.fastq.gz \
    ${FASTQ_DIR}/${SRA_FILE}_1_unpaired.fastq.gz \
    ${FASTQ_DIR}/${SRA_FILE}_2_paired.fastq.gz \
    ${FASTQ_DIR}/${SRA_FILE}_2_unpaired.fastq.gz

# Step 5: Index the reference genome using BWA
bwa index $GENOME_FILE

# Step 6: Align paired-end reads to the reference genome
bwa mem -t $THREADS $GENOME_FILE \
    ${FASTQ_DIR}/${SRA_FILE}_1_paired.fastq.gz \
    ${FASTQ_DIR}/${SRA_FILE}_2_paired.fastq.gz \
    > ${WORKING_DIR}/${SPECIES}/${SRA_FILE}.sam

# Step 7: Sort the aligned reads and save as BAM
samtools sort -@ $THREADS ${WORKING_DIR}/${SPECIES}/${SRA_FILE}.sam -o ${WORKING_DIR}/${SPECIES}/bwa.sorted.bam

# Step 8: Index the sorted BAM file
samtools index ${WORKING_DIR}/${SPECIES}/bwa.sorted.bam

# Step 9: Unzip and index the reference genome
gunzip -c $GENOME_FILE > ${GENOME_FILE%.gz}
samtools faidx ${GENOME_FILE%.gz}

# Step 10: Generate VCF using BCFtools
$BCFTOOLS/bcftools mpileup -C50 -f ${GENOME_FILE%.gz} ${WORKING_DIR}/${SPECIES}/bwa.sorted.bam | \
$BCFTOOLS/bcftools call -c -o ${WORKING_DIR}/${SPECIES}/output.vcf

# Step 11: Create a bed file with repeats
perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (sort(keys(%fa))){ print "$s\n$fa{$s}\n" }}' \
    ${GENOME_FILE%.gz} | \
perl -lne 'if(/^>(\S+)/){ $n=$1} else { while(/([a-z]+)/g){ printf("%s\t%d\t%d\n",$n,pos($_)-length($1),pos($_)) } }' \
    > ${WORKING_DIR}/${SPECIES}/repeats.bed

# Step 12: Remove repeats from VCF
bedtools subtract -header -a ${WORKING_DIR}/${SPECIES}/output.vcf -b ${WORKING_DIR}/${SPECIES}/repeats.bed > ${WORKING_DIR}/${SPECIES}/output_no_repeats.vcf
mv ${WORKING_DIR}/${SPECIES}/output_no_repeats.vcf ${WORKING_DIR}/${SPECIES}/output.vcf

# Step 13: Filter large chromosomes and prepare VCF files
awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' ${GENOME_FILE%.gz} | \
awk '$(NF) >= 5000000  {print $1}' > $TARGET_LEN

$BCFTOOLS/bcftools view ${WORKING_DIR}/${SPECIES}/output.vcf -Oz -o ${WORKING_DIR}/${SPECIES}/output.vcf.gz
tabix -p vcf ${WORKING_DIR}/${SPECIES}/output.vcf.gz

# Step 14: Prepare chromosome-specific VCF files
mkdir -p $VCF_DIR
cd $VCF_DIR
while read chr; do
    tabix -h ${WORKING_DIR}/${SPECIES}/output.vcf.gz $chr > ${chr}.vcf
done < $TARGET_LEN

# Step 15: Apply VCF filtering
for file in *.vcf; do
    vcftools --vcf $file --minDP 5 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out ${file%.vcf}_filtered
done

# Step 16: Compress and index filtered VCFs
for file in *_filtered.recode.vcf; do
    bgzip -@ $THREADS $file
    $BCFTOOLS/bcftools index ${file}.gz
done

# Step 17: Generate consensus fastq files for large chromosomes
while read chr; do
    samtools faidx ${GENOME_FILE%.gz} $chr | \
    $BCFTOOLS/bcftools consensus \
        -f ${GENOME_FILE%.gz} \
        ${chr}_filtered.recode.vcf.gz > ${chr}.fq
done < $TARGET_LEN

# Step 18: Concatenate and gzip final FASTQ
cat *.fq > diploid.fq
gzip diploid.fq
mv diploid.fq.gz ${WORKING_DIR}/${SPECIES}/
