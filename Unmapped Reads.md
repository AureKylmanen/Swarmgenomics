# Unmapped Reads
Unmapped reads are sequences from high-throughput sequencing data that do not align to a reference genome during the alignment process. These reads can reveal information on microbial communities to structural variations and novel genetic sequences. 

Utilizing tools like Kraken2 which can classify these reads, provide valuable insights into the broader genetic landscape of the samples. Kraken2 will create an output of a taxonomic classification of each read, helping to identify microbial species or other sources of the unmapped reads.

### Installations
Kraken2 is easily installed using miniconda with following steps.

```
# Install Kraken2 (version 2.1.3) to its own environment called kraken
conda create -n kraken -c bioconda kraken2=2.1.3

# If you encounter an error with solving environment, try:
conda config --env --add channels conda-forge

# Activate kraken
conda activate kraken

#Install krakentools
conda install bioconda::krakentools
```

You will also need to create or download a database. You may choose one from https://benlangmead.github.io/aws-indexes/k2 which suits your needs. Here we will download the PlusPFP-16 database.
```
# Make a directory for the database
mkdir kraken_database

# Download the database in the directory
cd kraken_database
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16gb_20240112.tar.gz
tar -xvzf k2_pluspf_16gb_20240112.tar.gz
```
Depending on the size of the database, this may take hours.

Download kreport2krona.py into your virtual machine.

### Classifying reads
For classifying the unmapped reads from a whole-genome sequencing project we use the kraken2 command. You need to replace $DBPATH with the path to where you have saved the database and the name you gave the database. The input file is a fasta file of the unmapped reads. You can also choose the number of threads you want kraken2 to use. The flags --classified-out and --unclassified-out are optional.

##### Creating the input file
```
# Unmapped reads
samtools view -b -f 4 bwa.sorted.bam > unmapped.bam

# Change to fastq
bedtools bamtofastq -i unmapped.bam -fq unmapped.fastq
```

##### Running kraken2
Activate kraken environment and run the script.
```
conda activate kraken

kraken2 \
--db ${DBPATH} \
--output ${OUTPUTNAME} \
--use-names \
--report ${REPORTNAME} \
--classified-out ${CLASSIFIEDNAME} \
--unclassified-out ${UNCLASSIFIEDNAME} \
--confidence 0.1 \
--threads ${THREADNUM} \
${INPUTNAME}

# Example
kraken2 \
--db /vol/storage/kraken_database/ \
--output acunmapped.kraken \
--use-names \
--report acunmapped.kreport \
--classified-out classifiedac \
--unclassified-out unclassifiedac \
--confidence 0.1 \
--threads 26 \
unmapped.fastq
```

### Visualising the results
You may then run the kreport2krona.py script which will produce an html file for visualising the results.
```
# Copy kreport2krona.py into your directory
# Run the plot
python3 kreport2krona.py -r REPORTNAME -o SAMPLE.krona
ktImportText SAMPLE.krona -o SAMPLE.krona.html

# Example
python3 kreport2krona.py -r acunmapped.kreport -o SAMPLE.krona
ktImportText SAMPLE.krona -o SAMPLE.krona.html
```
# Resources
For information on how to create your own database: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown, and on KrakenTools: https://github.com/jenniferlu717/KrakenTools/blob/master/README.md.
