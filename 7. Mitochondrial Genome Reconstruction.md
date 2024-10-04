# Mitochondrial Genome Reconstruction
The mitochondrial genome is a small, circular DNA found in mitochondria, distinct from the nuclear genome. It is maternally inherited and highly conserved, making it a valuable target for genetic studies.

Mitochrondial genome can be used for example in species identification, population genetics, forensic identification and biomedical purposes. Here we will use the constructed mitogenome in the identification of NUMTs in the following section.

## Installations
For the identification we will use GetOrganelle (https://github.com/Kinggerm/GetOrganelle).

```
# GetOrganelle installation
conda create -n getorganelle -c bioconda getorganelle
conda activate getorganelle
```

## Mitogenome assembly
The input files you will use are the fastq files you created in the previous steps. Change the number of threads as necessary.
```
# Create a directory for your results
mkdir mt

# Run get_organelle
get_organelle_config.py -a animal_mt 
$ get_organelle_from_reads.py -1 DESTINATION_PATH/*_1_paired.fastq.gz \
   -2 DESTINATION_PATH/*_2_paired.fastq.gz -F animal_mt \
   -o DESTINATION_PATH/mt/*_mitogenome -R 10 -t 26 -s

#Example
get_organelle_config.py -a animal_mt 
get_organelle_from_reads.py -1 /vol/storage/swarmGenomics/golden_eagle/fastq/ERR3316068_1_paired.fastq.gz  \
   -2 /vol/storage/swarmGenomics/golden_eagle/fastq/ERR3316068_2_paired.fastq.gz -F animal_mt \
   -o /vol/storage/swarmGenomics/golden_eagle/mt/ERR3316068_mitogenome -R 10 -t 26
```
The key output files include: \
  ``*.path_sequence.fasta``, each fasta file represents one type of genome structure \
  ``*.fastg``, the organelle related assembly graph to report for improvement and debug \
``*.selected_graph.gfa``, the organelle-only assembly graph \
``get_org.log.txt``, the log file 
