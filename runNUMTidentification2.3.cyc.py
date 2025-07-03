import os,sys
import csv
import pandas as pd
import re

## Example
# dir_work = '/vol/storage/swarmGenomics/golden_eagle/'
# file_mito = '/vol/storage/swarmGenomics/golden_eagle/mt/ERR3316068_mitogenome/animal_mt.K115.scaffolds.graph1.1.path_sequence.fasta'
# file_genome = '/vol/storage/swarmGenomics/golden_eagle/acreference.fna'
# file_bedR = '/vol/storage/numt/get_bed2.R'
# species = 'Aquila_Chrysaetos'

dir_work = '/path/working/directory'
file_mito = '/path/to/*.path_sequence.fasta'
file_genome = '/path/to/reference.fna'
file_bedR = '/path/to/get_bed2.R'
species = 'Scientific_Name'

os.system('makeblastdb -in ' + file_mito + ' -parse_seqids -dbtype nucl')
os.system('makeblastdb -in ' + file_genome + ' -parse_seqids -dbtype nucl')

os.system('blastn -db '+ file_genome + ' -query '+ file_mito +' -outfmt 7 -word_size 20 -num_threads 1 -out ' + dir_work + species +'.blast.out')
os.system("grep -v '#' " + dir_work + species +".blast.out | cut -d '\t' -f2,3,4,9,10 > " + dir_work + species + ".blast.clean.out")

os.system('ls -1 ' + dir_work + species + '.blast.clean.out > ' + dir_work + 'files_list.txt')
os.system('cp ' + file_bedR + ' ' + dir_work)
os.system('Rscript get_bed2.R')
os.system('rm get_bed2.R')

os.system('bedtools getfasta -fi '+ file_genome +' -bed ' + dir_work + species + '.blast.clean.out_AllFrags.bed -fo ' + dir_work + species + '.sep.fa')
