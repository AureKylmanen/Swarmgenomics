# Repeat Analysis
Repetitive regions were long believed to be just "junk DNA" but recent genomic projects, such as the Human Genome Project, revealed that these repeat elements play crucial roles in genome organization, function, and evolution.
Because of this, repeat analysis is an important step in genome annotation and involves identifying and characterizing repetitive elements within a genome. Two widely used tools for this purpose are RepeatModeler and RepeatMasker.

## Installations
RepeatModeler and RepeatMasker have many dependencies, so make sure you have everything below installed and configured.

### Previous installations
You should already have the following packages installed from step 0. \
SRA toolkit\
Samtools\
Bcftools\
Python3\
Bedtools\
Trimmomatic\
Fastqc

### Package Manager Installations
```
sudo apt-get install bzip2
sudo apt install gcc
sudo apt install make
sudo apt-get install libncurses5-dev libncursesw5-dev -y
sudo apt-get install libbz2-dev -y
sudo apt-get install libbz2-dev -y
sudo apt-get install liblzma-dev -y
sudo apt-get install g++
sudo apt install unzip
sudo apt install python3-pip
sudo apt install freebayes
sudo apt install bowtie2
sudo apt install default-jre 
pip install multiqc
pip install h5py
pip install biser
```
### Source Installations
Install these in /vol/storage/software
```
# Install BS-Snper
cd /vol/storage/software
wget https://github.com/hellbelly/BS-Snper/archive/refs/heads/master.zip
unzip master.zip
cd BS-Snper-master
sh BS-Snper.sh

# Install HMMER
cd /vol/storage/software
wget http://eddylab.org/software/hmmer/hmmer-3.3.2.tar.gz
tar -xzf hmmer-3.3.2.tar.gz
rm /vol/storage/software/hmmer-3.3.2.tar.gz
cd hmmer-3.3.2
./configure --prefix=/vol/storage/software/hmmer-3.3.2
make
make install

# Install RM BLAST
cd /vol/storage/software
wget http://www.repeatmasker.org/rmblast/rmblast-2.14.0+-x64-linux.tar.gz
tar -xzf rmblast-2.14.0+-x64-linux.tar.gz
rm rmblast-2.14.0+-x64-linux.tar.gz

# Install TRF
cd /vol/storage/software
wget https://github.com/Benson-Genomics-Lab/TRF/archive/refs/tags/v4.09.1.tar.gz
tar -xzf /vol/storage/software/v4.09.1.tar.gz
rm /vol/storage/software/v4.09.1.tar.gz
mkdir /vol/storage/software/TRF-4.09.1/build
cd /vol/storage/software/TRF-4.09.1/build 
/vol/storage/software/TRF-4.09.1/configure --prefix=/vol/storage/software/TRF-4.09.1
make
make install

# Install RepeatMasker
cd /vol/storage/software
wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.5.tar.gz
tar -xzf RepeatMasker-4.1.5.tar.gz
rm RepeatMasker-4.1.5.tar.gz

# Configure RepeatMasker
perl /vol/storage/software/RepeatMasker/configure

#Paths for RepeatMasker
/vol/storage/software/TRF-4.09.1/bin/trf
/vol/storage/software/rmblast-2.14.0/bin
/vol/storage/software/hmmer-3.3.2/bin

#Selections
Add a Search Engine:
   1. Crossmatch: [ Un-configured ]
   2. RMBlast: [ Configured ]
   3. HMMER3.1 & DFAM: [ Configured, Default ]
   4. ABBlast: [ Un-configured ]

   5. Done

# Install Bismark
cd /vol/storage/software
wget https://github.com/FelixKrueger/Bismark/archive/0.22.3.tar.gz
tar -xzf 0.22.3.tar.gz 
rm 0.22.3.tar.gz

# Install Autotrim
cd /vol/storage/software
git clone https://github.com/schellt/autotrim.git

# Download SeqTK
cd /vol/storage/software
git clone https://github.com/lh3/seqtk.git;
cd seqtk
make

# Install RECON
cd /vol/storage/software
wget http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
tar RECON-1.08.tar.gz
rm /vol/storage/software/RECON-1.08.tar.gz
cd /vol/storage/software/RECON-1.08/src
make
make install

# Install Ninja
cd /vol/storage/software
wget https://github.com/TravisWheelerLab/NINJA/archive/0.95-cluster_only.tar.gz
tar -xzf /vol/storage/software/0.95-cluster_only.tar.gz
rm /vol/storage/software/0.95-cluster_only.tar.gz
mv /vol/storage/software/NINJA-0.95-cluster_only/NINJA/Ninja_new /vol/storage/software/NINJA-0.95-cluster_only/NINJA/Ninja

# Install LTR Retriever
cd /vol/storage/software
wget https://github.com/oushujun/LTR_retriever/archive/v2.8.tar.gz
tar -xzf v2.8.tar.gz
rm v2.8.tar.gz

# Install Mafft
cd /vol/storage/software
wget https://mafft.cbrc.jp/alignment/software/mafft-7.505-with-extensions-src.tgz
tar -xzf /vol/storage/software/mafft-7.505-with-extensions-src.tgz
rm /vol/storage/software/mafft-7.505-with-extensions-src.tgz
cd /vol/storage/software/mafft-7.505-with-extensions/core
sed -i 's#PREFIX = /usr/local#PREFIX = /vol/storage/software/mafft-7.505-with-extensions#' /vol/storage/software/mafft-7.505-with-extensions/core/Makefile
sed -i 's#BINDIR = $(PREFIX)/bin#BINDIR = /vol/storage/software/mafft-7.505-with-extensions/bin#' /vol/storage/software/mafft-7.505-with-extensions/core/Makefile
make clean
make
make install

cd /vol/storage/software/mafft-7.505-with-extensions/extensions
sed -i 's#PREFIX = /usr/local#PREFIX = /vol/storage/software/mafft-7.505-with-extensions#' /vol/storage/software/mafft-7.505-with-extensions/extensions/Makefile
sed -i 's#BINDIR = $(PREFIX)/bin#BINDIR = /vol/storage/software/mafft-7.505-with-extensions/bin#' /vol/storage/software/mafft-7.505-with-extensions/extensions/Makefile
make clean
make 
make install

# Install CD-Hit
cd /vol/storage/software
wget https://github.com/weizhongli/cdhit/archive/refs/tags/V4.8.1.tar.gz
tar -xzf /vol/storage/software/V4.8.1.tar.gz
rm /vol/storage/software/V4.8.1.tar.gz
cd /vol/storage/software/cdhit-4.8.1
make 

#Install Genometools
cd /vol/storage/software
wget http://genometools.org/pub/genometools-1.6.2.tar.gz
tar -xzf /vol/storage/software/genometools-1.6.2.tar.gz
rm /vol/storage/software/genometools-1.6.2.tar.gz
cd /vol/storage/software/genometools-1.6.2
make prefix=/vol/storage/software/genometools-1.6.2 cairo=no

# Install RepeatModeler2
cd /vol/storage/software
wget http://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.4.tar.gz
tar -xvf /vol/storage/software/RepeatModeler-2.0.4.tar.gz
rm /vol/storage/software/RepeatModeler-2.0.4.tar.gz
cd /vol/storage/software/RepeatModeler-2.0.4

# Install RepeatScout
cd /vol/storage/software
wget http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
tar -xzf /vol/storage/software/RepeatScout-1.0.6.tar.gz
rm /vol/storage/software/RepeatScout-1.0.6.tar.gz
cd /vol/storage/software/RepeatScout-1.0.6
make

# Install USCS
rm -r /vol/storage/software/USCS
mkdir /vol/storage/software/USCS
cd /vol/storage/software/USCS
rsync -aP hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/ ./

cpan install JSON
#yes
#sudo
cpan install File::Which
cpan install URI
cpan install Devel::Size
cpan install LWP::UserAgent
cd /vol/storage/software/RepeatModeler-2.0.4
perl ./configure

# Configure RepeatModeler2 with Paths w/ yes ("y") to using LTR Retriever
*/vol/storage/software/RepeatMasker
*/vol/storage/software/RECON-1.08/bin
*/vol/storage/software/RepeatScout-1.0.6
*/vol/storage/software/TRF-4.09.1/bin
*/vol/storage/software/cdhit-4.8.1
*/vol/storage/software/USCS
*/vol/storage/software/rmblast-2.14.0/bin
#yes
*/vol/storage/software/genometools-1.6.2/bin
*/vol/storage/software/LTR_retriever-2.8
*/vol/storage/software/mafft-7.505-with-extensions/bin
*/vol/storage/software/NINJA-0.95-cluster_only/NINJA

# Install MELT
# Download from https://melt.igs.umaryland.edu/downloads.php
# Move downloaded file into denbi using FileZilla or scp
cd /vol/storage/software/
tar -xzf /vol/storage/software/MELTv2.2.2.tar.gz
```
## Running repeat analysis
You may download and edit the script RepeatAnnotator.bash to run it, or run it step by step as in the example below.
```
# Zip reference genome
gzip -c /vol/storage/swarmGenomics/golden_eagle/acreference.fna > /vol/storage/swarmGenomics/golden_eagle/acreference.fna.gz

# Create a directory for RepeatModeler2
mkdir -p /vol/storage/swarmGenomics/golden_eagle/RepeatModeler2

# Make database
cd /vol/storage/swarmGenomics/golden_eagle/RepeatModeler2
rnREF=`echo /vol/storage/swarmGenomics/golden_eagle/acreference.fna.gz | perl -pe 's#.*/##' | perl -pe 's#.gz$##'`
dbNAME=`echo /vol/storage/swarmGenomics/golden_eagle/acreference.fna.gz | perl -pe 's#.*/##' | perl -pe 's#\.f.*.gz$##'`
zcat /vol/storage/swarmGenomics/golden_eagle/acreference.fna.gz > /vol/storage/swarmGenomics/golden_eagle/RepeatModeler2/$rnREF

#Build Database
/vol/storage/software/RepeatModeler-2.0.4/BuildDatabase -engine ncbi -name $dbNAME /vol/storage/swarmGenomics/golden_eagle/RepeatModeler2/$rnREF

# Run RepeatModeler
# use "nohup ... &", this takes long
/vol/storage/software/RepeatModeler-2.0.4/RepeatModeler -threads 26 -database $dbNAME -engine ncbi -genomeSampleSizeMax 243000000 -LTRStruct

# Make directory for results
mkdir -p /vol/storage/swarmGenomics/golden_eagle/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7

#Copy the main output files
cd /vol/storage/swarmGenomics/golden_eagle/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7
cp /vol/storage/swarmGenomics/golden_eagle/RepeatModeler2/$dbNAME-families.fa /vol/storage/swarmGenomics/golden_eagle/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$dbNAME-families.fa
cp /vol/storage/swarmGenomics/golden_eagle/RepeatModeler2/$rnREF /vol/storage/swarmGenomics/golden_eagle/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$rnREF

#Run RepeatMasker
/vol/storage/software/RepeatMasker/RepeatMasker -e rmblast -gccalc -s -a -pa 26 -lib  /vol/storage/swarmGenomics/golden_eagle/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$dbNAME-families.fa /vol/storage/swarmGenomics/golden_eagle/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$rnREF

#Create graphs
GENOMESIZE=`grep -v ">" /vol/storage/swarmGenomics/golden_eagle/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$rnREF | perl -pe -chomp | wc -m`

perl /vol/storage/software/RepeatMasker/util/calcDivergenceFromAlign.pl -s /vol/storage/swarmGenomics/golden_eagle/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$rnREF.divsum -a  /vol/storage/swarmGenomics/golden_eagle/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$rnREF.GC-Adjusted.align /vol/storage/swarmGenomics/golden_eagle/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$rnREF.align
perl /vol/storage/software/RepeatMasker/util/createRepeatLandscape.pl -div /vol/storage/swarmGenomics/golden_eagle/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$rnREF.divsum -t "$dbNAME Repeat Landscape" -g $GENOMESIZE > /vol/storage/swarmGenomics/golden_eagle/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$dbNAME.repeat_landscape.html
echo 'RepeatMasker and Graphs Complete'
```
## Results
The final results are in an html file (e.g. acreference.repeat_landscape.html) which you will need to copy to your computer for viewing. The end result is a repeat landscape, with Kimura substitution level on the x-axis and percentage of repeat families on the y-axis. Lower Kimura values indicate fewer mutations and thus more recent insertions, while higher values indicate more mutations, suggesting older insertions. The shape and number of peaks can indicate the activity patterns of different repeat families over time. You may look up more information on the different repeat families and how to analyse these results.
