# PSMC 

Pairwise Sequentially Markovian Coalescent (PSMC) model is a method used in genetics to study the history of populations. It analyses the genetic variation within the DNA of a single individual to infer past population sizes over time. It is a way to look back in time and see how the number of ancestors of a species has changed, which can give insights into events like population bottlenecks or expansions.

Here we use software package from https://github.com/lh3/psmc

## Downloading the software

### PSMC

In your software directory (e.g. /vol/storage/software)
```
#clone the repository
git clone https://github.com/lh3/psmc.git

#make psmc
cd psmc
make

#make utils
cd utils
make
```
### PSMC plot
In your psmc directory (e.g. /vol/storage/software/psmc)
```
#clone the repository
git clone https://github.com/willyrv/ms-PSMC.git

#change permissions
chmod +x ./*
```

# Running PSMC

#### Convert fastq files to PSMC fasta format
```
/path/to/psmc/utils/fq2psmcfa -q20 /working_dir/diploid.fq.gz > /working_dir/diploid.psmcfa
```
#### Run the PSMC analysis
You may change the parameters according to your preferences. This step will take few hours, so you may wish to use "nohup" and "&".
```
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o diploid.psmc diploid.psmcfa

#example
nohup /vol/storage/psmc_plot_dir/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o /vol/storage/swarmGenomics/golden_eagle/diploid.psmc /vol/storage/swarmGenomics/golden_eagle/diploid.psmcfa &
```
#### Convert the output of PSMC analysis
```
utils/psmc2history.pl diploid.psmc | utils/history2ms.pl > ms-cmd.sh

#Example
/vol/storage/psmc_plot_dir/psmc/utils/psmc2history.pl /vol/storage/swarmGenomics/golden_eagle/diploid.psmc | /vol/storage/psmc_plot_dir/psmc/utils/history2ms.pl > /vol/storage/swarmGenomics/golden_eagle/ms-cmd.sh
```

# Plot the results
```
#Plot as pdf
utils/psmc_plot.pl diploid diploid.psmc

#Example
/vol/storage/psmc_plot_dir/psmc/utils/psmc_plot.pl -p /vol/storage/swarmGenomics/golden_eagle/diploid diploid.psmc

#Plot as png
#copy plot_results.py to working directory

cp plot_results.py ../golden_eagle/

python plot_results.py diploid diploid.psmc
```
