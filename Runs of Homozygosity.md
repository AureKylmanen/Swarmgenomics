# Runs of Homozygosity
Runs of Homozygosity (ROH) are continuous stretches of homozygous genotypes within an individual's genome, indicating that the segments are inherited from common ancestors. ROH can provide insights into the genetic history of populations, including levels of inbreeding, past population bottlenecks, and patterns of natural selection.
## Installations
You should already have bcftools installed, otherwise you will only need to download the roh_plot.py script for plotting.

## ROH analysis

### Preparations
```
# Make a directory for the roh results
mkdir roh
```

### Running the script
```
# Change directory to vcf
cd vcf

# Execute script
for file in *
do 
echo $file
/vol/storage/bcftools-1.19/bcftools roh -G30 --AF-dflt 0.4 $file > /vol/storage/swarmGenomics/golden_eagle/roh/$file.roh_chr.txt
python /vol/storage/roh_plot_dir/roh_plot.py $file.roh.png /vol/storage/swarmGenomics/golden_eagle/roh/$file.roh_chr.txt
done
```

### Output
The output file is a density plot (roh.png) to visualise the distribution of ROH lengths in the genome. The x-axis represents the length of ROH fragments and the y-axis represents the density of ROH fragments. Higher peaks in the density plot indicate lengths where ROH fragments are more frequently observed.
