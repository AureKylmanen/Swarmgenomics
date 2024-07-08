import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import glob

#Arguments: RoH-File and Name of Chromosome
#Example: python roh_density.py name.png /vol/storage/VulpesLagupos/results_with_repeats/roh_with_repeats.txt 


lst_file = glob.glob('/vol/storage/swarmGenomics/golden_eagle/roh/*.vcf.roh_chr.txt')

result=[]
for chr in lst_file:
    f=open(chr,"r")
    lines=f.readlines()
    counter=0
    for x in lines:
        counter+=1
        #only for the first chromosome
        if counter > 5 and x.split('\t')[0]=="RG":
            length=float(x.split('\t')[5].rstrip())
            result.append(length)
    f.close()




sns.kdeplot(result)
	#values, base = np.histogram(result, bins=4000)
	#cumulative = np.cumsum(values)
	#plt.plot(base[:-1], cumulative, c='blue')
	#sorted_data = np.sort(result)
	#plt.step(sorted_data, np.arange(sorted_data.size))  # From 0 to the number of data points-1
	#plt.step(sorted_data[::-1], np.arange(sorted_data.size))  # From the number of data points-1 to 0
    
plt.title('Density Plot')
plt.ylabel('Density of RoH fragments')
plt.xlabel('Length of RoH in bases')
#log scale
plt.xscale('log') 
plt.savefig('roh.png')



