import pandas as pd
import matplotlib.pyplot as plt

# Load the heterozygosity data
data = pd.read_csv("heterozygosity_results.txt", sep="\t")

# Create the bar plot
plt.figure(figsize=(10, 6))
plt.bar(data['Chromosome'], data['Heterozygosity'], color='skyblue')
plt.xlabel('Chromosome/Scaffold')
plt.ylabel('Heterozygosity')
plt.title('Heterozygosity Across Chromosomes/Scaffolds')
plt.xticks(rotation=90)
plt.tight_layout()

# Save the plot as an image file
plt.savefig("heterozygosity_plot.png")

# Display the plot
plt.show()

