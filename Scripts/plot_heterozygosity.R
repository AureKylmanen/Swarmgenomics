# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load the heterozygosity data
data <- read.delim("heterozygosity_results.txt", sep="\t")

# Check for missing values in Length and Heterozygosity columns and remove them
data <- data %>% filter(!is.na(Length) & !is.na(Heterozygosity))

# Get the top 20 longest scaffolds based on scaffold length
top_20_longest_scaffolds <- data %>%
  arrange(desc(Length)) %>%  # Sort by scaffold length in descending order
  head(20)  # Select the top 20 longest scaffolds

# Create the bar plot
ggplot(top_20_longest_scaffolds, aes(x = Chromosome, y = Heterozygosity)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +  # Adjust width to avoid floating
  labs(
    x = "Chromosome/Scaffold",
    y = "Heterozygosity",
    title = "Heterozygosity Across Top 20 Longest Scaffolds"
  ) +
  theme_classic() +  # White background
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Adjust label size and angle for clarity
    axis.ticks.x = element_line(linewidth = 0.5),  # Shorter tick marks
    panel.grid = element_blank(),  # Remove grid lines
    plot.title = element_text(hjust = 0.5)  # Center title
  ) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1)) +  # Adjust dodging for x-axis labels (1 dodge for clearer view)
  coord_cartesian(clip = 'off') +  # Ensure bars are plotted properly without clipping
  scale_y_continuous(expand = c(0, 0))  # Remove space between bars and x-axis

# Save the plot as an image file
ggsave("heterozygosity_top_20_longest_plot.png", width = 10, height = 6)
