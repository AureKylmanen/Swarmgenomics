args <- commandArgs(trailingOnly = TRUE)
het_file <- args[1]            # heterozygosity_results.txt
top_n <- as.numeric(args[2])   # e.g., 20
bar_color <- args[3]            # Color passed from params
out_file <- args[4]            # heterozygosity_plot.png

# Load the heterozygosity data
data <- read.delim(het_file, sep="\t")

# Remove rows with missing values
data <- data %>% filter(!is.na(Length) & !is.na(Heterozygosity))

# Select top N longest scaffolds
top_scaffolds <- data %>%
  arrange(desc(Length)) %>%
  head(top_n)

# Create bar plot
ggplot(top_scaffolds, aes(x = Chromosome, y = Heterozygosity)) +
  geom_bar(stat = "identity", fill = bar_color, width = 0.7) +   # Use variable here
  labs(
    x = "Chromosome/Scaffold",
    y = "Heterozygosity",
    title = paste("Heterozygosity Across Top", top_n, "Longest Scaffolds")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.ticks.x = element_line(linewidth = 0.5),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(expand = c(0, 0))

# Save plot
ggsave(out_file, width = 10, height = 6)
