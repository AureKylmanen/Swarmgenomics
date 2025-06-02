# Load required libraries 
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(readr)
library(RColorBrewer)

# Avoid scientific notation for heterozygosity
options(scipen = 999)

# Input files
idxstats_file <- "idxstats_clean.csv"
het_file <- "heterozygosity_results.txt"  # Expecting Chromosome, Heterozygosity, Length
read_length <- 100

# Read data
idx <- read.csv(idxstats_file)
het <- read_tsv(het_file, col_types = cols())

# Sort chromosomes by descending length
het_sorted <- het %>%
  arrange(desc(Length)) %>%
  mutate(chrom = factor(Chromosome, levels = Chromosome))

# Merge heterozygosity info into idxstats
idx_merged <- idx %>%
  inner_join(het_sorted, by = c("chrom" = "Chromosome")) %>%
  arrange(desc(Length)) %>%
  mutate(
    chrom = factor(chrom, levels = chrom),
    unmapped_per_mbp = unmapped / (length / 1e6),
    coverage = (mapped * read_length) / length
  )

# Custom theme for all plots
custom_theme <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.tag = element_text(size = 16, face = "bold")
  )

# Define colors
cols <- brewer.pal(8, "Set2")

# Plot A: Chromosome Length
p1 <- ggplot(idx_merged, aes(x = chrom, y = length / 1e6)) +
  geom_bar(stat = "identity", fill = cols[1], width = 0.6) +
  labs(title = "Scaffold Length", x = "Scaffold", y = "Length (Mbp)", tag = "A") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  custom_theme

# Plot B: Unmapped Reads per Mbp
p2 <- ggplot(idx_merged, aes(x = chrom, y = unmapped_per_mbp)) +
  geom_bar(stat = "identity", fill = cols[2], width = 0.6) +
  labs(title = "Unmapped Reads per Mbp", x = "Scaffold", y = "Unmapped / Mbp", tag = "B") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  custom_theme

# Plot C: Coverage
p3 <- ggplot(idx_merged, aes(x = chrom, y = coverage)) +
  geom_bar(stat = "identity", fill = cols[4], width = 0.6) +
  labs(title = paste("Coverage (Read Length:", read_length, "bp)"), x = "Scaffold", y = "Coverage (X)", tag = "C") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  custom_theme

# Plot D: Heterozygosity
p4 <- ggplot(idx_merged, aes(x = chrom, y = Heterozygosity)) +
  geom_bar(stat = "identity", fill = cols[5], width = 0.6) +
  labs(title = "Heterozygosity per Scaffold", x = "Scaffold", y = "Heterozygosity", tag = "D") +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.00001), expand = expansion(mult = c(0, 0.1))) +
  custom_theme

# Arrange in 2x2 grid
final_plot <- grid.arrange(
  p1, p2, p3, p4,
  ncol = 2
)

# Save to file
ggsave("idxstats_summary.png", final_plot, width = 14, height = 10)

