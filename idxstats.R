args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript idxstats_plot.R <input_csv> <read_length>")
}

input_file <- args[1]
read_length <- as.numeric(args[2])

library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)

# Load and filter data
data <- read.csv(input_file)

# Prepare data with unmapped reads and coverage
top_chr <- data %>%
  arrange(desc(length)) %>%
  slice(1:20) %>%
  mutate(
    chrom = factor(chrom, levels = chrom),  # preserve order
    unmapped_per_mbp = unmapped / (length / 1e6),  # Unmapped reads per Mbp
    coverage = (mapped * read_length) / length,   # Coverage per chromosome
    mapped_percentage = mapped / (mapped + unmapped),   # Mapped percentage
    unmapped_percentage = unmapped / (mapped + unmapped)  # Unmapped percentage
  )

# Define a consistent color palette
cols <- brewer.pal(8, "Set2")  # distinctive & colorblind-friendly

# Plot 1: Chromosome Length Bar Plot
p1 <- ggplot(top_chr, aes(x = chrom, y = length / 1e6)) +  # Convert length to Mbp
  geom_bar(stat = "identity", fill = cols[1], width = 0.7) +
  theme_classic() +  # Using classic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black", linewidth = 0.5),  # Use linewidth instead of size
        plot.margin = margin(10, 10, 10, 10)) +  # Add margin to avoid clipping
  labs(title = "Top 20 Scaffold Lengths",
       x = NULL, y = "Length (Mbp)") +
  scale_y_continuous(expand = c(0, 0))

# Plot 2: Unmapped Reads per Mbp per Chromosome
p2 <- ggplot(top_chr, aes(x = chrom, y = unmapped_per_mbp)) +
  geom_bar(stat = "identity", fill = cols[2], width = 0.7) +
  theme_classic() +  # Using classic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black", linewidth = 0.5),  # Use linewidth instead of size
        plot.margin = margin(10, 10, 10, 10)) +  # Add margin to avoid clipping
  labs(title = "Unmapped Reads per Mbp",
       x = NULL, y = "Unmapped Reads / Mbp") +
  scale_y_continuous(expand = c(0, 0))

# Plot 3: Coverage Bar Plot
p3 <- ggplot(top_chr, aes(x = chrom, y = coverage)) +
  geom_bar(stat = "identity", fill = cols[4], width = 0.7) +
  theme_classic() +  # Using classic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black", linewidth = 0.5),  # Use linewidth instead of size
        plot.margin = margin(10, 10, 10, 10)) +  # Add margin to avoid clipping
  labs(title = paste("Estimated Coverage (Read Length:", read_length, "bp)"),
       x = "Chromosome", y = "Coverage (X)") +
  scale_y_continuous(expand = c(0, 0))

# Combine all plots into one PNG (stacked bar plot + other plots)
combined_plot <- p1 / p2 / p3 + plot_layout(ncol = 1)

# Save the combined plot as a PNG
ggsave("idxstats_summary.png", combined_plot, width = 12, height = 12)
