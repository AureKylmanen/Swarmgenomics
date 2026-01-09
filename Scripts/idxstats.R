args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript idxstats_plot.R <input_csv> <read_length> <top_scaffolds> <palette> [custom_colors]")
}

input_file    <- args[1]
read_length   <- as.numeric(args[2])
top_n         <- as.numeric(args[3])
palette_name  <- args[4]
custom_colors <- if (length(args) >= 5 && args[5] != "") strsplit(args[5], ",")[[1]] else NULL

library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)

# Load and filter data
data <- read.csv(input_file)

# Prepare data with unmapped reads and coverage
top_chr <- data %>%
  arrange(desc(length)) %>%
  slice(1:top_n) %>%
  mutate(
    chrom = factor(chrom, levels = chrom),  # preserve order
    unmapped_per_mbp = unmapped / (length / 1e6),  # Unmapped reads per Mbp
    coverage = (mapped * read_length) / length,   # Coverage per chromosome
    mapped_percentage = mapped / (mapped + unmapped),   # Mapped percentage
    unmapped_percentage = unmapped / (mapped + unmapped)  # Unmapped percentage
  )

# Define a consistent color palette
if (!is.null(custom_colors)) {
  cols <- custom_colors
} else {
  cols <- brewer.pal(max(3, min(8, top_n)), palette_name)
}

# Plot 1: Chromosome Length Bar Plot
p1 <- ggplot(top_chr, aes(x = chrom, y = length / 1e6)) +
  geom_bar(stat = "identity", fill = cols[1], width = 0.7) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = paste("Top", top_n, "Scaffold Lengths"),
    x = NULL,
    y = "Length (Mbp)"
  ) +
  scale_y_continuous(expand = c(0, 0))


# Plot 2: Unmapped Reads per Mbp per Chromosome
p2 <- ggplot(top_chr, aes(x = chrom, y = unmapped_per_mbp)) +
  geom_bar(stat = "identity", fill = cols[2], width = 0.7) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Unmapped Reads per Mbp",
    x = NULL,
    y = "Unmapped Reads / Mbp"
  ) +
  scale_y_continuous(expand = c(0, 0))


# Plot 3: Coverage Bar Plot
p3 <- ggplot(top_chr, aes(x = chrom, y = coverage)) +
  geom_bar(stat = "identity", fill = cols[4], width = 0.7) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = paste("Estimated Coverage (Read Length:", read_length, "bp)"),
    x = "Chromosome",
    y = "Coverage (X)"
  ) +
  scale_y_continuous(expand = c(0, 0))


# Save each plot individually
ggsave("idxstats_scaffold_lengths.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("idxstats_unmapped_reads_per_mbp.png", plot = p2, width = 10, height = 6, dpi = 300)
ggsave("idxstats_coverage.png", plot = p3, width = 10, height = 6, dpi = 300)
