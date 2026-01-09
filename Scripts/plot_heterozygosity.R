#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ============================
# Arguments
# ============================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript plot_heterozygosity.R <het_file> <top_n> <bar_color> <out_file>")
}

het_file  <- args[1]
top_n     <- as.numeric(args[2])
bar_color <- args[3]
out_file  <- args[4]

# ============================
# Load data
# ============================
data <- read.delim(het_file, sep = "\t")

data <- data %>%
  filter(!is.na(Length), !is.na(Heterozygosity))

top_scaffolds <- data %>%
  arrange(desc(Length)) %>%
  head(top_n)

# ============================
# Plot
# ============================
p <- ggplot(top_scaffolds, aes(x = Chromosome, y = Heterozygosity)) +
  geom_bar(stat = "identity", fill = bar_color, width = 0.7) +
  labs(
    x = "Chromosome / Scaffold",
    y = "Heterozygosity",
    title = paste("Heterozygosity Across Top", top_n, "Longest Scaffolds")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_continuous(expand = c(0, 0))

ggsave(out_file, plot = p, width = 10, height = 6, dpi = 300)
