# ============================
# Libraries
# ============================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(magrittr)
})

# ============================
# Arguments
# ============================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript heterozygosity_plot.R <het_file> <top_n> <bar_color> <out_file>")
}

het_file  <- args[1]
top_n     <- as.numeric(args[2])
bar_color <- args[3]
out_file  <- args[4]

# ============================
# Load data
# ============================
data <- read.delim(het_file, sep = "\t")

# Remove rows with missing values
data <- data %>%
  filter(!is.na(Length) & !is.na(Heterozygosity))

# Select top N longest scaffolds
top_scaffolds <- data %>%
  arrange(desc(Length)) %>%
  head(top_n)

# ============================
# Plot
# ============================
p <- ggplot(top_scaffolds, aes(x = Chromosome, y = Heterozygosity)) +
  geom_bar(stat = "identity", fill = bar_color, width = 0.7) +
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
  scale_x_discrete(guide = guide_axis(n.dodg_
