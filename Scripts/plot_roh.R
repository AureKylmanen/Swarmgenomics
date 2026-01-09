suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
})

# ============================
# Arguments
# ============================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: Rscript plot_roh.R <RG.txt> <bins> <labels> <count_color> <length_color> <out_png>")
}

rg_file      <- args[1]
bins         <- as.numeric(strsplit(args[2], ",")[[1]])
labels       <- strsplit(args[3], ",")[[1]]
count_color  <- args[4]
length_color <- args[5]
out_png      <- args[6]

# ============================
# Load data
# ============================
roh <- read.table(rg_file, header = FALSE, sep = "\t")

colnames(roh) <- c(
  "Type", "Sample", "Chromosome",
  "Start", "End", "Length_bp",
  "Num_markers", "Quality"
)

roh$Length_mb <- roh$Length_bp / 1e6

roh$Size_category <- cut(
  roh$Length_mb,
  breaks = bins,
  labels = labels,
  right = FALSE
)

summary <- roh %>%
  group_by(Size_category) %>%
  summarise(
    Total_count = n(),
    Total_length_mb = sum(Length_mb),
    .groups = "drop"
  ) %>%
  complete(
    Size_category = factor(labels, levels = labels),
    fill = list(Total_count = 0, Total_length_mb = 0)
  )

# ============================
# Theme
# ============================
custom_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18)
  )

# ============================
# Plots
# ============================
p_count <- ggplot(summary, aes(Size_category, Total_count)) +
  geom_bar(stat = "identity", fill = count_color) +
  labs(
    title = "Number of Runs of Homozygosity",
    x = "RoH Size Category",
    y = "Count"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  custom_theme

p_length <- ggplot(summary, aes(Size_category, Total_length_mb)) +
  geom_bar(stat = "identity", fill = length_color) +
  labs(
    title = "Total Length of Runs of Homozygosity",
    x = "RoH Size Category",
    y = "Total Length (Mb)"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  custom_theme

ggsave(
  out_png,
  grid.arrange(p_count, p_length, ncol = 2),
  width = 12,
  height = 6,
  dpi = 300
)
