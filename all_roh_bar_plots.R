# Load necessary libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)

# Step 1: Load the data
roh_data <- read.table("RG.txt", header = FALSE, comment.char = "#", sep = "\t")

# Assign column names (assuming the 6th column is Length_bp)
colnames(roh_data) <- c("Type", "Sample", "Chromosome", "Start", "End", "Length_bp", "Num_markers", "Quality")

# Step 2: Convert Length (bp) to Length (Mb)
roh_data$Length_mb <- roh_data$Length_bp / 1e6

# Step 3: Categorize ROHs into new size classes
bins <- c(0.01, 0.1, 1, 3, Inf)  # Updated breakpoints to start at 0.01
labels <- c("10kbp - 0.1 Mbp", "0.1 - 1 Mbp", "1 - 3 Mbp", ">3 Mbp")  # Updated categories
roh_data$Size_category <- cut(roh_data$Length_mb, breaks = bins, labels = labels, right=FALSE) 

# Step 4: Aggregate data by size category (ensuring all categories exist)
all_categories <- c("10kbp - 0.1 Mbp", "0.1 - 1 Mbp", "1 - 3 Mbp", ">3 Mbp")  # Define all categories

summary <- roh_data %>%
    group_by(Size_category) %>%
    summarize(
        Total_count = n(),
        Total_length_mb = sum(Length_mb)
    ) %>%
    ungroup() %>%
    complete(Size_category = factor(all_categories, levels = all_categories),  
             fill = list(Total_count = 0, Total_length_mb = 0)) %>%
    filter(!is.na(Size_category))  # Remove NA rows

# Define a common theme with increased font sizes and axis lines using linewidth
custom_theme <- theme_minimal() +
    theme(
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        axis.line = element_line(linewidth=0.5), # Add axis lines with thickness of 0.8 using linewidth
        axis.text.x = element_text(angle=45, hjust=1, size=14), # Rotate x-axis labels and increase font size
        axis.text.y = element_text(size=14),                   # Increase y-axis tick label font size
        axis.title.x = element_text(size=16),                  # Increase x-axis title font size
        axis.title.y = element_text(size=16),                  # Increase y-axis title font size
        plot.title = element_text(hjust=0.5, size=18)          # Center-align title and increase its font size
    )

# Bar plot for number of ROHs with custom theme applied and dynamic y-axis adjustment
count_plot <- ggplot(summary, aes(x=Size_category, y=Total_count)) +
    geom_bar(stat="identity", fill="#48c9b0") + 
    scale_y_continuous(expand=expansion(mult=c(0, 0.1))) +   # Automatically add ~10% space above bars dynamically 
    labs(title="Number of RoHs",
         x="RoH Size Category",
         y="Number of RoHs") +
    custom_theme

# Bar plot for total length of ROHs in Mb with custom theme applied and dynamic y-axis adjustment
length_plot <- ggplot(summary, aes(x=Size_category, y=Total_length_mb)) +
    geom_bar(stat="identity", fill="#f5b041") + 
    scale_y_continuous(expand=expansion(mult=c(0, 0.1))) +   # Automatically add ~10% space above bars dynamically 
    labs(title="Total Length of RoHs (Mb)",
         x="RoH Size Category",
         y="Total Length (Mb)") +
    custom_theme

# Save plots as PNG file side by side using gridExtra package:
png("roh_bar_plots.png", width=1200, height=600)
grid.arrange(count_plot, length_plot, ncol=2) # Arrange two plots side by side
dev.off()
