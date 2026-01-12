### NUMT visualization
### numt_visual.04.R
### 2025/09/26
### Visualize NUMT on mitochondrial genome
### v3 2024/11/22
### Final version
### v2 2024/11/15
### Different packages
### v1 2024/11/15
### Visualize NUMT on mito genome
### @author:Chen Yu Chi

# BiocManager::install("Biostrings")
# BiocManager::install("ggbio")
library(data.table)
library(Biostrings)
library(dplyr)
library(stringr)
library(circlize)

dir_work <- '.'
setwd(dir_work)

# numt
fa_file <- list.files(".", pattern = "\\.sep\\.fa$", full.names = TRUE)
numt <- fread(fa_file, header = FALSE)
numt <- numt[grepl("^>", V1), .(id = sub("^>", "", V1))]
setorder(numt, id)

# Read BLAST output and remove comment lines
blast_file <- list.files(".", pattern = "\\.blast\\.out$", full.names = TRUE)
blast_lines <- readLines(blast_file)
blast_clean <- blast_lines[!grepl("^#", blast_lines)]
blast <- fread(text = blast_clean, header = FALSE)

# Vectorized normalization and location ID
blast_chk <- blast[, .(
  loc = paste0(V2, ":", pmin(V9, V10), "-", pmax(V9, V10)),
  mt_start = V7,
  mt_end   = V8
)]
blast_chk <- unique(blast_chk, by = "loc")

# Merge with NUMTs
numt_loc <- merge(numt, blast_chk, by.x = "id", by.y = "loc")
colnames(numt_loc) <- c("numt", "mt_start", "mt_end")
numt_loc <- numt_loc %>%
  mutate(
    mt_start = pmin(mt_start, mt_end),
    mt_end   = pmax(mt_start, mt_end)
  ) %>%
  arrange(mt_start)

# mt annotation
mt_files <- list.files(path = ".", pattern = "_sequence\\.fasta$", full.names = TRUE, recursive = TRUE)
seq_lengths <- sapply(mt_files, function(f) sum(width(readDNAStringSet(f))))
mt_file <- mt_files[which.max(seq_lengths)]
mt_seq <- readDNAStringSet(mt_file)
len <- sum(width(mt_seq))

# Get all CSV files recursively in mitogenome directory
# Search recursively for CSV files in mitogenome folder
mito_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
csv_candidates <- unlist(lapply(mito_dirs, function(dir) {
  list.files(path = dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
}))

# Keep only CSVs with a 'details' column
valid_csvs <- Filter(function(f) {
  tryCatch({
    # Use tab as separator
    df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    "details" %in% colnames(df)
  }, error = function(e) FALSE)
}, csv_candidates)

if(length(valid_csvs) == 0){
  stop("❌ No valid CSV with 'details' column found in mitogenome folder")
}

# Pick the CSV with most rows (optional heuristic)
mt_info <- read.table(valid_csvs[which.max(sapply(valid_csvs, function(f) {
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  nrow(df)
}))], header = TRUE, sep = "\t", stringsAsFactors = FALSE)

mt_region <- mt_info$details


# Parse annotation regions
if (length(mt_region) > 1) {
  # Count the number of matches per line
  n_matches <- sapply(mt_region, function(line) {
    nrow(str_match_all(line, "([a-zA-Z0-9]+)\\((\\d+)-(\\d+)")[[1]])
  })
  
  # Select the line with the most matches
  mt_region <- mt_region[which.max(n_matches)]
}

# Parse annotation regions
pattern <- "([a-zA-Z0-9]+)\\((\\d+)-(\\d+)"
matches <- str_match_all(mt_region, pattern)[[1]]

genes  <- matches[, 2]
starts <- as.numeric(matches[, 3])
ends   <- as.numeric(matches[, 4])

mt_annot <- data.frame(gene = genes, start = starts, end = ends)

# Initialize mt_start and mt_end
mt_start <- 1   # mitochondrial coordinates usually start at 1
mt_end   <- len

# Create gaps safely
gap_starts <- c(mt_start, mt_annot$end + 1)
gap_ends   <- c(mt_annot$start - 1, mt_end)

# Only keep valid gaps
valid_idx <- which(gap_starts <= gap_ends)

gaps <- data.frame(
  gene  = paste0("intron", seq_along(valid_idx)),
  start = gap_starts[valid_idx],
  end   = gap_ends[valid_idx]
)


# Remove invalid gaps
gaps <- gaps[gaps$start < gaps$end, ]


combined <- rbind(
  mt_annot,
  data.frame(gene = gaps$gene, start = gaps$start, end = gaps$end)
)


# Example categories
# Define categories and colors
combined$category <- ifelse(grepl("^intron", combined$gene), "Noncoding",
                            ifelse(grepl("^rrn", combined$gene), "rRNA",
                                   ifelse(grepl("^nad|cox|cytb|atp", combined$gene), "Protein-Coding", "Other")))

category_colors <- c("Protein-Coding" = "aquamarine3", 
                     "rRNA" = "skyblue", 
                     "Noncoding" = "gray", 
                     "Other" = "pink")

combined$color <- category_colors[combined$category]
combined <- combined[order(combined$start), ]


# Initialize parameters for spacing
# Save plot directly in RESULTS_DIR
results_dir <- file.path("results")
if(!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
png(file.path(results_dir, "numt.png"), width = 4000, height = 4000, res = 300)

circos.clear()
circos.par(cell.padding = c(0, 0, 0, 0), 
           track.margin = c(0.05, 0.05))  # Adjust track margin for spacing

# Initialize circos plot for gene data
circos.initialize(factors = combined$gene, xlim = cbind(combined$start, combined$end))

# Gene track with category colors (first track)
circos.trackPlotRegion(
  ylim = c(0, 1),  # Y-limits for the gene track
  bg.border = NA,  # No border for the background
  track.height = 0.2,  # Reduce track height for gene track
  panel.fun = function(x, y) {
    sector_index <- get.cell.meta.data("sector.index")
    color <- combined$color[combined$gene == sector_index]
    
    # Draw gene rectangles
    circos.rect(combined$start[combined$gene == sector_index],
                ybottom = 0,
                combined$end[combined$gene == sector_index],
                ytop = 1,
                col = color,
                border = NA)
    
    # Add gene labels for non-intron genes
    if (!grepl("^intron", sector_index)) {
      mid_point <- (combined$start[combined$gene == sector_index] + 
                      combined$end[combined$gene == sector_index]) / 2
      circos.text(mid_point, 1.2, sector_index, cex = 2, col = "black", 
                  facing = "outside", niceFacing = TRUE, font = 2)
    }
  }
)

n = length(numt_loc$numt)
# NUMT regions track (second track)
circos.trackPlotRegion(
  ylim = c(0.5, n + 0.5),  # Y-limits for the NUMT track (space for the track)
  bg.border = NA,  
  track.height = 0.5,  
  panel.fun = function(x, y) {
    numt_sectors <- numt_loc$numt
    rectangle_height <- 0.5
    
    # Numt as sector
    for (i in seq_along(numt_sectors)) {
      # Get the corresponding NUMT data
      start_pos <- numt_loc$mt_start[i]
      end_pos <- numt_loc$mt_end[i]

      # rectangular method
      circos.rect(start_pos,
                  ytop = n - i + 1 + rectangle_height / 2,
                  end_pos,
                  ybottom = n - i + 1 - rectangle_height / 2, col = "black")
      # lines method
      circos.lines(c(start_pos, end_pos),
                   c(n - i + 1, n - i + 1),
                   col = "black",
                   lwd = 0.8)
      
    }
  }
)

# Clear the plot at the end
circos.clear()

dev.off()
