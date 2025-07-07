### NUMT visualization
### numt_visual.03.R
### v3 2024/11/22
### Final version
### v2 2024/11/15
### Different packages
### v1 2024/11/15
### Visualize NUMT on mito genome
### @author:Chen Yu Chi

library(Biostrings)
library(circlize)
library(stringr)

dir_work <- '/path/to/working/dir/'
setwd(dir_work)

# mito genes annotations
mt_seq <- readDNAStringSet('**.path_sequence.fasta')
len <- width(mt_seq)

mt_info <- read.table('*.csv', header = T, stringsAsFactors = F)
mt_region <- mt_info$details

numt_loc <- read.table('mt_loc.lst', header = F, stringsAsFactors = F)
colnames(numt_loc) <- c('numt', 'mt_start', 'mt_end')
numt_loc$mt_start <- ifelse(numt_loc$mt_start > numt_loc$mt_end, 
                            numt_loc$mt_end, numt_loc$mt_start)
numt_loc$mt_end <- ifelse(numt_loc$mt_start > numt_loc$mt_end, 
                          numt_loc$mt_start, numt_loc$mt_end)
numt_loc <- numt_loc[order(numt_loc$mt_start), ]

# Get genes, starts and ends
pattern <- "([a-zA-Z0-9]+)\\((\\d+)-(\\d+)"
matches <- str_match_all(mt_region, pattern)

genes <- matches[[1]][, 2]
starts <- as.numeric(matches[[1]][, 3])
ends <- as.numeric(matches[[1]][, 4])

mt_annot <- data.frame(gene = genes, start = starts, end = ends)
#mt_annot$color <- viridis(n = nrow(mt_annot))



# Example of circular plot
mt_start <- 0
mt_end <- len

gaps <- data.frame(
  gene = paste0("intron", seq_len(nrow(mt_annot) + 1)),
  start = c(mt_start, mt_annot$end + 1),
  end = c(mt_annot$start - 1, mt_end)
)

# filter out non-exsit gaps
gaps <- gaps[gaps$start <= gaps$end, ]

# Combine genes and gaps
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
png("circos_plot.line.png", width = 1200, height = 1200, res = 150)
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
      circos.text(mid_point, 1.2, sector_index, cex = 0.8, col = "black", 
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
      # circos.rect(start_pos, 
      #             ytop = n - i + 1 + rectangle_height / 2, 
      #             end_pos,
      #             ybottom = n - i + 1 - rectangle_height / 2, col = "black")
      # lines method
      circos.lines(c(start_pos, end_pos),
                   c(n - i + 1, n - i + 1),
                   col = "black",
                   lwd = 1.5)
      
    }
  }
)

# Clear the plot at the end
circos.clear()

dev.off()
