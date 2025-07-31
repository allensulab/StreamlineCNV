library(HMMcopy)
library(IRanges)

rfile <- "input.wig"
gfile <- "/data/GRCh38_gc_500k.wig"
mfile <- "/data/GRCh38_map_40mer_500k.wig"
x1_untrimmed_reads <- wigsToRangedData(rfile,gfile,mfile)
corrected <- correctReadcount(x1_untrimmed_reads)
corrected$chr <- as.factor(corrected$chr)
segments <- HMMsegment(corrected)

# Find appropriate offset for each chr to plot all points on same cartesian plane
lengths <- corrected[, .(chrlen = max(end)), by = chr]
lengths <- lengths[c(1, 2, 3, 4, 5, 6, 7, 9, 10, 12, 11, 13, 14, 15, 16, 17, 18, 19, 21, 20, 24, 23, 8, 22), ] # Just sorting
lengths$offset <- c(0, cumsum(as.numeric(lengths$chrlen))[-nrow(lengths)])
corrected$state <- segments$state

# Mergig offset information back into original corrected data frame
data <- merge(corrected, lengths, by = "chr")


# Desired order
chr_order <- c(as.character(1:22), "X", "Y")

# Keep only chromosomes that exist in your data
chr_order_present <- chr_order[chr_order %in% data$chr]

# Subset and reorder by space using RangedData::subset
data <- data[chr_order_present]

# Reproducing internals of plotSegment for maximum customization

cols <- c("#32CD32", "#0909bc","#0909bc", "#0909bc", "#FF0000","#FF0000");
segs <- merge(segments$segs, lengths[, c(1, 3)], all.x = TRUE)


# Ensure chr column is character
segs$chr <- as.character(segs$chr)

# Turn into factor with correct levels, then order
segs$chr <- factor(segs$chr, levels = chr_order)
segs <- segs[order(segs$chr), ]


k <- 0
data$state <- NA
for (i in 1:nrow(segs)) {
    for (j in seq(segs$start[i], segs$end[i], by = 500000)) {
        k <- k + 1
	if (k > nrow(data)) break
        if (segs$median[i] < 0.195 && segs$median[i] > -0.13) {
            data$state[k] <- 3
        } else if (segs$median[i] >= 0.195) {
            data$state[k] <- 5
        } else {
            data$state[k] <- 1
        }

        if (data$chr[k] == "chrY" || data$chr[k]=="chrX") {
            data$state[k] <- 1
        }
    }
}


pdf ("single.fixed.recolor.pdf", height = 5, width = 30)
par (mar= c(5.1, 6.1, 4.1, 2.1), c(3.2, 1, 0))
plot(x = data$start + data$offset, y = data$copy, col = cols[data$state], pch = 20, ylim = c(-1.2,1.2), ylab="Relative Copy Number (log2)", xlab="",  cex.axis=2, cex.lab=2,xaxt='n')
abline(v = lengths$offset) # Separate chromosomes
text(lengths$offset, 1, labels = gsub("chr", "", lengths$chr), adj = -0.1, cex=2) # Label chromosomes


# Inserting the segmentation lines if desired...
segments(x0 = segs$start + segs$offset, y0 = segs$median, x1 = segs$end + segs$offset, col = "#0909bc", lwd = 3) # Note, segments is a function... not the output
dev.off();
