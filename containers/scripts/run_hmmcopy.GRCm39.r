# the scripts were modifed from code provided by HMMcopy author Daniel Lai (communication with Xiaosai Yao)

args = commandArgs(TRUE)

evalue=as.numeric(args[1]);  

library(HMMcopy)
rfile <- "input.wig"
gfile <- "/data/GRCm39_gc_500k.wig"
mfile <- "/data/GRCm39_map_40mer_500k.wig"
x1_untrimmed_reads <- wigsToRangedData(rfile,gfile,mfile)
corrected <- correctReadcount(x1_untrimmed_reads)
corrected <- corrected[!is.na(corrected$copy) & is.finite(corrected$copy)]
corrected$chr <- as.factor(corrected$chr)


param <- HMMsegment(corrected, getparam = TRUE)
param$e = evalue 
#param$e= 0.995; Use this line to change the e value in HMMcopy if necessary, un-comment this line and change the value.
evalue

segments <- HMMsegment(corrected, param = param )

chromosomes <- levels(corrected$chr) # List of chromosomes...

lengths <- corrected[, .(chrlen = max(end)), by = chr]

lengths <- lengths[c(1,2,4,5,6,8,7,9,12,10,13,14,15,11,16,17,18,19,20,3,21), ] # Just sorting
lengths$offset <- c(0, cumsum(as.numeric(lengths$chrlen))[-nrow(lengths)])

corrected$state <- segments$state
#corrected$start <- start(corrected)
#corrected$end <- end(corrected)


# Merging offset information back into original corrected data frame
data <- merge(corrected, lengths, by = "chr")

pdf ("single.fixed.y.pdf", height = 5, width = 30)
# Reproducing internals of plotSegment for maximum customization
cols <- stateCols() # Predefined colours

range <- quantile(data$copy, na.rm = TRUE, prob = c(0.01, 0.99)) # Optional range, tweak this as desired
plot(x = data$start + data$offset, y = data$copy,xaxt='n', xlab="Chromosome", ylab="CopyNumber", col = cols[data$state], pch = 20, ylim = c(-1,1) )
title(xlab="", mgp=c(0.5,1,0) )
abline(v = lengths$offset) # Separate chromosomes
#text(lengths$offset, 1, labels = lengths$chr, adj = -0.1) # Label chromosomes
text(lengths$offset, 1, labels = gsub("chr", "", lengths$chr), adj = -0.2) # Label chromosomes

# Inserting the segmentation lines if desired...
segs <- merge(segments$segs, lengths[, c(1, 3)], all.x = TRUE)
segments(x0 = segs$start + segs$offset, y0 = segs$median, x1 = segs$end + segs$offset, col = "green", lwd = 3) # Note, segments is a function... not the output, unfortunately naming...
dev.off();

#rangedDataToSeg(corrected, "out.seg", column = "copy",verbose = TRUE)
write.table(segments$segs, file="segments.txt", sep="\t")
n=matrix(c(as.character(corrected$space),corrected$copy),ncol=2)
write.table(n, file="data.txt", sep="\t",row.names=FALSE, col.names=FALSE)
