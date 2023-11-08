# plot hic matrix using ggplot2
# setwd
setwd("P:/Yutangchen/Kyuss_data/kyuss_nextpolish/hicexplorer")

# import matrix data and bin bed file
hicm <- read.table("hic/hic_pro_matrix", header = F, stringsAsFactors = F)
binbed <- read.table("hic/dpnii.bed", header = F, stringsAsFactors = F)

# select chr
binbed_chr <- binbed[binbed$V1 %in% paste("chr", 1:7, sep = ""), ]

# select interchr matrix 
hicm_chr <- hicm[hicm$V1 <= max(binbed_chr$V4) & hicm$V2 <= max(binbed_chr$V4), ]

# use the midpoint of each bin as the coordinate
binbed_chr$V5 <- (binbed_chr$V2 + binbed_chr$V3)/2

# add chr, pos information from binbed table to hicm table
# binx is either colum 1 or 2 in binbed_chr
findchr <- function(x, df, binx){
  return(df[df$V4 == x[binx], 1])
}

findpos <- function(x, df, binx){
  return(df[df$V4 == x[binx], 5])
}

# find corresponding chr for bins in the first two columns
bin1_chr <- unlist(apply(hicm_chr, 1, findchr, df = binbed_chr, binx = 1))
bin2_chr <- unlist(apply(hicm_chr, 1, findchr, df = binbed_chr, binx = 2))

# find corresonding pos for bins in the first two columns
bin1_pos <- unlist(apply(hicm_chr, 1, findpos, df = binbed_chr, binx = 1))
bin2_pos <- unlist(apply(hicm_chr, 1, findpos, df = binbed_chr, binx = 2))

hicm_chr$V4 <- bin1_chr
hicm_chr$V5 <- bin2_chr
hicm_chr$V6 <- bin1_pos
hicm_chr$V7 <- bin2_pos

colnames(hicm_chr) <- c("bin1", "bin2", "count", "chrbin1", "chrbin2", "chrpos1", "chrpos2")

# delete these vectors
rm(list = c("bin1_chr", "bin2_chr", "bin1_pos", "bin2_pos"))

# make a full matrix
hicm_chr_r <- cbind.data.frame(hicm_chr$bin2,
                               hicm_chr$bin1,
                               hicm_chr$count,
                               hicm_chr$chrbin2,
                               hicm_chr$chrbin1,
                               hicm_chr$chrpos2,
                               hicm_chr$chrpos1)

colnames(hicm_chr_r) <- colnames(hicm_chr)

hicm_chr_full <- rbind.data.frame(hicm_chr, hicm_chr_r)

write.table(hicm_chr_full, "hicm_chr_full.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# draw the HiC contact matrix using ggplot2
library(ggplot2)
library(scales)

# inter chr interaction plot
png('interchr_log2.png', 1500, 1500)
ggplot(data = hicm_chr_full, aes(x = chrpos1, y = chrpos2, color = log2(count))) + 
         geom_point(pch = 15) +
         facet_grid(factor(chrbin2, levels = c("chr7", "chr6", "chr5", "chr4", "chr3", "chr2", "chr1")) ~ chrbin1, scales="free", space="free") +
         labs(x = "", y = "") +
         scale_colour_gradient2(low = "white", high = "red") +
         theme(axis.line = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position='none',
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              strip.text = element_text(face = "bold", size = rel(1.5)),
              strip.background = element_blank(),
              panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid")) +
         scale_y_continuous(expand = c(0.02,0.02)) + 
         scale_x_continuous(expand = c(0.02,0.02))        
dev.off()         

