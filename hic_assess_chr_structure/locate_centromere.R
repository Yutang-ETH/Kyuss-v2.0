# estimate the centromere position based on the inter-chromosome contact

# setwd
setwd("P:/Yutangchen/Kyuss_data/kyuss_nextpolish/hicexplorer")

library(ggplot2)

hicm_chr <- read.table('hicm_chr_full.txt', header = T, stringsAsFactors = F)

####
ggplot(xx, aes(x = chrpos1, y = count)) + 
  geom_point(pch = 15) +
  facet_grid(rows = vars(chrbin2),scales="free", space="free") 

test <- inter_chr1[inter_chr1$bin1 == 103, ]

centrox <- NULL
for(i in paste('chr', 1:7, sep = '')){
  
  inter_chrx <- hicm_chr[hicm_chr$chrbin1 == i & hicm_chr$chrbin2 != i, ]
  
  # get bins with the highest count
  binx <- c()
  for(x in unique(inter_chrx$chrbin2)){
    
    xx <- inter_chrx[inter_chrx$chrbin2 == x, ]
    xx <- xx[xx$chrpos1 >= 20000000 & xx$chrpos1 <= max(xx$chrpos1) - 20000000, ]
    binx <- c(binx, xx[which.max(xx$count), 1])
    
  }
  
  # try to rule out noise with the distribution of the bins
  # the most frequent bins should be the bins in centromere
  yy <- hist(binx, breaks = seq(min(binx), max(binx) + 20, by = 20), plot = F)
  
  # find bins that are most frequent
  zz <- binx[binx >= yy$breaks[which.max(yy$counts)] & binx <= yy$breaks[which.max(yy$counts)] + 20]
  
  # now this is the lead bin, which is the representative of centromere
  lead <- sort(zz)[ceiling(length(zz) / 2 + 0.1)]
  
  # now need to define a region for centromere
  
  for(x in unique(inter_chrx$chrbin2)){
    
    xx <- inter_chrx[inter_chrx$chrbin2 == x, ]
    xx0 <- which(xx$count - quantile(xx$count)[4] <= 0)
    leadrownumber <- which(xx$count == max(xx[xx$bin1 == lead, 3]))
    
    centrox <- rbind.data.frame(centrox, xx[tail(xx0[which(xx0 <= leadrownumber)], 1) : head(xx0[which(xx0 >= leadrownumber)], 1), ])
    
  }
  
  
}

mycentro <- NULL
for(x in unique(centrox$chrbin2)){
  
  xx <- centrox[centrox$chrbin2 == x, ]
  xx <- xx[order(xx$bin2), ]
  xx <- xx[!duplicated(xx$bin2), ]
  mycentro <- rbind.data.frame(mycentro, xx)
}

output_centro <- NULL
for(x in sort(unique(mycentro$chrbin2))){
  
  xx <- mycentro[mycentro$chrbin2 == x, ]
  bed <- data.frame(chr = x, start = min(xx$chrpos2), end = max(xx$chrpos2))
  output_centro <- rbind.data.frame(output_centro, bed)
  
}

write.table(output_centro, 'kyuss_v2_chr_centro.bed', quote = F, sep = '\t', col.names = F,row.names = F)
