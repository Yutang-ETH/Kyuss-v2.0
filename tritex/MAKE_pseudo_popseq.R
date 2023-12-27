
args = commandArgs(trailingOnly=TRUE)

# args[1] the single copy bed
# args[2] fai file of the guide map
# args[3] prefix of the output file

# setwd("P:/Yutangchen/Rabiosa_canu/hic")
# make my own pseudo-pop-seq table
# first import the single-copy bed file
mybed <- read.table(args[1], header = F, stringsAsFactors = F)
myfai <- read.table(args[2], header = F, stringsAsFactors = F)

css_contig <- mybed$V4
popseq_cM <- mybed$V2/10^6
sorted_genome <- rep("", dim(mybed)[1])
css_contig_length <- mybed$V3 - mybed$V2
popseq_alphachr <- paste(sub("chr", "", mybed$V1), sorted_genome, sep = "")
sorted_arm <- rep(NA, dim(mybed)[1])
popseq_chr <- as.numeric(sub("chr", "", mybed$V1))
sorted_alphachr <- popseq_alphachr
sorted_chr <- popseq_chr
sorted_lib <- popseq_alphachr
flip <- rep(FALSE, dim(mybed)[1])

# sorted_genome <- popseq_chr

mymax <- c()
for(i in 1:length(unique(mybed$V1))){
  mymax <- c(mymax, rep(myfai[myfai$V1 == unique(mybed$V1)[i], 2]/10^6, length(mybed$V1[mybed$V1 == unique(mybed$V1)[i]])))
}

mb <- popseq_cM

Rabiosa_popseq <- cbind.data.frame(css_contig, popseq_cM, sorted_genome, css_contig_length,
                                   popseq_alphachr, sorted_arm, popseq_chr, sorted_alphachr,
                                   sorted_chr, sorted_lib, flip, mymax, mb)
colnames(Rabiosa_popseq)[12] <- "max"

Rabiosa_popseq <- Rabiosa_popseq[order(Rabiosa_popseq$popseq_chr), ]
saveRDS(Rabiosa_popseq, file = paste(args[3], "_popseq.Rds", sep = ""))

# import the save Rabiosa_popseq file and check
# readRDS('Rabiosa_popseq.Rds')->mypopseq
