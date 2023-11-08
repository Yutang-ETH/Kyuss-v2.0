# install package and load the library
# install.packages("circlize")
library(circlize)

# sed the working directory
setwd("P:/Yutangchen/Kyuss_data/kyuss_nextpolish/circos")

# import chromosome information and make the required input file for circlize
myfai <- read.table("kyuss.nextdenovo.juicer.chr.fasta.fai", header = F, stringsAsFactors = F)

mychrinfo <- cbind.data.frame(myfai$V1, rep(0, nrow(myfai)), myfai$V2)
colnames(mychrinfo) <- c("chr", "start", "end")
mychrinfo$chr <- factor(mychrinfo$chr, levels = mychrinfo$chr)
head(mychrinfo)

# import gene and repeat 
gene <- read.table("gene.bed", header = F, stringsAsFactors = F)
gene$V4 <- floor((gene$V3 + gene$V2)/2)

REPEAT <- read.table("repeat.bed", header = F, stringsAsFactors = F)
REPEAT$V4 <- floor((REPEAT$V3 + REPEAT$V2)/2)
REPEAT$V5 <- floor(REPEAT$V3 - REPEAT$V2)

# GC content
GC<- read.table("kyuss_v2_chr.gc", header = F, stringsAsFactors = F)
GC$V1 <- sub("_sliding:.*", "", GC$V1)

# short read alignment depth
depth <- read.table("chr.1Mb.depth.bed", header = F, stringsAsFactors = F)

#-------------------------------------------------------------------------------------------------------#

# make the circos plot

# add gap between chromosomes
circos.par(gap.after = c(rep(3, 6), 15), 
           track.margin = c(0.05, 0.05),
           start.degree = 90)

# initialize the circos plot
circos.initialize(mychrinfo$chr, xlim = mychrinfo[, 2:3])

# first track, highlight each assembly
circos.track(sector = mychrinfo$chr, ylim = c(0, 1), cell.padding = c(0, 0, 0, 0), track.height = mm_h(1),
             bg.col = rep("gray90", 7),
             bg.border = "black" ) # c(rep("gray60", 7), rep("tomato", 7), rep("lightblue", 7)))

for(i in 1:nrow(mychrinfo)){
  
  circos.text((mychrinfo[i, 3] - mychrinfo[i, 2])/2, mm_y(5), 
              labels = gsub(".*chr", "", levels(mychrinfo$chr)[i]), 
              sector.index = levels(mychrinfo$chr)[i],
              track.index = 1,
              facing = "bending.inside",
              cex = 1.2)
  
  circos.axis(h = "top", 
              major.at = seq(0, mychrinfo[i, 3] + 50*10^6, by = 50*10^6),
              labels = seq(0, floor((mychrinfo[i, 3] + 50*10^6)/10^6), by = 50),
              sector.index = levels(mychrinfo$chr)[i],
              track.index = 1,
              labels.cex = 0.5,
              lwd = 1,
              minor.ticks = 0,
              major.tick.length = mm_y(0.2),
              direction = "outside",
              labels.facing = "outside",
              labels.pos.adjust = F)
              # col = c(rep("gray80", 7), rep("tomato", 7), rep("lightblue", 7))[i])

}


# second track, distribution of gene, heat map or hist gram

circos.track(sector = mychrinfo$chr, ylim = c(0, 100), cell.padding = c(0, 0, 0, 0), track.height = mm_h(4), bg.border = NA)

for(i in 1:nrow(mychrinfo)){
  
  # 10^6, window size 1 Mb
  mywindow <- seq(mychrinfo$start[i], mychrinfo$end[i], by = 10^6)
  
  # count how many genes in every window
  mynumber <- c()
  for(j in 1:(length(mywindow)-1)){
    mynumber <- c(mynumber, sum(gene[gene$V1 == mychrinfo[i, 1], 4] >= mywindow[j] & gene[gene$V1 == mychrinfo[i, 1], 4] < mywindow[j+1]))
  }
  
  # plot bar
  circos.rect(xleft = mywindow[-1],
              xright = mywindow[-1],
              ybottom = rep(0, length(mynumber)),
              ytop = mynumber,
              col = "lightblue4",
              border = "lightblue4",
              track.index = 2,
              sector.index = levels(mychrinfo$chr)[i])

}

# add axis
circos.yaxis(
  side = "left",
  at = seq(0, 100, 25),
  labels = seq(0, 100, 25),
  labels.cex = 0.5,
  tick.length = convert_x(0.1, "mm", levels(mychrinfo$chr)[1], 2),
  sector.index = levels(mychrinfo$chr)[1])


# third track, distribution of repeat

circos.track(sector = mychrinfo$chr, ylim = c(0, 1), cell.padding = c(0, 0, 0, 0), track.height = mm_h(3), bg.border = NA)

# use a line to show the abundance of repeats

for(i in 1:nrow(mychrinfo)){
  
  # 10^6, window size 1 Mb
  mywindow <- seq(mychrinfo$start[i], mychrinfo$end[i], by = 10^6)
  mydf <- cbind.data.frame(mywindow[1:(length(mywindow)-1)], mywindow[2:length(mywindow)])
  names(mydf) <- c("win1", "win2")
  
  # count how many repeats in every window
  myfun <- function(myrepeat = myrepeat, x){
    return(round(sum(myrepeat[myrepeat[ ,4] >= x[1] & myrepeat[ ,4] < x[2], 5])/10^6, 2))
  }
  
  mypercent <- unlist(apply(mydf, 1, myfun, myrepeat = REPEAT[REPEAT$V1 == mychrinfo[i, 1], ]))
  
  # plot bar
  circos.rect(xleft = mywindow[-1],
              xright = mywindow[-1],
              ybottom = rep(0, length(mypercent)),
              ytop = mypercent,
              col = "gold2",
              border = "gold2",
              track.index = 3,
              sector.index = levels(mychrinfo$chr)[i])
 
}

# add axis
circos.yaxis(
  side = "left",
  at = seq(0, 1, 0.25),
  labels = seq(0, 1, 0.25),
  labels.cex = 0.5,
  tick.length = convert_x(0.1, "mm", levels(mychrinfo$chr)[1], 3),
  sector.index = levels(mychrinfo$chr)[1])

# fourth track, distribution of GC contene

circos.track(sector = mychrinfo$chr, ylim = c(0, 1), cell.padding = c(0, 0, 0, 0), track.height = mm_h(3), bg.border = NA)

for(i in 1:nrow(mychrinfo)){
  
  # 10^6, window size 1 Mb
  mywindow <- seq(mychrinfo$start[i], mychrinfo$end[i], by = 10^6)
  
  # count how many GC in every window
  mypercent <- round(GC[GC$V1 == mychrinfo[i, 1], 2]/100, 4)
  
  # plot line 
  circos.lines(x = mywindow,
               y = mypercent,
               sector.index = levels(mychrinfo$chr)[i],
               col = "violetred4",
               type = "s",
               lwd = 2)
}

# add axis
circos.yaxis(
  side = "left",
  at = seq(0, 1, 0.25),
  labels = seq(0, 1, 0.25),
  labels.cex = 0.5,
  tick.length = convert_x(0.1, "mm", levels(mychrinfo$chr)[1], 4),
  sector.index = levels(mychrinfo$chr)[1])

# fifth track, short read alignment coverage

circos.track(sector = mychrinfo$chr, ylim = c(0, 100), cell.padding = c(0, 0, 0, 0), track.height = mm_h(3), bg.border = NA)

for(i in 1:nrow(mychrinfo)){
  
  # 10^6, window size 1 Mb
  mywindow <- seq(mychrinfo$start[i], mychrinfo$end[i], by = 10^6)
  
  # count depth
  mydepth <- depth[depth$V1 == mychrinfo[i, 1], 4]
  
  # plot line 
  circos.lines(x = mywindow,
               y = mydepth,
               sector.index = levels(mychrinfo$chr)[i],
               col = "orange4",
               type = "s",
               lwd = 2)
}

# add axis
circos.yaxis(
  side = "left",
  at = seq(0, 100, 25),
  labels = seq(0, 100, 25),
  labels.cex = 0.5,
  tick.length = convert_x(0.1, "mm", levels(mychrinfo$chr)[1], 5),
  sector.index = levels(mychrinfo$chr)[1])

circos.clear()

# add labels to each track
text(0, 0.95, labels = "A", pos = 2, offset = 1.5)
text(0, 0.83, labels = "B", pos = 2, offset = 1.5)
text(0, 0.67, labels = "C", pos = 2, offset = 1.5)
text(0, 0.53, labels = "D", pos = 2, offset = 1.5)
text(0, 0.38, labels = "E", pos = 2, offset = 1.5)


text(-0.2, 0.2, labels = "A: chromosome", cex = 0.8, font = 2, pos = 4)
text(-0.2, 0.1, labels = "B: gene density", cex = 0.8, font = 2, pos = 4)
text(-0.2, 0, labels = "C: repeat density", cex = 0.8, font = 2, pos = 4)
text(-0.2, -0.1, labels = "D: GC content", cex = 0.8, font = 2, pos = 4)
text(-0.2, -0.2, labels = "E: short-read alignment depth", cex = 0.8, font = 2, pos = 4)
