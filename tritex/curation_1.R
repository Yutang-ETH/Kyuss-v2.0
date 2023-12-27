# working directory is"/home/yutachen/public/Yutangchen/Kyuss_data/tritex"

#load TRITEX R code
source('/home/yutachen/public/Yutangchen/Kyuss_data/tritex/bitbucket/pseudomolecule_construction.R')

# #read pseudo-POPSEQ data
# # The file is provided in the download folder: pseudo-POPSEQ 1 Mb = 1 pseudo-cM
readRDS('Rabiosa_popseq.Rds')->popseq
# convert data.frame to data table
setDT(popseq)

# # 
# # #read contig lengths
f <- 'kyuss/kyuss.nextpolish.asm.rename.fasta.fai'
fread(f, head=F, select=1:2, col.names=c("scaffold", "length"))->fai
# # 
# # #read "POPSEQ" (guide map) marker alignment
f <- 'kyuss/kyuss_rabiosa.paf.gz'
read_morexaln_minimap(paf=f, popseq=popseq, minqual=30, minlen=500, prefix=F)->morexaln
# # 
# # #read Hi-C ("fragment") pairs  
dir <- 'kyuss'
fread(paste('find', dir, '| grep "fragment_pairs.tsv.gz$" | xargs zcat'),
      header=F, col.names=c("scaffold1", "pos1", "scaffold2", "pos2"))->fpairs
# # 
# # #init and save assembly
init_assembly(fai=fai, cssaln=morexaln, fpairs=fpairs) -> assembly
anchor_scaffolds(assembly = assembly, popseq=popseq, species="lolium") -> assembly
add_hic_cov(assembly, binsize=1e4, binsize2=1e6, minNbin=50, innerDist=3e5, cores=20)->assembly
saveRDS(assembly, file="kyuss_assembly.Rds")
# # 
# # #create diagnostic plots for contigs longer than 1 Mb and putative chimeras with strong drops in Hi-C coverage
assembly$info[length >= 1e6, .(scaffold, length)][order(-length)] -> s
plot_chimeras(assembly=assembly, scaffolds=s,  species="lolium", refname="Rabiosa", autobreaks=F, mbscale=1,
              file="kyuss_assembly_1Mb.pdf", cores=30)
# # 
assembly$info[length >= 1e6 & mri <= -1, .(scaffold, length)][order(-length)] -> s
plot_chimeras(assembly=assembly, scaffolds=s,  species="lolium", refname="Rabiosa", autobreaks=F, mbscale=1,
              file="kyuss_assembly_chimeras.pdf", cores=30)

# # ll. 38 and 33 are identical so that the pages in the chimeras PDF can be referenced
assembly$info[length >= 1e6 & mri <= -1, .(scaffold, length)][order(-length)] -> ss
i=1
ss[i]$scaffold -> s
assembly$cov[s, on='scaffold'][bin >= 1e6 & bin <= 10e6][order(r)][1, .(scaffold, bin)] -> b

i=16
ss[i]$scaffold -> s
rbind(b, assembly$cov[s, on='scaffold'][bin >= 0.5e6 & bin <= 1.5e6][order(r)][1, .(scaffold, bin)])->b

i=16
ss[i]$scaffold -> s
rbind(b, assembly$cov[s, on='scaffold'][bin >= 2e6 & bin <= 3e6][order(r)][1, .(scaffold, bin)])->b

i=16
ss[i]$scaffold -> s
rbind(b, assembly$cov[s, on='scaffold'][bin >= 4e6 & bin <= 4.5e6][order(r)][1, .(scaffold, bin)])->b

i=17
ss[i]$scaffold -> s
rbind(b, assembly$cov[s, on='scaffold'][bin >= 0 & bin <= 1e6][order(r)][1, .(scaffold, bin)])->b

i=17
ss[i]$scaffold -> s
rbind(b, assembly$cov[s, on='scaffold'][bin >= 2e6 & bin <= 3e6][order(r)][1, .(scaffold, bin)])->b

i=17
ss[i]$scaffold -> s
rbind(b, assembly$cov[s, on='scaffold'][bin >= 4e6 & bin <= 4.3e6][order(r)][1, .(scaffold, bin)])->b

i=17
ss[i]$scaffold -> s
rbind(b, assembly$cov[s, on='scaffold'][bin >= 4.3e6 & bin <= 4.7e6][order(r)][1, .(scaffold, bin)])->b

i=18
ss[i]$scaffold -> s
rbind(b, assembly$cov[s, on='scaffold'][bin >= 0.8e6 & bin <= 1e6][order(r)][1, .(scaffold, bin)])->b

i=18
ss[i]$scaffold -> s
rbind(b, assembly$cov[s, on='scaffold'][bin >= 3e6 & bin <= 3.7e6][order(r)][1, .(scaffold, bin)])->b

i=19
ss[i]$scaffold -> s
rbind(b, assembly$cov[s, on='scaffold'][bin >= 0 & bin <= 0.2e6][order(r)][1, .(scaffold, bin)])->b

i=19
ss[i]$scaffold -> s
rbind(b, assembly$cov[s, on='scaffold'][bin >= 0.4e6 & bin <= 0.6e6][order(r)][1, .(scaffold, bin)])->b

#rename column and plot again to double-check
setnames(b, "bin", "br")
plot_chimeras(assembly=assembly, scaffolds=b, br=b, species="lolium", refname="Rabiosa",  mbscale=1,
              file="kyuss_assembly_chimeras_final.pdf", cores=30)

#break the scaffolds
break_scaffolds(b, assembly, prefix="contig_corrected_v1_", slop=1e4, cores=30, species="lolium") -> assembly_v1
saveRDS(assembly_v1, file="kyuss_assembly_v1.Rds")

#load position of restriction fragments
fed <- "kyuss/kyuss.nextpolish.asm.rename_DpnII_fragments_30bp.bed"
read_fragdata(info=assembly_v1$info, file=fed)->frag_data
# # # #
#consider contigs >= 500 kb
frag_data$info[!is.na(hic_chr) & length >= 5e5, .(scaffold, nfrag, length, chr=hic_chr, cM=popseq_cM)]->hic_info
# # # #

hic_map(info=hic_info, assembly=assembly_v1, frags=frag_data$bed, species="lolium", ncores=30,
        min_nfrag_scaffold=30, max_cM_dist = 1000,
        binsize=2e5, min_nfrag_bin=10, gap_size=100)->hic_map_v1

saveRDS(hic_map_v1, file="kyuss_hic_map_v1.Rds")

snuc <- 'kyuss/kyuss.nextpolish.asm.rename_DpnII_fragments_30bp_split.nuc.txt'

hic_plots(rds="kyuss_hic_map_v1.Rds", assembly=assembly_v1, cores=30, species="lolium", nuc=snuc) -> hic_map_v2

#### after checking the Hi-C contact using the map_inspector.Rmd Shinny app, break more chimeric contigs
# load assembly object
readRDS("kyuss_assembly_v1.Rds") -> assembly_v1

data.table(scaffold=c('contig_corrected_v1_15', 
                      'contig_corrected_v1_71', 'contig_corrected_v1_77',
                      'contig_corrected_v1_20',
                      'contig_corrected_v1_42', 
                      'contig_corrected_v1_81',
                      'contig_corrected_v1_83', 'contig_corrected_v1_16',
                      'contig_corrected_v1_33'
)) -> ss
plot_chimeras(assembly=assembly_v1, scaffolds=ss,  species="lolium", refname="Rabiosa", autobreaks=F, mbscale=1,
              file="kyuss_assembly_v1_chimeras.pdf", cores=30)

i=1
ss[i]$scaffold -> s
assembly_v1$cov[s, on='scaffold'][bin >= 12.5e6 & bin <= 13e6][order(r)][1, .(scaffold, bin)] -> b

i=2
ss[i]$scaffold -> s
rbind(b, assembly_v1$cov[s, on='scaffold'][bin >= 70e6 & bin <= 71.3e6][order(r)][1, .(scaffold, bin)])->b

i=3
ss[i]$scaffold -> s
rbind(b, assembly_v1$cov[s, on='scaffold'][bin >= 48.7e6 & bin <= 49e6][order(r)][1, .(scaffold, bin)])->b

i=4
ss[i]$scaffold -> s
rbind(b, assembly_v1$cov[s, on='scaffold'][bin >= 120e6 & bin <= 120.3e6][order(r)][1, .(scaffold, bin)])->b

i=5
ss[i]$scaffold -> s
rbind(b, assembly_v1$cov[s, on='scaffold'][bin >= 57.3e6 & bin <= 57.4e6][order(r)][1, .(scaffold, bin)])->b

i=6
ss[i]$scaffold -> s
rbind(b, assembly_v1$cov[s, on='scaffold'][bin >= 1e6 & bin <= 1.3e6][order(r)][1, .(scaffold, bin)])->b

i=6
ss[i]$scaffold -> s
rbind(b, assembly_v1$cov[s, on='scaffold'][bin >= 2.3e6 & bin <= 2.4e6][order(r)][1, .(scaffold, bin)])->b

i=7
ss[i]$scaffold -> s
rbind(b, assembly_v1$cov[s, on='scaffold'][bin >= 1e6 & bin <= 3e6][order(r)][1, .(scaffold, bin)])->b

i=8
ss[i]$scaffold -> s
rbind(b, assembly_v1$cov[s, on='scaffold'][bin >= 15e6 & bin <= 15.8e6][order(r)][1, .(scaffold, bin)])->b

i=9
ss[i]$scaffold -> s
rbind(b, assembly_v1$cov[s, on='scaffold'][bin >= 0.5e6 & bin <= 1e6][order(r)][1, .(scaffold, bin)])->b

###
setnames(b, "bin", "br")
plot_chimeras(assembly=assembly_v1, scaffolds=b, br=b, species="lolium", refname="Rabiosa",  mbscale=1,
              file="kyuss_assembly_v1_chimeras_final.pdf", cores=30)

break_scaffolds(b, assembly_v1, prefix="contig_corrected_v2_", slop=1e4, cores=30, species="lolium") -> assembly_v2
saveRDS(assembly_v2, file="kyuss_assembly_v2.Rds")

fbed <- 'kyuss/kyuss.nextpolish.asm.rename_DpnII_fragments_30bp.bed'
read_fragdata(info=assembly_v2$info, file=fbed)->frag_data

#consider contigs >= 500 kb
frag_data$info[!is.na(hic_chr) & length >= 5e5, .(scaffold, nfrag, length, chr=hic_chr, cM=popseq_cM)]->hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="lolium", ncores=30,
        min_nfrag_scaffold=30, max_cM_dist = 1000,
        binsize=2e5, min_nfrag_bin=10, gap_size=100)->hic_map_v2

saveRDS(hic_map_v2, file="kyuss_hic_map_v2.Rds")

snuc <- 'kyuss/kyuss.nextpolish.asm.rename_DpnII_fragments_30bp_split.nuc.txt'
hic_plots(rds="kyuss_hic_map_v2.Rds", assembly=assembly_v2, cores=30, species="lolium",
          nuc=snuc) -> hic_map_v2

#### another round of checking chimera with the shinny app
# load assembly object
readRDS("kyuss_assembly_v2.Rds") -> assembly_v2

data.table(scaffold=c('contig_corrected_v2_116',
                      'contig_corrected_v2_120',
                      'contig_corrected_v2_19', 
                      'contig_corrected_v2_117'
                       )) -> ss
plot_chimeras(assembly=assembly_v2, scaffolds=ss,  species="lolium", refname="Rabiosa", autobreaks=F, mbscale=1,
              file="kyuss_assembly_v2_chimeras.pdf", cores=30)

i=1
ss[i]$scaffold -> s
assembly_v2$cov[s, on='scaffold'][bin >= 0.5e6 & bin <= 1.5e6][order(r)][1, .(scaffold, bin)] -> b

i=2
ss[i]$scaffold -> s
rbind(b, assembly_v2$cov[s, on='scaffold'][bin >= 119e6 & bin <= 120e6][order(r)][1, .(scaffold, bin)])->b

i=3
ss[i]$scaffold -> s
rbind(b, assembly_v2$cov[s, on='scaffold'][bin >= 103e6 & bin <= 103.3e6][order(r)][1, .(scaffold, bin)])->b

i=4
ss[i]$scaffold -> s
rbind(b, assembly_v2$cov[s, on='scaffold'][bin >= 15e6 & bin <= 15.1e6][order(r)][1, .(scaffold, bin)])->b

setnames(b, "bin", "br")
plot_chimeras(assembly=assembly_v2, scaffolds=b, br=b, species="lolium", refname="Rabiosa",  mbscale=1,
              file="kyuss_assembly_v2_chimeras_final.pdf", cores=30)

break_scaffolds(b, assembly_v2, prefix="contig_corrected_v3_", slop=1e4, cores=30, species="lolium") -> assembly_v3
saveRDS(assembly_v3, file="kyuss_assembly_v3.Rds")

fbed <- 'kyuss/kyuss.nextpolish.asm.rename_DpnII_fragments_30bp.bed'
read_fragdata(info=assembly_v3$info, file=fbed)->frag_data

#consider contigs >= 500 kb
frag_data$info[!is.na(hic_chr) & length >= 5e5, .(scaffold, nfrag, length, chr=hic_chr, cM=popseq_cM)]->hic_info

hic_map(info=hic_info, assembly=assembly_v3, frags=frag_data$bed, species="lolium", ncores=30,
        min_nfrag_scaffold=30, max_cM_dist = 1000,
        binsize=2e5, min_nfrag_bin=20, gap_size=100)->hic_map_v3

saveRDS(hic_map_v3, file="kyuss_hic_map_v3.Rds")

snuc <- 'kyuss/kyuss.nextpolish.asm.rename_DpnII_fragments_30bp_split.nuc.txt'
hic_plots(rds="kyuss_hic_map_v3.Rds", assembly=assembly_v3, cores=30, species="lolium",
          nuc=snuc) -> hic_map_v3

# another round of checking chimera
# load assembly object
readRDS("kyuss_assembly_v3.Rds") -> assembly_v3

data.table(scaffold=c('contig_corrected_v3_38', 'contig_corrected_v3_67',
                      'contig_corrected_v3_153', 
                      'contig_corrected_v3_138',
                      'contig_corrected_v3_147'
                       )) -> ss

plot_chimeras(assembly=assembly_v3, scaffolds=ss,  species="lolium", refname="Rabiosa", autobreaks=F, mbscale=1,
              file="kyuss_assembly_v3_chimeras.pdf", cores=30)

i=1
ss[i]$scaffold -> s
assembly_v3$cov[s, on='scaffold'][bin >= 80.3e6 & bin <= 80.9e6][order(r)][1, .(scaffold, bin)] -> b

i=2
ss[i]$scaffold -> s
rbind(b, assembly_v3$cov[s, on='scaffold'][bin >= 1.6e6 & bin <= 1.9e6][order(r)][1, .(scaffold, bin)])->b

i=3
ss[i]$scaffold -> s
rbind(b, assembly_v3$cov[s, on='scaffold'][bin >= 102e6 & bin <= 103e6][order(r)][1, .(scaffold, bin)])->b

i=4
ss[i]$scaffold -> s
rbind(b, assembly_v3$cov[s, on='scaffold'][bin >= 0.5e6 & bin <= 1e6][order(r)][1, .(scaffold, bin)])->b

i=5
ss[i]$scaffold -> s
rbind(b, assembly_v3$cov[s, on='scaffold'][bin >= 14.5e6 & bin <= 14.7e6][order(r)][1, .(scaffold, bin)])->b

setnames(b, "bin", "br")
plot_chimeras(assembly=assembly_v3, scaffolds=b, br=b, species="lolium", refname="Rabiosa",  mbscale=1,
              file="kyuss_assembly_v3_chimeras_final.pdf", cores=30)

break_scaffolds(b, assembly_v3, prefix="contig_corrected_v4_", slop=1e4, cores=30, species="lolium") -> assembly_v4
saveRDS(assembly_v4, file="kyuss_assembly_v4.Rds")

fbed <- 'kyuss/kyuss.nextpolish.asm.rename_DpnII_fragments_30bp.bed'
read_fragdata(info=assembly_v4$info, file=fbed)->frag_data

#consider contigs >= 500 kb
frag_data$info[!is.na(hic_chr) & length >= 5e5, .(scaffold, nfrag, length, chr=hic_chr, cM=popseq_cM)]->hic_info

hic_map(info=hic_info, assembly=assembly_v4, frags=frag_data$bed, species="lolium", ncores=21,
        min_nfrag_scaffold=30, max_cM_dist = 1000,
        binsize=2e5, min_nfrag_bin=20, gap_size=100)->hic_map_v4

saveRDS(hic_map_v4, file="kyuss_hic_map_v4.Rds")

snuc <- 'kyuss/kyuss.nextpolish.asm.rename_DpnII_fragments_30bp_split.nuc.txt'
hic_plots(rds="kyuss_hic_map_v4.Rds", assembly=assembly_v4, cores=30, species="lolium",
          nuc=snuc) -> hic_map_v4

# another round of checking chimera
# load assembly object
# readRDS("kyuss_assembly_v4.Rds") -> assembly_v4

data.table(scaffold=c('contig_corrected_v4_168', 'contig_corrected_v4_170',
                      'contig_corrected_v4_159'
                      )) -> ss

plot_chimeras(assembly=assembly_v4, scaffolds=ss,  species="lolium", refname="Rabiosa", autobreaks=F, mbscale=1,
              file="kyuss_assembly_v4_chimeras.pdf", cores=30)

i=1
ss[i]$scaffold -> s
assembly_v4$cov[s, on='scaffold'][bin >= 0.5e6 & bin <= 1e6][order(r)][1, .(scaffold, bin)] -> b

i=2
ss[i]$scaffold -> s
rbind(b, assembly_v4$cov[s, on='scaffold'][bin >= 1.7e6 & bin <= 1.9e6][order(r)][1, .(scaffold, bin)])->b

i=2
ss[i]$scaffold -> s
rbind(b, assembly_v4$cov[s, on='scaffold'][bin >= 4.5e6 & bin <= 5e6][order(r)][1, .(scaffold, bin)])->b

i=2
ss[i]$scaffold -> s
rbind(b, assembly_v4$cov[s, on='scaffold'][bin >= 7.3e6 & bin <= 7.6e6][order(r)][1, .(scaffold, bin)])->b

i=2
ss[i]$scaffold -> s
rbind(b, assembly_v4$cov[s, on='scaffold'][bin >= 7.8e6 & bin <= 8e6][order(r)][1, .(scaffold, bin)])->b

i=3
ss[i]$scaffold -> s
rbind(b, assembly_v4$cov[s, on='scaffold'][bin >= 13e6 & bin <= 14e6][order(r)][1, .(scaffold, bin)])->b

setnames(b, "bin", "br")
plot_chimeras(assembly=assembly_v4, scaffolds=b, br=b, species="lolium", refname="Rabiosa",  mbscale=1,
              file="kyuss_assembly_v4_chimeras_final.pdf", cores=30)

break_scaffolds(b, assembly_v4, prefix="contig_corrected_v5_", slop=1e4, cores=30, species="lolium") -> assembly_v5
saveRDS(assembly_v5, file="kyuss_assembly_v5.Rds")

fbed <- 'kyuss/kyuss.nextpolish.asm.rename_DpnII_fragments_30bp.bed'
read_fragdata(info=assembly_v5$info, file=fbed)->frag_data

#consider contigs >= 500 kb
frag_data$info[!is.na(hic_chr) & length >= 5e5, .(scaffold, nfrag, length, chr=hic_chr, cM=popseq_cM)]->hic_info

hic_map(info=hic_info, assembly=assembly_v5, frags=frag_data$bed, species="lolium", ncores=30,
        min_nfrag_scaffold=30, max_cM_dist = 1000,
        binsize=2e5, min_nfrag_bin=20, gap_size=100)->hic_map_v5

saveRDS(hic_map_v5, file="kyuss_hic_map_v5.Rds")

snuc <- 'kyuss/kyuss.nextpolish.asm.rename_DpnII_fragments_30bp_split.nuc.txt'
hic_plots(rds="kyuss_hic_map_v5.Rds", assembly=assembly_v5, cores=30, species="lolium",
          nuc=snuc) -> hic_map_v5

# output excel file to manually curate contigs
# load assembly object
# assembly_v5 <- readRDS("kyuss_assembly_v5.Rds")

#Export Hi-C map to excel
write_hic_map(rds="kyuss_hic_map_v5.Rds", file="kyuss_hic_map_v5.xlsx", species="lolium")

# after editing the excel, now reconstruct the map
# load assembly object
# assembly_v5 <- readRDS("kyuss_assembly_v5.Rds")

#make edits in Excel and import
read_hic_map(rds="kyuss_hic_map_v5.Rds", file="kyuss_hic_map_v5_edited.xlsx") -> nmap
diff_hic_map(rds="kyuss_hic_map_v5.Rds", nmap, species="lolium") 
hic_map(species="lolium", agp_only=T, map=nmap)->hic_map_v6

saveRDS(hic_map_v6, file="kyuss_hic_map_v6.Rds")

#final Hi-C plots
snuc <- 'kyuss/kyuss.nextpolish.asm.rename_DpnII_fragments_30bp_split.nuc.txt'
hic_plots(rds="kyuss_hic_map_v6.Rds", assembly=assembly_v5, cores=30, species="lolium", nuc=snuc) -> hic_map_v6

# another round output excel file to manually curate contigs
# load assembly object
# assembly_v5 <- readRDS("kyuss_assembly_v5.Rds")

#Export Hi-C map to excel
write_hic_map(rds="kyuss_hic_map_v6.Rds", file="kyuss_hic_map_v6.xlsx", species="lolium")

# after editing the excel, now reconstruct the map

#load TRITEX R code
source('/home/yutachen/public/Yutangchen/Kyuss_data/tritex/bitbucket/pseudomolecule_construction.R')

# load assembly object
assembly_v5 <- readRDS("kyuss_assembly_v5.Rds")

#make edits in Excel and import
read_hic_map(rds="kyuss_hic_map_v6.Rds", file="kyuss_hic_map_v6_edited.xlsx") -> nmap
diff_hic_map(rds="kyuss_hic_map_v6.Rds", nmap, species="lolium") 
hic_map(species="lolium", agp_only=T, map=nmap)->hic_map_v7

saveRDS(hic_map_v7, file="kyuss_hic_map_v7.Rds")

#final Hi-C plots
snuc <- 'kyuss/kyuss.nextpolish.asm.rename_DpnII_fragments_30bp_split.nuc.txt'
hic_plots(rds="kyuss_hic_map_v7.Rds", assembly=assembly_v5, cores=30, species="lolium", nuc=snuc) -> hic_map_v7

############################################################################################
### almost there, now construct pseudomolecules
#load TRITEX R code
source('/home/yutachen/public/Yutangchen/Kyuss_data/tritex/bitbucket/pseudomolecule_construction.R')

library(xlsx)

# load assembly object
assembly_v5 <- readRDS("kyuss_assembly_v5.Rds")

# load hic map
hic_map_v7 <- readRDS("kyuss_hic_map_v7.Rds")

#pseudmolecule stats
hic_map_v7$agp[agp_chr != "chrUn" & gap == F][, .("Ncontig" = .N, "N50 (Mb)"=round(n50(scaffold_length/1e6), 1),
                                                   "max_length (Mb)"=round(max(scaffold_length/1e6),1),
                                                   "min_length (kb)"=round(min(scaffold_length/1e3),1)), key=agp_chr] -> res
hic_map_v7$chrlen[, .(agp_chr, "length (Mb)"=round(length/1e6, 1))][res, on="agp_chr"] -> res
setnames(res, "agp_chr", "chr")

hic_map_v7$agp[gap == F & agp_chr != "chrUn"][, .("Ncontig" = .N, "N50 (Mb)"=round(n50(scaffold_length/1e6), 1),
                                                   "max_length (Mb)"=round(max(scaffold_length/1e6),1),
                                                   "min_length (kb)"=round(min(scaffold_length/1e3),1))][, agp_chr := "1-7"] -> res2
hic_map_v7$chrlen[, .(agp_chr="1-7", "length (Mb)"=round(sum(length)/1e6, 1))][res2, on="agp_chr"] -> res2
setnames(res2, "agp_chr", "chr")

hic_map_v7$agp[gap == F & agp_chr == "chrUn"][, .(chr="un", "Ncontig" = .N, "N50 (Mb)"=round(n50(scaffold_length/1e6), 1),
                                                   "max_length (Mb)"=round(max(scaffold_length/1e6),1),
                                                   "min_length (kb)"=round(min(scaffold_length/1e3),1),
                                                   "length (Mb)"=sum(scaffold_length/1e6))] -> res3

rbind(res, res2, res3) -> res

write.xlsx(res, file="230215_kyuss_hic_map_v7_pseudomolecule_stats.xlsx")

#compile pseudomolecules
fasta <- "kyuss/kyuss.nextpolish.asm.rename.fasta"
sink("230215_kyuss_nextdenovo_pseudomolecules_v1.log")
compile_psmol(fasta=fasta, output="230215_kyuss_new_pseudomolecules_v1",
              hic_map=hic_map_v7, assembly=assembly_v5, cores=31)
sink()

print("pseudomolecules compiling finished")

################################################################################################################################
# make fina hic contact matrix

# load assembly object
assembly_v5 <- readRDS("kyuss_assembly_v5.Rds")

#final Hi-C plots
snuc <- 'kyuss/kyuss.nextpolish.asm.rename_DpnII_fragments_30bp_split.nuc.txt'
hic_plots(rds="kyuss_hic_map_v7.Rds", assembly=assembly_v5, cores=30, species="lolium", nuc=snuc) -> hic_map_v7

#plot intra-chromosomal map in 7x3 grid
hmap <- hic_map_v7
ncol=100
colorRampPalette(c("white", "red"))(ncol)->whitered
links <- hmap$hic_1Mb$norm
chrs <- unique(links$chr1)
links[chr1 %in% chrs, .(chr=chr1, bin1, bin2, l=log10(nlinks_norm))]->z
chrNames(agp=T, species="lolium")[z, on="chr"]->z
binsize <- min(links[dist > 0]$dist)
z[, col := whitered[cut(l, ncol, labels=F)]]

png(file="230215_kyuss_nextdenovo_Hi-C_intrachromosomal_contact_matrices.png", res=200, height=2000, width=4000)
par(mfrow=c(2,4), cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
par(mar=c(5,5,3,3))
lapply(chrs, function(i){
  z[chr == i, plot(0, type='n', bty='l', xlim=range(bin1/1e6), ylim=range(bin2/1e6), xlab="position (Mb)", ylab="position (Mb)", col=0, 
                   main=paste("Chr", chrNames(species="lolium")[chr == i, alphachr], sep = ""))]
  # hmap$agp[gap == T & agp_chr == chrNames(agp=T, species="lolium")[chr == i, agp_chr], abline(lwd=1, col='gray', v=(agp_start+agp_end)/2e6)]
  # hmap$agp[gap == T & agp_chr == chrNames(agp=T, species="lolium")[chr == i, agp_chr], abline(lwd=1, col='gray', h=(agp_start+agp_end)/2e6)]
  z[chr == i, rect((bin1-binsize)/1e6, (bin2-binsize)/1e6, bin1/1e6, bin2/1e6, col=col, border=NA)]
})
dev.off()

##################################################################################################################################
### dot plot between kyuss pseudo molecule and rabiosa guide map
#load TRITEX R code
source('/home/yutachen/public/Yutangchen/HiCproject/bitbucket/pseudomolecule_construction.R')

library(xlsx)

# load hic map
# hic_map_v7 <- readRDS("kyuss_hic_map_v7.Rds")

# read GMI423 to psmol PAF and make alignment plots
pf <- 'kyuss/rabiosa_to_kyuss_pmolecule.paf.gz'
read_paf(pf)->paf
fread('pmolecule/rabiosa_unphased_chr_singlecopy_100bp.bed', head=F)->b
b[, V2 := V2 + 1]
setnames(b, c("agp_chr_qry", "agp_start_qry", "agp_end_qry", "query"))
b[paf, on="query"] -> p
#remove short and non-unique alignment
p[alnlen > 50 & mapq >= 30] -> p

#calculate offsets, 10 Mb between chromosomes
hic_map_v7$chrlen -> yy
yy[, plot_offset := cumsum(c(0, length[-.N] + 1e7))]

#read GMI423 Hi-C map and calculate offsets
# readRDS('MPB_FM13_hic_map_v4.Rds')$chrlen->xx
# xx[, plot_offset := cumsum(c(0, length[-.N] + 1e7))]

fai <- '/home/yutachen/public/Yutangchen/Sikem_hifi/tritex_utg/rabiosa_chr/rabiosa_unphased_chr.fa.fai'
fread(head=F, sel=1:2, col.names=c("agp_chr", "length"), fai)->xx
rbind(xx[1:7][order(substr(agp_chr, 5,5))], xx[22])->xx
xx[, plot_offset := cumsum(c(0, length[-.N] + 1e7))]


#merge AGPs with PAF and calculate plot positions
p[reference != "chrUn" & agp_chr_qry != "chrUn"] -> zz
xx[, .(agp_chr_qry=agp_chr, offq=plot_offset)][zz, on="agp_chr_qry"] -> zz
yy[, .(reference=agp_chr, offr=plot_offset)][zz, on="reference"] -> zz
zz[, ppos_r := offr + (reference_start + reference_end)/2]
zz[, ppos_q := offq + (agp_start_qry + agp_end_qry)/2]

png("230215_rabiosa_to_kyuss_psmol.png", res=300, height=11*300, width=11*300)
par(mar=c(5,5,2,1))
zz[, plot(xpd=NA, ppos_r, xaxs='i', yaxs='i',
          ylim=c(-2e7, xx[8, plot_offset + 2e7]),
          xlim=c(-2e7, yy[8, plot_offset + 2e7]),
          type='n', xlab="", bty='n', ylab="", xaxt='n', yaxt='n', ppos_q)]
title(ylab="position in Rabiosa pseudomolecules")
title(xlab="position in Kyuss pseudomolecules") #change names for other genotypes
yy[1:7, axis(1, plot_offset + length/2, srt=90, sub("chr", "", agp_chr), tick=F)]
yy[, abline(v=plot_offset - 5e6, col='blue')]
xx[1:7,  axis(2, plot_offset + length/2, las=1, sub("chr", "", agp_chr), tick=F)]
xx[, abline(h=plot_offset - 5e6, col='blue')]
zz[, points(pch=".", col="#00000011", ppos_r, cex=1, ppos_q)]
dev.off()

#intrachromosomal plot
p[reference == agp_chr_qry] -> qq
png("230215_rabiosa_to_kyuss_psmol_intra.png", height=3*300, width=21*300, res=300)
par(mar=c(5,5,3,3))
par(mfcol=c(1,7))
lapply(sort(unique(xx[1:7]$agp_chr)), function(i){
  qq[i, on="reference"][, plot(pch=".", cex=1,  bty='l', las=1, xlim=c(0,  hic_map_v7$chrlen[agp_chr == i]$length/1e6),
                               col="#00000033", xlab="Rabiosa_new (Mb)", ylab="Rabiosa_old (Mb)",
                               reference_start/1e6, agp_start_qry/1e6, main=i)]
  yy[agp_chr == i, abline(v=c(0, length/1e6), col='gray')]
  xx[agp_chr == i, abline(h=c(0, length/1e6), col='gray')]
})
dev.off()


