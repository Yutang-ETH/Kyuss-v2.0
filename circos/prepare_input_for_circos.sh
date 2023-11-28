#!/bin/bash

# 01.04.2023

# calculate GC content in every 1 Mb window using seqkit
cat kyuss.nextdenovo.juicer.chr.fasta | seqkit sliding -s 1000000 -W 1000000 -g | seqkit fx2tab -n -g > kyuss_v2_chr.gc &

# make a gene bed from the gff file
grep 'chr' liftoff_annotation_kyuss_v2.gff | grep 'gene' | cut -f1,4,5 > gene.bed

# make a repeat bed from the repeatmask gff file
# grep 'chr' kyuss.nextdenovo.juicer.fasta.out.gff | cut -f1,4,5 > repeat.bed
seqkit locate -P -r -p 'N+' --bed -G ../repeatmask/kyuss.nextdenovo.juicer.fasta.masked | cut -f1-3 > repeat.bed

# make a faidx file for kyuss v2 chromosomes
samtools faidx kyuss.nextdenovo.juicer.chr.fasta

# find gaps in the pseudo-chromosomes
seqkit locate -p N+ -G -P -r -i --bed kyuss_v1_chr.fasta > kyuss_v1_chr_gap.bed
seqkit locate -p N+ -G -P -r -i --bed kyuss_v2_chr.fasta > kyuss_v2_chr_gap.bed

