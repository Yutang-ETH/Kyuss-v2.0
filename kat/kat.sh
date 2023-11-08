#!/bin/bash

# 06.02.2023

# KAT kmer analysis

# using 23 mer
myread="/scratch/yutang/kyuss_nextpolish"
wd="/scratch/yutang/kyuss_nextpolish/KAT"
assembly="/scratch/yutang/kyuss_nextpolish/kyuss_polish/genome.nextpolish.fasta"
 
kat comp -t 30 -m 23 -h -H 50000000000 -o ${wd}/ "${myread}/20201113.B-Kyuss_39_R?.fastq.gz" ${assembly} 
kat plot spectra-cn -x 200 -o ${wd}/kyuss_nextpolish_kat.png ${wd}/-main.mx
