#!/bin/bash

# 01.04.2023

# using polca to find homozygous SNPs in kyuss assembly

# copy assembly here
# cp ~/public/Yutangchen/Kyuss_data/kyuss_nextpolish/hicexplorer/assembly/kyuss.nextdenovo.juicer.fasta ./

# run pocal
polca.sh -a kyuss.nextdenovo.juicer.fasta -r '/scratch/yutang/kyuss_nextpolish/20201113.B-Kyuss_39_R1.fastq.gz /scratch/yutang/kyuss_nextpolish/20201113.B-Kyuss_39_R2.fastq.gz' -t 48 -m 3G
