#!/bin/bash

# 02.04.2023

# using polca to find homozygous SNPs in kyuss assembly

# copy assembly here
# cp ~/public/Yutangchen/Kyuss_data/kyuss_published_v1/Kyuss_1697_assembly.fa ./

# run pocal
polca.sh -a Kyuss_1697_assembly.fa -r '/scratch/yutang/kyuss_nextpolish/20201113.B-Kyuss_39_R1.fastq.gz /scratch/yutang/kyuss_nextpolish/20201113.B-Kyuss_39_R2.fastq.gz' -t 48 -m 3G
