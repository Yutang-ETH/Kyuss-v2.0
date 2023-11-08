#!/bin/bash

# using liftofftools to compare the genome annotation between Kyuss v1.0 and Kyuss v2.0

# 01.04.2023

# reference is Kyuss v1.0
ref=Kyuss_1697_assembly.fa

# reference annotation
refgff=Kyuss_1697_KYUS.gff

# target is Kyuss v2.0
target=kyuss.nextdenovo.juicer.fasta

# target liftoff annotation
targff=liftoff_annotation_kyuss_v2.gff

# run liftoff tools to analyze synteny between two versions
liftofftools synteny -r ${ref} -t ${target} -rg ${refgff} -tg ${targff} -r-sort rorder.txt -t-sort torder.txt
