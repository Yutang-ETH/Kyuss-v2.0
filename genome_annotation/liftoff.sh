#!/bin/bash

# lift genes in kyuss v1 to kyuss v2 using liftoff

# 01.04.2023

# installed liftoff from source using setup.py
# conda and pip installation failed

# reference is kyuss v1
ref=Kyuss_1697_assembly.fa

# gff is the annotation file of kyuss v1
gff=Kyuss_1697_KYUS.gff

# target is kyuss v2
target=kyuss.nextdenovo.juicer.fasta

# run liftoff
liftoff -g ${gff} -o liftoff_annotation.gff ${target} ${ref}
