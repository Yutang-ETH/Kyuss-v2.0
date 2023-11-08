#!/bin/bash

# 09.02.2023

# busco of kolumbus polished assembly

assembly="/scratch/yutang/kyuss_nextpolish/kyuss_polish/genome.nextpolish.fasta"
database="/scratch/yutang/busco_database/embryophyta_odb10"

mkdir unphased_genome
busco -i ${assembly} \
      -l ${database} \
      -m genome \
      -o unphased_genome \
      -f \
      -c 20 \
      --offline


