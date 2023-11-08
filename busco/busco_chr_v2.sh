#!/bin/bash

# 05.02.2023

# busco of kolumbus polished assembly

assembly="/scratch/yutang/kyuss_polca/kyuss_v2/kyuss.nextdenovo.juicer.chr.fasta"
database="/scratch/yutang/busco_database/embryophyta_odb10"

mkdir v2_genome
busco -i ${assembly} \
      -l ${database} \
      -m genome \
      -o v2_genome \
      -f \
      -c 24 \
      --offline
