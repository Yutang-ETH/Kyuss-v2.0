#!/bin/bash

# 05.02.2023

# busco of kolumbus polished assembly

assembly="/scratch/yutang/kyuss_polca/kyuss_v1/Kyuss_1697_assembly_chr.fa"
database="/scratch/yutang/busco_database/embryophyta_odb10"

mkdir v1_genome
busco -i ${assembly} \
      -l ${database} \
      -m genome \
      -o v1_genome \
      -f \
      -c 24 \
      --offline
