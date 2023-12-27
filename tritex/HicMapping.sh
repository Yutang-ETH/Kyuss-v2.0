#!/bin/bash

# mkdir tmp
ref="kyuss/kyuss.nextpolish.asm.rename.fasta"
bed="kyuss/kyuss.nextpolish.asm.rename_DpnII_fragments_30bp.bed"
linker="GATCGATC"
TMPDIR="./tmp"

/home/yutachen/public/Yutangchen/Kyuss_data/tritex/bitbucket/run_hic_mapping.zsh --threads 48 --mem '100G' --ref $ref --linker $linker --bed $bed --tmp $TMPDIR kyuss
