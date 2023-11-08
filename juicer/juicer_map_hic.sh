#!/bin/bash

# map hic using juicer.sh

# actually should use --assembly with -e because -e doesn't creat noredup.txt

wd=/scratch/yutang/manual_curation/kyuss
ref=kyuss.nextdenovo.fasta
site=kyuss_DpnII.txt
ID="kyuss"
enzyme="DpnII"

${wd}/scripts/juicer.sh -e --assembly -t 60 -T 60 -g ${ID} -d ${wd}/ -p ${wd}/references/chrom.sizes -s ${enzyme} -y ${wd}/restriction_sites/${site} -z ${wd}/references/${ref} -D ${wd}/ > myjuicer.stdout 2> myjuicer.stderr
