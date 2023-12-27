#!/usr/bin/env zsh

scaf="kyuss/kyuss.nextpolish.asm.fasta"
ref='kyuss/kyuss.nextpolish.asm.rename.fasta'
echo $ref:r.n50

awk '/^>/ {print "\n>contig_"++n; next} {print $0}' $scaf | awk NF > $ref && {samtools faidx $ref && ./bitbucket/n50 $ref.fai > $ref:r.n50}
