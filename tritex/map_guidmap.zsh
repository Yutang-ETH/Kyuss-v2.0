#!/usr/bin/env zsh 

# align Asativa_GMI423 single-copy regions to contigs

minimap='/home/yutachen/anaconda3/envs/tritex/bin/minimap2' 
ref="kyuss/kyuss.nextpolish.asm.rename.fasta" 
qry='pmolecule/rabiosa_unphased_chr_singlecopy_100bp.fasta' 
size='5G' 
mem='5G' 
threads=20  
{
 $minimap -t $threads -2 -I $size -K $mem -x asm5 $ref $qry | bgzip > kyuss/kyuss_rabiosa.paf.gz
} 2> kyuss/kyuss_rabiosa.paf.err 
