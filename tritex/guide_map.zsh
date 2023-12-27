#get single-copy sequences by maskin 31 mer occuring at least twice
#https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/miscellaneous/mask_assembly.zsh
# fa='rabiosa_chr/rabiosa_unphased_chr.fa'
# mask='/home/yutachen/public/Yutangchen/Sikem_hifi/tritex_utg/bitbucket/mask_assembly.zsh'
# zsh $mask --mem 300G --fasta $fa --mincount 2 --out "." &

#extract BED and FASTA files for single-copy sequences >= 100 bp
ref='rabiosa_chr/rabiosa_unphased_chr.fa'
bed='rabiosa_unphased_chr_masked_noGaps.bed'
out="${ref:r}_singlecopy_100bp"

awk '$3 - $2 >= 100' $bed | grep -v chrUn | awk '{print $0"\tseq_"NR}' \
  | tee $out.bed | bedtools getfasta -fi $ref -bed /dev/stdin -name -fo $out.fasta

sed -i "s/:.*//g" rabiosa_chr/rabiosa_unphased_chr_singlecopy_100bp.fasta
