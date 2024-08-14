
* copy kyuss nextploish asm to Kyuss_data/tritex/kyuss/

* now change the contig name

    `./ChangeScaffoldsName.zsh`

* extract single-copy kmers or regions

    `zsh guide_map.zsh`

* insilico digestion of the assembly

    `./IN_SILICO_digest.sh`

* copy hic file to tritex folder from /scratch/yutang/manual_curation/kyuss/fastq

    `cp Kyuss_hic_R1.fastq.gz ~/public/Yutangchen/Kyuss_data/tritex/kyuss/Hic_reads_R1.fastq.gz &`
    `cp Kyuss_hic_R2.fastq.gz ~/public/Yutangchen/Kyuss_data/tritex/kyuss/Hic_reads_R2.fastq.gz &`

* map hic data

    `HicMapping.sh`

* align to the guide map (reference), rabiosa assembly, https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_030979885.1/

    `map_guidmap.zsh`

* make popseq object, arg[1] single-copy kmer bed, arg[2] the fai file of the guide map (this case Rabiosa), arg[3] the prefix of the output file, supposed to be the name of the guide map 

* Tritex does not provide this MAKE_pseudo_popseq.R script, I wrote it

    `Rscript MAKE_pseudo_popseq.R /home/yutachen/public/Yutangchen/Kyuss_data/tritex/pmolecule/rabiosa_unphased_chr_singlecopy_100bp.bed /home/yutachen/public/Yutangchen/Sikem_hifi/tritex_utg/rabiosa_chr/rabiosa_unphased_chr.fa.fai Rabiosa`

* run the Rscript interactively, I mean copy the block of code to the R interface and run the code

    `curation_1.R`

* for more details of how to run Tritex, please contact Martin Mascher, IPK.
