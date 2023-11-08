#!/bin/bash

# 01.04.2023

# mask Kyuss v2.0 using repeatmask with the repeat database Dario produced based on Rabiosa NRGene's assembly

RepeatMasker -pa 60 -lib Rabiosa_repeats_centroids_codes_180830.fa -qq -norna -no_is -gff -cutoff 250 -gccalc -engine ncbi -dir . kyuss.nextdenovo.juicer.fasta
