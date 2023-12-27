#!/bin/bash

#run in silico digest, digest_emboss.zsh is on the TRITEX bitbucket

ref='kyuss/kyuss.nextpolish.asm.rename.fasta'

/home/yutachen/public/Yutangchen/Kyuss_data/tritex/bitbucket/digest_emboss.zsh --ref $ref --enzyme 'DpnII' --sitelen 4 --minlen 30 
