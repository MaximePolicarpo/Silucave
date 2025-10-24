#!/bin/bash

LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

eval "$(conda shell.bash hook)"
conda activate miniprot

curr_OGG=$1


rm Coding_sequences_alignments/$curr_OGG.prot.fai
samtools faidx Coding_sequences_alignments/$curr_OGG.prot

longest_gene=`sort -k2 -n Coding_sequences_alignments/$curr_OGG.prot.fai | tail -1 | cut -f1`

samtools faidx Coding_sequences_alignments/$curr_OGG.prot $longest_gene | sed "s/>.*/>$curr_OGG/g"  > Coding_sequences_alignments/$curr_OGG.reprez_seq.prot

rm Coding_sequences_alignments/$curr_OGG.prot.fai