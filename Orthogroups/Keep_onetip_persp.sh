#!/bin/bash

LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge ; module load R/4.4.1-foss-2023b


current_tree_rooted=$1 
vector_good_tips=$2
curr_OGG=$3 

#Keep only one tip per spercies
Rscript Keep_onetip_persp.R $current_tree_rooted $vector_good_tips


module purge ; module load SAMtools

#Extract sequences
xargs samtools faidx Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln < Coding_sequences_alignments/$curr_OGG.OnetipPersp.txt > Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.OneTipPerSp
xargs samtools faidx Coding_sequences_alignments/$curr_OGG.cds.trimmed < Coding_sequences_alignments/$curr_OGG.OnetipPersp.txt > Coding_sequences_alignments/$curr_OGG.cds.trimmed.OneTipPerSp

#Rename sequences
for species in `cat species_list.txt` ; do
	sed -i "s/>$species.*/>$species/g" Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.OneTipPerSp
	sed -i "s/>$species.*/>$species/g" Coding_sequences_alignments/$curr_OGG.cds.trimmed.OneTipPerSp
done


#Estimate ML tree with phangorn

module purge ; module load R/4.4.1-foss-2023b

Rscript Estimate_ML_tree_Phangorn.R Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.OneTipPerSp Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.OneTipPerSp.nwk

echo "Nwk file generated"
