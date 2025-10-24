#!/bin/bash


#SBATCH --job-name=abSREL   # Job name


eval "$(conda shell.bash hook)"
conda activate HyPhy_Env

OGG_alignment=$1
OGG_tree=$2
curr_OGG=$3

#Assign cavefish genes as test branches

cp $OGG_tree $OGG_tree.marked
for cave_species in `cat cavefish_species.txt` ; do grep ">$cave_species" $OGG_alignment | sed 's/>//g' ; done > Coding_sequences_alignments/$curr_OGG.caveID
for branch in `cat Coding_sequences_alignments/$curr_OGG.caveID` ; do  sed -i "s/$branch/$branch {test}/g" $OGG_tree.marked ; done
rm Coding_sequences_alignments/$curr_OGG.caveID

#launch abSREL
hyphy absrel --alignment $OGG_alignment --tree $OGG_tree.marked --code Universal --branches test ENV=TOLERATE_NUMERICAL_ERRORS=1

#launch RELAX
hyphy relax --alignment $OGG_alignment --tree $OGG_tree.marked --code Universal --test test --reference foreground ENV=TOLERATE_NUMERICAL_ERRORS=1
