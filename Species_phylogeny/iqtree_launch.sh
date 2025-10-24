#!/bin/bash


#SBATCH --job-name=ML_model_iqtree   # Job name

eval "$(conda shell.bash hook)"

module load IQ-TREE

gene_alignment=$1 

iqtree -s $gene_alignment --seqtype AA -m TEST -mrate G4 -nt 10 -bb 1000