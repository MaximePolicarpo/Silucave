#!/bin/bash


#SBATCH --job-name=Collapse_low_support   # Job name



#conda = miniprot ; cd-hit ; seqkit


LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
module load R/4.4.1-foss-2023b

original_tree=$1
collapsed_tree=$2

Rscript collapse_lowsupport.R $original_tree $collapsed_tree