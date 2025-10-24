#!/bin/bash

#SBATCH --job-name=OrthoFinder 

eval "$(conda shell.bash hook)" 
conda activate OrthoFinder_env

species_tree=$1
proteome=$2

orthofinder -t 50 -a 15 -s $species_tree -f $proteome