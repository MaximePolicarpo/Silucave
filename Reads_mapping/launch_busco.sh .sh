#!/bin/bash


#SBATCH --job-name=BUSCO   # Job name


eval "$(conda shell.bash hook)"
conda activate BUSCO_5_6_1

species=$1
assembly_name=$2

busco -i $species/$assembly_name -m genome -c 40 -l actinopterygii_odb10 -o BUSCO_$species --offline 

