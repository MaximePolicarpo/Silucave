#!/bin/bash


#SBATCH --job-name=Grantham   # Job name


module purge ; module load R/4.4.1-foss-2023b


curr_file=$1
curr_output=$2

Rscript compute_grantham_dist.R $1 $2