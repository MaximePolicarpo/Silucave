#!/bin/bash

#SBATCH --job-name=trimm  # Job name

module load trimal

i=$1

trimal -in $i -out $i.trimal -automated1