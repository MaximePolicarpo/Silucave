#!/bin/bash

#SBATCH --job-name=test   # Job name


Genome=$1
prot_Seq=$2
species=$3
working_dir=$4

export SLURM_TMPDIR=/tmp/
export TMPDIR=/tmp/


singularity exec -B $HOME:$HOME braker3.sif braker.pl --genome=$Genome --prot_seq=$prot_Seq --species=$species --workingdir=$working_dir --threads 25 --gff3 --AUGUSTUS_CONFIG_PATH=/home/.augustus/
