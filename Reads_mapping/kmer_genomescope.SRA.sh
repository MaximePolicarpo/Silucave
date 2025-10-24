#!/bin/bash


#SBATCH --job-name=Scope   # Job name

eval "$(conda shell.bash hook)"
conda activate GenomeScope

species_ID=$1

if ls -l | grep -q "$species_ID.R1.cleaned.fastq.gz" ; then
	echo "$species_ID.R1.cleaned.fastq.gz" > FILES.$species_ID
	echo "$species_ID.R2.cleaned.fastq.gz" >> FILES.$species_ID
else 
	echo "$species_ID.cleaned.fastq.gz" > FILES.$species_ID
fi


mkdir temp_dir.$species_ID/
kmc -k21 -t10 -m64 -ci1 -cs10000 @FILES.$species_ID $species_ID.reads /scicore/home/salzburg/polica0000/SiluCave/Silucave_data/temp_dir.$species_ID/
kmc_tools transform $species_ID.reads histogram $species_ID.reads.histo -cx10000
genomescope2 -i $species_ID.reads.histo -o $species_ID.GenomeScope_rslt -k 21
