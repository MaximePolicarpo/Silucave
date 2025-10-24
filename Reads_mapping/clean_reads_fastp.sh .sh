#!/bin/bash


#SBATCH --job-name=clean_reads  # Job name


module purge ; module load fastp

species=$1

if ls -l | grep -q "$species.R1.fastq.gz" ; then
	fastp -i $species.R1.fastq.gz -I $species.R2.fastq.gz -o $species.R1.cleaned.fastq.gz -O $species.R2.cleaned.fastq.gz -h $species.html -j $species.json
	rm $species.R1.fastq.gz ; rm $species.R2.fastq.gz
else
	fastp -i $species.fastq.gz -o $species.cleaned.fastq.gz
	rm $species.fastq.gz
fi
