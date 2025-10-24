#!/bin/bash


#SBATCH --job-name=Scope   # Job name

eval "$(conda shell.bash hook)"
conda activate GenomeScope

species_ID=$1
kmer_size=$2

rm $species_ID.reads.histo
rm $species_ID.reads.kmc_pre
rm $species_ID.reads.kmc_suf
rm -rf $species_ID.GenomeScope_rslt/
rm -rf temp_dir.$species_ID/

echo "$species_ID.R1.cleaned.fastq.gz" > FILES.$species_ID
echo "$species_ID.R2.cleaned.fastq.gz" >> FILES.$species_ID

mkdir temp_dir.$species_ID/
kmc -k$kmer_size-t10 -m64 -ci1 -cs10000 @FILES.$species_ID $species_ID.reads temp_dir.$species_ID/
kmc_tools transform $species_ID.reads histogram $species_ID.reads.histo -cx10000
genomescope2 -i $species_ID.reads.histo -o $species_ID.GenomeScope_rslt -k $kmer_size
