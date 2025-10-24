#!/bin/bash


#SBATCH --job-name=read_maps  # Job name

ref_genome=$1
sample_name=$2

module purge ; module load bwa-mem2

bwa-mem2 mem -t 40 $ref_genome $sample_name.R1.cleaned.fastq.gz $sample_name.R2.cleaned.fastq.gz > $sample_name.sam

module purge ; module load SAMtools

samtools sort $sample_name.sam -o $sample_name.sorted.bam ; samtools index $sample_name.sorted.bam ; rm $sample_name.sam 

echo "Computing mean depth ..."
samtools depth -a $sample_name.sorted.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' #the -a flag allows to also take into account regions with no reads mapped

