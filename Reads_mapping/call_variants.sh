#!/bin/bash


#SBATCH --job-name=read_maps  # Job name

ref_genome=$1
sample_name=$2


module purge ; module load BCFtools

#Perform variant calling with BCFtools
bcftools mpileup -Ou --threads 20 -f $ref_genome $sample_name.sorted.bam | bcftools call --threads 20 -mv -Oz -o $sample_name.vcf.gz 

#Normalize indels
bcftools norm -f $ref_genome $sample_name.vcf.gz -Ob -o $sample_name.norm.bcf --threads 20 ; rm $sample_name.vcf.gz

#Filter adjacent indels withing 5bp
bcftools filter --threads 8 --IndelGap 5 $sample_name.norm.bcf -Ob -o $sample_name.norm.flt-indels.bcf ; rm $sample_name.norm.bcf

#Index file
bcftools index $sample_name.norm.flt-indels.bcf --threads 20
