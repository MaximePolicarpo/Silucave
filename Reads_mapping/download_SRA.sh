#!/bin/bash


species_name=$1 

grep "$species_name" SRA_data_filt.txt | cut -f2- -d "," | tr ',' '\n' > $species_name.SRA.id

for ID in `cat $species_name.SRA.id` ; do
	fasterq-dump --split-files --skip-technical $ID
	gzip $ID\_1.fastq ; gzip $ID\_2.fastq
done

for ID in `cat $species_name.SRA.id` ; do
	cat $ID\_1.fastq.gz >> $species_name.R1.fastq.gz
	cat $ID\_2.fastq.gz >> $species_name.R2.fastq.gz
done


for ID in `cat $species_name.SRA.id` ; do
	rm $ID\_1.fastq.gz
	rm $ID\_2.fastq.gz
done