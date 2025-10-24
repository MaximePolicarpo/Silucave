#!/bin/bash

species_name=$1 

grep "$species_name" SRA_data_filt.txt | cut -f2- -d "," | tr ',' '\n' > $species_name.SRA.id

for ID in `cat $species_name.SRA.id` ; do
	fasterq-dump --split-files --skip-technical $ID
	gzip $ID.fastq
done

for ID in `cat $species_name.SRA.id` ; do
	cat $ID.fastq.gz >> $species_name.fastq.gz
done


for ID in `cat $species_name.SRA.id` ; do
	rm $ID.fastq.gz
done