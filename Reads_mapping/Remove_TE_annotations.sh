#!/bin/bash


#SBATCH --job-name=Remove_TE   # Job name


#conda = miniprot ; cd-hit ; seqkit

eval "$(conda shell.bash hook)"
conda activate RepeatMasker_Env
module purge ; module load SAMtools

species=$1

rm -r $species.RepeatMasker_dir ; mkdir $species.RepeatMasker_dir
RepeatMasker -pa 40 -s -e ncbi -a -no_is -gff -species Actinopterygii Coding_sequences/$species.cds -dir $species.RepeatMasker_dir >& $species.repeatmasker.out

grep -v "Simple_repeat\|Low_complexity\|Satellite\| tRNA \| rRNA \| snRNA " $species.RepeatMasker_dir/$species.cds.out  | grep -v "begin.*end" | grep -v "position in query" | sed 's/  */	/g' | cut -f6 | sort | uniq | sed '/^[[:space:]]*$/d' > $species.RepeatMasker_dir/Genes_TE.id

grep ">" Coding_sequences/$species.cds | sed 's/>//g' | sort > $species.RepeatMasker_dir/all_id.txt
comm -23 <(sort $species.RepeatMasker_dir/all_id.txt) <(sort $species.RepeatMasker_dir/Genes_TE.id) > $species.RepeatMasker_dir/all_id_woTE.txt

xargs samtools faidx Coding_sequences/$species.cds < $species.RepeatMasker_dir/all_id_woTE.txt > Coding_sequences_woTE/$species.cds
xargs samtools faidx Proteins/$species.prot < $species.RepeatMasker_dir/all_id_woTE.txt > Proteins_woTE/$species.prot
