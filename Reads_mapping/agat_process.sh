#!/bin/bash

#SBATCH --job-name=agat_parsing   # Job name

eval "$(conda shell.bash hook)"
conda activate myagat_env 

module purge ; module load EMBOSS 


species=$1

GFF_file=`ls -l $species | grep ".gff" | sed 's/.* //g'`
genome_file=`ls -l $species | grep ".fna$" | sed 's/.* //g'`

agat_sp_keep_longest_isoform.pl --gff $species/$GFF_file -o $species/$GFF_file.LI #keep longest isoform
agat_sp_extract_sequences.pl -gff $species/$GFF_file.LI --fasta $species/$genome_file -t cds -o $species/$species.cds #extract the CDS
sed 's/ .*//g' $species/$species.cds | sed 's/>.*|/>/g' > $species/$species.renamed.cds
sed -i 's/:/_/g' $species/$species.renamed.cds #remove ":" special character
sed -i "s/>/>$species---/g" $species/$species.renamed.cds #Add species name to headers
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' $species/$species.renamed.cds > $species/$species.temp.cds ; mv $species/$species.temp.cds $species/$species.renamed.cds #Remove empty sequences
transeq $species/$species.renamed.cds $species/$species.prot ; sed -i 's/_1$//g' $species/$species.prot #translate all CDS
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $species/$species.prot | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\*/)' | tr "\t" "\n" > $species/$species.prot.nostop #Remove seq with internal stop codons
