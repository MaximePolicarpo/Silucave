#!/bin/bash


#SBATCH --job-name=Parse_BRAKER   # Job name

species=$1

echo "$species"

cd $species.BRAKER/

samtools faidx braker.codingseq 

#First removal of non supported transcripts
singularity exec braker3.sif selectSupportedSubsets.py --fullSupport FULLSUPPORT --anySupport ANYSUPPORT --noSupport NOSUPPORT braker.gtf hintsfile.gff
sed 's/.*transcript_id//g' NOSUPPORT  | grep -v "#" | sed 's/ \"//g' | sed 's/\".*//g' | sort | uniq > genes_to_remove.txt
grep ">" braker.codingseq | sed 's/>//g' | sort | uniq > all_genes_id.txt
comm -23 all_genes_id.txt genes_to_remove.txt > all_genes_id_supported.txt

#Second keep longest transcript per gene
sed 's/\..*//g' all_genes_id_supported.txt | sort | uniq > uniq_genes_id.txt

rm uniq_transcript_id.txt
for gene in `cat uniq_genes_id.txt` ; do 

	grep "^$gene\." braker.codingseq.fai | sort -k2 -nr | head -1 | cut -f1 >> uniq_transcript_id.txt

done

xargs samtools faidx braker.codingseq  < uniq_transcript_id.txt > $species.codingseq.LI

transeq $species.codingseq.LI $species.prot.LI ; sed -i 's/_1$//g' $species.prot.LI 

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $species.prot.LI | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\*/)' | tr "\t" "\n" > $species.prot.LI.nostop

grep ">" $species.prot.LI.nostop | sed 's/>//g' | sort | uniq > good_id.txt
xargs samtools faidx $species.codingseq.LI < good_id.txt > $species.codingseq.LI.nostop ; rm good_id.txt

