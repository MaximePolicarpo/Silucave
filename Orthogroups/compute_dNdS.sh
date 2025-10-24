#!/bin/bash


LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
module load R/4.2.1-foss-2022a
module load MAFFT


curr_OGG=$1 
ortho_file=$2

IFS=$'\n'



dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"

echo "Computing stats on $curr_OGG"



##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
## === CODING SEQUENCES
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


module purge ; module load SAMtools

#Extract the list of sequences present in the OGG
grep "$curr_OGG	" $ortho_file | sed 's/, /,/g' | cut -f4- | tr ',' '\n' | tr '\t' '\n' | grep -v "$curr_OGG" | sed '/^$/d' | sed 's/_1$//g' > $curr_OGG.list
sed $'s/[^[:print:]\t]//g' $curr_OGG.list > $curr_OGG.list.reformat ; rm $curr_OGG.list

#Extract the corresponding coding sequences 
xargs samtools faidx concatenated_CDS.fa < $curr_OGG.list.reformat > Coding_sequences_alignments/$curr_OGG.cds


module purge ; conda activate HyPhy_Env


#Remove stop codons
hyphy cln Universal Coding_sequences_alignments/$curr_OGG.cds "No/No" Coding_sequences_alignments/$curr_OGG.clean.cds ; sed -i 's/?//g' Coding_sequences_alignments/$curr_OGG.clean.cds
mv Coding_sequences_alignments/$curr_OGG.clean.cds Coding_sequences_alignments/$curr_OGG.cds
sed -i 's/---$//g' Coding_sequences_alignments/$curr_OGG.cds
sed -i 's/-$//g' Coding_sequences_alignments/$curr_OGG.cds
sed -i 's/--$//g' Coding_sequences_alignments/$curr_OGG.cds


conda activate dN_dS_env


#Translate coding sequences into proteins
transeq Coding_sequences_alignments/$curr_OGG.cds Coding_sequences_alignments/$curr_OGG.prot ; sed -i 's/_1$//g' Coding_sequences_alignments/$curr_OGG.prot

#Remove sequences with X characters -- create problems when performing alignment
module purge ; module load SAMtools

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Coding_sequences_alignments/$curr_OGG.prot | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\X/)' | tr "\t" "\n" > Coding_sequences_alignments/$curr_OGG.prot.clean
mv Coding_sequences_alignments/$curr_OGG.prot.clean Coding_sequences_alignments/$curr_OGG.prot
grep ">" Coding_sequences_alignments/$curr_OGG.prot | sed 's/>//g' > Coding_sequences_alignments/$curr_OGG.id
xargs samtools faidx Coding_sequences_alignments/$curr_OGG.cds < Coding_sequences_alignments/$curr_OGG.id  > Coding_sequences_alignments/$curr_OGG.cds.cleaned ; mv Coding_sequences_alignments/$curr_OGG.cds.cleaned Coding_sequences_alignments/$curr_OGG.cds
rm Coding_sequences_alignments/$curr_OGG.cds.fai

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Coding_sequences_alignments/$curr_OGG.cds | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\-/)' | tr "\t" "\n" > Coding_sequences_alignments/$curr_OGG.cds.cleaned
mv Coding_sequences_alignments/$curr_OGG.cds.cleaned Coding_sequences_alignments/$curr_OGG.cds
grep ">" Coding_sequences_alignments/$curr_OGG.cds | sed 's/>//g' > Coding_sequences_alignments/$curr_OGG.id
xargs samtools faidx Coding_sequences_alignments/$curr_OGG.prot < Coding_sequences_alignments/$curr_OGG.id  > Coding_sequences_alignments/$curr_OGG.prot.cleaned ; mv Coding_sequences_alignments/$curr_OGG.prot.cleaned Coding_sequences_alignments/$curr_OGG.prot
rm Coding_sequences_alignments/$curr_OGG.prot.fai

#Count the number of sequences remaining and the number of different species included in the sequence file

grep ">" Coding_sequences_alignments/$curr_OGG.prot | sed 's/>//g' | sort > $curr_OGG.list.reformat
nb_seq=`grep -c ">" Coding_sequences_alignments/$curr_OGG.cds`


#Lets make an alignment if there are at-least three sequences. If there are more than 50 sequences, then use mafft to make a template alignment first

if [ $nb_seq -ge 3 ] ; then

	if [ $nb_seq -le 50 ] ; then 

		muscle5.1 -align Coding_sequences_alignments/$curr_OGG.prot -output Coding_sequences_alignments/$curr_OGG.prot.aln
		trimal -in Coding_sequences_alignments/$curr_OGG.prot.aln -cons 100 -backtrans Coding_sequences_alignments/$curr_OGG.cds -out Coding_sequences_alignments/$curr_OGG.cds.aln
		trimal -in Coding_sequences_alignments/$curr_OGG.prot.aln -backtrans Coding_sequences_alignments/$curr_OGG.cds -automated1 -out Coding_sequences_alignments/$curr_OGG.cds.trimmed
		trimal -in Coding_sequences_alignments/$curr_OGG.prot.aln -automated1 -out Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln
	
	else 
	
		grep ">" Coding_sequences_alignments/$curr_OGG.prot | sed 's/>//g' | shuf -n50 | sort > template.$curr_OGG.id 
		xargs samtools faidx Coding_sequences_alignments/$curr_OGG.prot < template.$curr_OGG.id > Coding_sequences_alignments/$curr_OGG.template.prot
		comm -23 $curr_OGG.list.reformat template.$curr_OGG.id  > remaining.$curr_OGG.id ; rm template.$curr_OGG.id
		xargs samtools faidx Coding_sequences_alignments/$curr_OGG.prot < remaining.$curr_OGG.id > Coding_sequences_alignments/$curr_OGG.toadd.prot
		rm remaining.$curr_OGG.id
		
		muscle5.1 -align Coding_sequences_alignments/$curr_OGG.template.prot -output Coding_sequences_alignments/$curr_OGG.template.aln
		
		module purge ; module load MAFFT

		mafft --add Coding_sequences_alignments/$curr_OGG.toadd.prot --keeplength Coding_sequences_alignments/$curr_OGG.template.aln > Coding_sequences_alignments/$curr_OGG.prot.aln
		rm Coding_sequences_alignments/$curr_OGG.toadd.prot ; rm Coding_sequences_alignments/$curr_OGG.template.aln ; rm Coding_sequences_alignments/$curr_OGG.template.prot
		trimal -in Coding_sequences_alignments/$curr_OGG.prot.aln -cons 100 -backtrans Coding_sequences_alignments/$curr_OGG.cds -out Coding_sequences_alignments/$curr_OGG.cds.aln
		trimal -in Coding_sequences_alignments/$curr_OGG.prot.aln -backtrans Coding_sequences_alignments/$curr_OGG.cds -automated1 -out Coding_sequences_alignments/$curr_OGG.cds.trimmed
		trimal -in Coding_sequences_alignments/$curr_OGG.prot.aln -automated1 -out Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln
	
	fi
fi



#Now lets make a maximum likelihood tree from the alignment (with the best model and four gamma rates)

nb_seq=`grep -c ">" Coding_sequences_alignments/$curr_OGG.prot.aln`

if [ $nb_seq -ge 3 ] ; then
	if [ $nb_seq -ge 4 ] ; then 
		iqtree -s Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln --seqtype AA -m TEST -mrate G4 -nt 8 -bb 1000
	else
		iqtree -s Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln --seqtype AA -m TEST -mrate G4 -nt 8
	fi
fi



#Now compute the dN/dS per branch using HyPhy


module purge ; conda activate HyPhy_Env

#FitMG94.bf can be downloaded from HyPhy github
hyphy FitMG94.bf --alignment Coding_sequences_alignments/$curr_OGG.cds.trimmed --tree Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.treefile --code Universal --type local ENV=TOLERATE_NUMERICAL_ERRORS=1




grep "LB\":" Coding_sequences_alignments/$curr_OGG.cds.trimmed.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > Coding_sequences_alignments/$curr_OGG.curr_LB_values.txt
grep "MLE\":" Coding_sequences_alignments/$curr_OGG.cds.trimmed.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > Coding_sequences_alignments/$curr_OGG.curr_MLE_values.txt
grep "UB\":" Coding_sequences_alignments/$curr_OGG.cds.trimmed.FITTER.json | sed 's/.*://g' | sed 's/,.*//g' > Coding_sequences_alignments/$curr_OGG.curr_UB_values.txt
grep "\"dN\"" Coding_sequences_alignments/$curr_OGG.cds.trimmed.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > Coding_sequences_alignments/$curr_OGG.curr_dN_values.txt
grep "\"dS\"" Coding_sequences_alignments/$curr_OGG.cds.trimmed.FITTER.json  | grep -v "{" | sed 's/.*://g' | sed 's/,//g' > Coding_sequences_alignments/$curr_OGG.curr_dS_values.txt
grep -B2 "LB\":" Coding_sequences_alignments/$curr_OGG.cds.trimmed.FITTER.json | grep -v "\-\-" | grep -v "Confidence" | grep -v "LB\":"  | sed 's/\"//g' | sed 's/:.*//g' | sed 's/^ *//g' > Coding_sequences_alignments/$curr_OGG.curr_labels
paste -d "," Coding_sequences_alignments/$curr_OGG.curr_labels Coding_sequences_alignments/$curr_OGG.curr_LB_values.txt Coding_sequences_alignments/$curr_OGG.curr_MLE_values.txt Coding_sequences_alignments/$curr_OGG.curr_UB_values.txt Coding_sequences_alignments/$curr_OGG.curr_dN_values.txt Coding_sequences_alignments/$curr_OGG.curr_dS_values.txt | sed "s/^/$curr_OGG,/g"  > Coding_sequences_alignments/$curr_OGG.omega.csv 
rm Coding_sequences_alignments/$curr_OGG.curr_LB_values.txt ; rm Coding_sequences_alignments/$curr_OGG.curr_MLE_values.txt ; rm Coding_sequences_alignments/$curr_OGG.curr_UB_values.txt ; rm Coding_sequences_alignments/$curr_OGG.curr_labels ; rm Coding_sequences_alignments/$curr_OGG.curr_dN_values.txt ; rm Coding_sequences_alignments/$curr_OGG.curr_dS_values.txt



#Remove unecessary files


rm Coding_sequences_alignments/$curr_OGG.id
rm Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.bionj
rm Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.ckp.gz
rm Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.contree
rm Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.iqtree
rm Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.mldist
rm Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.model.gz
rm Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.splits.nex
rm Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.uniqueseq.phy
rm Coding_sequences_alignments/$curr_OGG.prot.trimmed.fai
rm $curr_OGG.list.reformat



echo "Orthogroups dN/dS computed correctly"

dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"

