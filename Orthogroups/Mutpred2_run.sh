#!/bin/bash


#SBATCH --job-name=Mutpred2   # Job name



curr_OGG=$1

#First, reconstruct the ancestral amino acid sequences in the trees. /scicore/home/salzburg/polica0000/SiluCave/OrthoFinder/aaml_control_file.ctl  => Empricial+F + LG substitution matrix

eval "$(conda shell.bash hook)"
conda activate PAML
module purge ; module load R/4.4.1-foss-2023b

Rscript Remove_nodelabels.R Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.treefile Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.treefile.noNodeLabels


rm -rf temp_dir.$curr_OGG/ ; mkdir temp_dir.$curr_OGG 

cp Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln temp_dir.$curr_OGG/
cp Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.treefile.noNodeLabels temp_dir.$curr_OGG/
sed "s/alignment.aln/$curr_OGG.prot.trimmed.aln/g" aaml_control_file.ctl | sed "s/treefile.nwk/$curr_OGG.prot.trimmed.aln.treefile.noNodeLabels/g" | sed "s/output.aaml/$curr_OGG.aaml/g" > temp_dir.$curr_OGG/aaml_control_file.$curr_OGG.ctl

cd temp_dir.$curr_OGG/ ; codeml aaml_control_file.$curr_OGG.ctl ; cd ../


mv temp_dir.$curr_OGG/rst Coding_sequences_alignments/$curr_OGG.aaml.rst
mv temp_dir.$curr_OGG/$curr_OGG.aaml Coding_sequences_alignments/$curr_OGG.aaml

rm -rf temp_dir.$curr_OGG/



##Now convert AAML results to Mutpred2 inputs with R

Rscript Prepare_mutpred_input.R Coding_sequences_alignments/$curr_OGG.aaml.rst Coding_sequences_alignments/$curr_OGG.aaml.summary


#remove branches with no substitutions 

grep -v ",," Coding_sequences_alignments/$curr_OGG.aaml.summary > Coding_sequences_alignments/$curr_OGG.aaml.summary.temp ; mv Coding_sequences_alignments/$curr_OGG.aaml.summary.temp Coding_sequences_alignments/$curr_OGG.aaml.summary


#Make an mutpred2 input file and verify the correctness of substitutions positions

[ -f "Coding_sequences_alignments/$curr_OGG.mutpred.input" ] && rm "Coding_sequences_alignments/$curr_OGG.mutpred.input"

IFS=$'\n' 
for line in `cat Coding_sequences_alignments/$curr_OGG.aaml.summary` ; do
	tip_label=`echo "$line" | cut -f1 -d ","`
	list_AA_subst=`echo "$line" | cut -f2 -d "," | sed 's/\///g'`
	anc_seq_label=`echo "$line" | cut -f3 -d ","`
	anc_seq=`grep "node #$anc_seq_label" Coding_sequences_alignments/$curr_OGG.aaml.rst | tail -1 | sed "s/node #$anc_seq_label//g" | sed "s/^ *//g" | sed 's/ //g'`

	echo "$list_AA_subst" | tr ' ' "\n"  | grep "[A-z]" | sed 's/\([A-Z]\)\([0-9]\+\)\([A-Z]\)/\1,\2,\3/' > Coding_sequences_alignments/$curr_OGG.currsubt

	for curr_subst in `cat Coding_sequences_alignments/$curr_OGG.currsubt` ; do
		curr_AA=`echo "$curr_subst" | cut -f1 -d ","`
		curr_pos=`echo "$curr_subst" | cut -f2 -d ","`
		curr_pos_min1=$(( curr_pos - 1 ))
		curr_pos_AA=${anc_seq:$curr_pos_min1:1}

		if [ $curr_pos_AA == $curr_AA ] ; then
			echo "$curr_subst" >> Coding_sequences_alignments/$curr_OGG.currsubt.filt
		fi

	done

	list_AA_subst_filt=`sed 's/,//g' Coding_sequences_alignments/$curr_OGG.currsubt.filt | tr "\n" " " | sed 's/ $/\n/'`

	echo ">$tip_label $list_AA_subst_filt" >> Coding_sequences_alignments/$curr_OGG.mutpred.input
	echo "$anc_seq" >> Coding_sequences_alignments/$curr_OGG.mutpred.input

	rm Coding_sequences_alignments/$curr_OGG.currsubt ; rm Coding_sequences_alignments/$curr_OGG.currsubt.filt

done




##Finally, run Mutpred2 ! :) 

module purge 
conda deactivate
conda activate mutpred2

run_mutpred2.sh -i Coding_sequences_alignments/$curr_OGG.mutpred.input -o Coding_sequences_alignments/$curr_OGG.mutpred.results

