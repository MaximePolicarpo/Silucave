#!/bin/bash


#SBATCH --job-name=RELAX   # Job name


eval "$(conda shell.bash hook)"
conda activate HyPhy_Env

IFS=$'\n'
curr_line=$1
curr_test_branch=`echo "$curr_line" | cut -f2 -d ","`
curr_OGG=`echo "$curr_line" | cut -f1 -d ","`

mkdir Relax_Targeted_Results/$curr_test_branch/

cp Coding_sequences_alignments/$curr_OGG.cds.trimmed Relax_Targeted_Results/$curr_test_branch/$curr_OGG.cds.trimmed
cp Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.treefile Relax_Targeted_Results/$curr_test_branch/$curr_OGG.prot.trimmed.aln.treefile


module purge ; module load R/4.4.1-foss-2023b

Rscript assign_node_as_ref.R Relax_Targeted_Results/$curr_test_branch/$curr_OGG.prot.trimmed.aln.treefile Relax_Targeted_Results/$curr_test_branch/$curr_OGG.prot.trimmed.aln.treefile.marked


grep "$curr_OGG" aBSREL_df_significant.csv | cut -f2 -d "," > Relax_Targeted_Results/$curr_test_branch/$curr_OGG.POS.tiplabs
grep ">" Relax_Targeted_Results/$curr_test_branch/$curr_OGG.cds.trimmed | sed 's/>//g' > Relax_Targeted_Results/$curr_test_branch/$curr_OGG.ALL.tiplabs

for pos_gene in `cat Relax_Targeted_Results/$curr_test_branch/$curr_OGG.POS.tiplabs` ; do 
	grep -v "$pos_gene" Relax_Targeted_Results/$curr_test_branch/$curr_OGG.ALL.tiplabs > Relax_Targeted_Results/$curr_test_branch/temp 
	mv Relax_Targeted_Results/$curr_test_branch/temp Relax_Targeted_Results/$curr_test_branch/$curr_OGG.ALL.tiplabs
done

for all_gene in `cat Relax_Targeted_Results/$curr_test_branch/$curr_OGG.ALL.tiplabs` ; do
	sed -i "s/$all_gene/$all_gene {reference}/g" Relax_Targeted_Results/$curr_test_branch/$curr_OGG.prot.trimmed.aln.treefile.marked
done


sed -i "s/$curr_test_branch/$curr_test_branch {test}/g" Relax_Targeted_Results/$curr_test_branch/$curr_OGG.prot.trimmed.aln.treefile.marked

hyphy relax --alignment Relax_Targeted_Results/$curr_test_branch/$curr_OGG.cds.trimmed --tree Relax_Targeted_Results/$curr_test_branch/$curr_OGG.prot.trimmed.aln.treefile.marked --code Universal --test test --reference reference ENV=TOLERATE_NUMERICAL_ERRORS=1
