###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###=============================== Launch OrthoFinder  =============================================================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


#Put all proteome fasta files in the same folder
mkdir Siluriformes_Proteomes

for species in `cat species_list.txt` ; do
	cp Proteins_woTE/$species.prot Siluriformes_Proteomes/$species.fa
done

#Run OrthoFinder usign the species phylogeny computed with ASTRAL

sbatch --qos=1week -c 50 --mem=200G --job-name=Ortho -e error_ortho.out -o slurm_ortho.out launch_orthofinder.sh AMAS_concatenated_alignment_BUSCO.fa.timetree.pruned.nwk Siluriformes_Proteomes


###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###=============================== Align all orthogroups  =============================================================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


#First with only RefSeq species

mkdir Coding_sequences_alignments
mkdir Log_files


ortho_file=Siluriformes_Proteomes/OrthoFinder/Results_Oct11/Phylogenetic_Hierarchical_Orthogroups/N1.tsv
cut -f1 $ortho_file | tail -n+2 > list_orthogroups.txt
split -l 10000 --numeric-suffixes list_orthogroups.txt splitted_list_orthogroups_


cat Coding_sequences_woTE/*.cds > concatenated_CDS.fa


#Align all orthogroups, make maximum likelihood phylogenies, and compute dn/ds per branch

for curr_OGG in `cat list_orthogroups.txt` ; do sbatch --qos=6hours -c 8 --mem=20G compute_dNdS.sh $curr_OGG Siluriformes_Proteomes/OrthoFinder/Results_Oct11/Phylogenetic_Hierarchical_Orthogroups/N1.tsv ; done


###Extract and combine dN/dS results

rm dNdS_per_orthogroup.csv
for curr_OGG in `cat list_orthogroups.txt` ; do 
	grep -v ",Node" Coding_sequences_alignments/$curr_OGG.omega.csv >> dNdS_per_orthogroup.csv
done



###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###=============================== Reduce orthogroups to keep only one gene per species  ============================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================

#To use RERconverge, we need orthogroups with only one sequence per species
#We thus need to reduce orthogroups sizes.
#For that, one solution is to keep only the gene per species with the lowest phylogenetic distance from the root (see https://academic.oup.com/mbe/article/37/7/2052/5804990)

### First, root at midpoint every ML trees computed with IQ-tree

for curr_OGG in `cat list_orthogroups.txt` ; do 
	if [ -f Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.treefile ]; then echo "$curr_OGG" >> list_orthogroups.toroot.txt ; fi
done


for curr_OGG in `cat list_orthogroups.toroot.txt` ; do 
	sbatch -qos=30min -c 4 --mem=8G root_at_midpoint.sh Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.treefile Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.treefile.midpointRooted
done


#For species with more than 1 gene, remove extra gene(s) from the alignment, and compute a ML tree with the same topology as the species tree using phangorn


for curr_OGG in `cat list_orthogroups.toroot.txt` ; do 
	sbatch -qos=30min -c 4 --mem=8G Keep_onetip_persp.sh Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.treefile.midpointRooted Coding_sequences_alignments/$curr_OGG.OnetipPersp.txt $curr_OGG
done



###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###=============================== Launch RERconverge  ==============================================================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================

#First combine all newick files into a single file

for curr_OGG in `cat list_orthogroups.toroot.txt` ; do
	newick_tree=$(cat "Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.OneTipPerSp.nwk")
	echo -e "$curr_OGG\t$newick_tree" >> AllTrees.Phangorn.nwk
done

awk -F'\t' '$2 != ""' AllTrees.Phangorn.nwk > temp ; mv temp AllTrees.Phangorn.nwk


## Assign goterms to each orthogroup

#First, keep the longest gene of each orthogroup

for curr_OGG in `cat list_orthogroups.toroot.txt` ; do sbatch -qos=30min -c 4 --mem=8G Extract_reprez.sh $curr_OGG ; done

cat Coding_sequences_alignments/*.reprez_seq.prot > Representative_proteins_orthogroups.fa 


#Run interproscan with goterm database

mkdir out_interpro
module load Java/21.0.2

./interproscan-5.70-102.0/interproscan.sh -i Representative_proteins_orthogroups.fa  -d out_interpro/ -t p -appl PANTHER --goterms -f TSV -cpu 20 -T temp.interpro/ --verbose
rm -r temp.interpro

cut -f1,14 out_interpro/Representative_proteins_orthogroups.fa.tsv | grep "GO:" | awk -F'\t' '{split($2, go_terms, "|"); for (i in go_terms) { print $1 "\t" go_terms[i]; }}' | sort | uniq | sed 's/(.*//g' > OGG_GOterms.tsv
cut -f1,5,13 out_interpro/Representative_proteins_orthogroups.fa.tsv > OGG_Panther_Family.tsv

#Assign gene names based on the best BLASTP match against a database of ray-finned fishes genes from RefSeq

for species in `cat RefSeq_species_list.txt` ; do
	cat $species/$species.prot.nostop >> RefSeq_genes.prot #The prot.nostop files are generated from GFF3 files of RefSeq species using agat. See the "Reads mapping and processing" folder
done

module purge ; module load DIAMOND

diamond makedb --in RefSeq_genes.prot -d RefSeq_genes --threads 20

diamond blastp -p 20 --query Representative_proteins_orthogroups.fa  --max-target-seqs 1 --out ReprezSeqs_vs_RefSeq.blastp --db RefSeq_genes --outfmt 6

grep ">" Representative_proteins_orthogroups.fa | sed 's/>//g' > Representative_proteins_orthogroups.id

for gene in `cat Representative_proteins_orthogroups.id` ; do

	evalue=`grep -m1 "$gene" ReprezSeqs_vs_RefSeq.blastp | cut -f11`
	best_match=`grep -m1 "$gene" ReprezSeqs_vs_RefSeq.blastp | cut -f2`
	species_best_match=`echo "$best_match" | sed 's/---.*//g'`
	best_match_gene=`echo "$best_match" | sed 's/.*---//g'`
	GFF_file=`ls -l $species_best_match/ | grep ".gff$" | sed 's/.* //g'`
	
	product=`grep "$best_match_gene" $species_best_match/$GFF_file | grep -m1 "	mRNA	\|	gene	\|	V_gene_segment 	\|	C_gene_segment	" | sed 's/.*product=//g' | sed 's/.*standard_name=//g' | sed 's/;.*//g'`
	gene_name=`grep "$best_match_gene" $species_best_match/$GFF_file | grep -m1 "	mRNA	\|	gene	\|	V_gene_segment	\|	C_gene_segment	" | sed 's/.*Parent=//g' | sed 's/;.*//g'`

	echo "$gene	$gene_name	$product	$evalue" >> OGG_Names_Product.tsv
done


#Assign KEGG pathways to each gene

#=> https://www.kegg.jp/ghostkoala/ ==> submit sequences and download annotations results
awk '{if (NF == 1) print $0, "\tnoKO"; else print $0}' user_ko.txt > user_ko.filled.txt


#Make a nice table og GO terms

cut -f2 OGG_GOterms.tsv | sort | uniq > uniq_GO.txt

rm goterms_annot.gmt
for curr_GO in `cat uniq_GO.txt` ; do

	if grep -q "$curr_GO" goterm_desc.csv ; then curr_desc=`grep "$curr_GO" goterm_desc.csv | cut -f2 | sed 's/\///g'` ; else echo "$curr_GO" ; fi
	if grep -q "$curr_GO" goterm_class.csv  ; then curr_class=`grep "$curr_GO" goterm_class.csv | cut -f2` ; else echo "$curr_GO" ; fi
	grep "$curr_GO" OGG_GOterms.tsv | cut -f1 | tr "\n" "\t" | sed '$a\' | sed "s/^/$curr_GO.$curr_desc	$curr_class	/g" >> goterms_annot.gmt
done

grep "biological_process" goterms_annot.gmt > goterms_annot.gmt.BP
grep "cellular_component" goterms_annot.gmt > goterms_annot.gmt.CP
grep "molecular_function" goterms_annot.gmt > goterms_annot.gmt.MF

grep "^GO" goterm_desc.csv  > temp ; mv temp goterm_desc.csv
sed -i "s/'//g" goterm_desc.csv
sed -i "s/(//g" goterm_desc.csv
sed -i "s/)//g" goterm_desc.csv


#Make a nice table for PANTHER families

cut -f2 OGG_Panther_Family.tsv | sort | uniq > uniq_PANTHER.txt

for curr_panther in `cat uniq_PANTHER.txt` ; do

	curr_desc=`grep -m1 "$curr_panther" OGG_Panther_Family.tsv | cut -f3 | sed 's/\///g'`

	grep "$curr_panther" OGG_Panther_Family.tsv | cut -f1 | tr "\n" "\t" | sed '$a\' | sed "s/^/$curr_panther.$curr_desc	PANTHER	/g" >> panther_annot.gmt
done



#Make a nice GMT to use in gProfiler2

cp goterms_annot.gmt.BP goterms_annot.BP.gprofiler.gmt
cp goterms_annot.gmt.CP goterms_annot.CC.gprofiler.gmt
cp goterms_annot.gmt.MF goterms_annot.MF.gprofiler.gmt
cp panther_annot.gmt panther_annot.gprofiler.gmt
cp KEGG.all.HOG.gmt KEGG.gprofiler.gmt


#Now launch RERconverge !

Rscript RER_siluriformes.R

###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###=============================== Launch aBSREL and RELAX  =========================================================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================

#Keep trees with at-least one cavefish + one surface fish branch 


for curr_OGG in `cat list_orthogroups.toroot.txt` ; do 
	nb_cave_branch=`grep ">" Coding_sequences_alignments/$curr_OGG.cds.trimmed | grep "Prietella_phreatophila\|CHM6\|CSV83\|CUL4\|CUL9\|Trichomycterus_rosablanca" | wc -l`
	nb_surface_branch=`grep ">" Coding_sequences_alignments/$curr_OGG.cds.trimmed | grep -v "Prietella_phreatophila\|CHM6\|CSV83\|CUL4\|CUL9\|Trichomycterus_rosablanca" | wc -l`
	if [ $nb_cave_branch -ge 1 ] && [ $nb_surface_branch -ge 1 ] ; then echo "$curr_OGG" >> list_orthogroups.absrel_relax.txt ; fi
done


for curr_OGG in `cat list_orthogroups.absrel_relax.txt` ; do 
	sbatch -qos=6hours -c 4 --mem=8G aBSREL_RELAX.sh Coding_sequences_alignments/$curr_OGG.cds.trimmed Coding_sequences_alignments/$curr_OGG.prot.trimmed.aln.treefile $curr_OGG
done


### Extract aBSREL results

IFS=$'\n'
rm aBSREL_results.csv
for curr_OGG in `cat list_orthogroups.absrel_relax.txt` ; do 


	if [ -f Coding_sequences_alignments/$curr_OGG.cds.trimmed.ABSREL.json ]; then

		grep ':\"test\"' Coding_sequences_alignments/$curr_OGG.cds.trimmed.ABSREL.json | sed 's/:.*//g' | sed 's/^ *//g' | sed 's/"//g' > test_branches.txt 
		
		for branch in `cat test_branches.txt` ; do 
			LRT=`grep -A15 "$branch\":{" Coding_sequences_alignments/$curr_OGG.cds.trimmed.ABSREL.json | grep "LRT" | sed "s/.*://g" | sed 's/,$//g'`
			corr_pvalue=`grep -A15 "$branch\":{" Coding_sequences_alignments/$curr_OGG.cds.trimmed.ABSREL.json | grep "Corrected P-value" | sed "s/.*://g" | sed 's/,$//g'`
			uncorr_pvalue=`grep -A15 "$branch\":{" Coding_sequences_alignments/$curr_OGG.cds.trimmed.ABSREL.json | grep "Uncorrected P-value" | sed "s/.*://g" | sed 's/,$//g'`
		
			echo "$curr_OGG,$branch,$LRT,$corr_pvalue,$uncorr_pvalue" >> aBSREL_results.csv
		
		done

	fi
done


### Extract RELAX results

IFS=$'\n'
rm RELAX_results.csv
for curr_OGG in `cat list_orthogroups.absrel_relax.txt` ; do 

	if [ -f Coding_sequences_alignments/$curr_OGG.cds.trimmed.RELAX.json ]; then

		LRT=`grep "LRT" Coding_sequences_alignments/$curr_OGG.cds.trimmed.RELAX.json | sed 's/.*://g' | sed "s/,//g"`
		K_param=`grep "relaxation or intensification parameter" Coding_sequences_alignments/$curr_OGG.cds.trimmed.RELAX.json | sed 's/.*://g'`

		echo "$curr_OGG,$LRT,$K_param" >> RELAX_results.csv
	fi


done


###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###=============================== Launch targeted RELAX analysis ===================================================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


#For each branch detected under episodic diversification with aBSREL, we will test if they are under relaxed selection with RELAX

rm -rf Log_RELAX ;  mkdir Log_RELAX
rm -rf Relax_Targeted_Results ;  mkdir Relax_Targeted_Results


for curr_line in `cat aBSREL_df_significant.csv` ; do 
	sbatch -qos=6hours -c 4 --mem=8G Targeted_RELAX.sh $curr_line
done




### Extract Cand RELAX results
IFS=$'\n'
rm RELAX_results.Cand.csv
for curr_line in `cat aBSREL_df_significant.csv` ; do
	curr_test_branch=`echo "$curr_line" | cut -f2 -d ","`
	curr_OGG=`echo "$curr_line" | cut -f1 -d ","`
	LRT=`grep "LRT" Relax_Targeted_Results/$curr_test_branch/$curr_OGG.cds.trimmed.RELAX.json | sed 's/.*://g' | sed "s/,//g"`
	K_param=`grep "relaxation or intensification parameter" Relax_Targeted_Results/$curr_test_branch/$curr_OGG.cds.trimmed.RELAX.json | sed 's/.*://g'`
	echo "$curr_OGG,$curr_test_branch,$LRT,$K_param" >> RELAX_results.Cand.csv
done



###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###=============================== Launch cSUBST  ===================================================================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


conda activate csubst_env

#Run csubst on every orthogroups to detect convergence between cave species. We will run it with both the gene tree AND species tree

#Only keep OGG with at-least two cave branches and two surface branches

for curr_OGG in `cat list_orthogroups.toroot.txt` ; do 
	nb_cave_branch=`grep ">" Coding_sequences_alignments/$curr_OGG.cds.trimmed | grep "Prietella_phreatophila\|CHM6\|CSV83\|CUL4\|CUL9\|Trichomycterus_rosablanca" | wc -l`
	nb_surface_branch=`grep ">" Coding_sequences_alignments/$curr_OGG.cds.trimmed | grep -v "Prietella_phreatophila\|CHM6\|CSV83\|CUL4\|CUL9\|Trichomycterus_rosablanca" | wc -l`
	if [ $nb_cave_branch -ge 2 ] && [ $nb_surface_branch -ge 2 ] ; then echo "$curr_OGG" >> list_orthogroups.csubst.txt ; fi
done

#Launch Csubst ! 

mkdir Csubst_results_GeneTrees
mkdir Csubst_results_SpeciesTrees


for curr_OGG in `cat list_orthogroups.csubst.txt` ; do sbatch -qos=6hours -c 4 --mem=8G csubst_launch.sh $curr_OGG ; done


#Extract Csubst results

rm all_OGG_csubst_cb_2.GeneTree.tsv ; rm all_OGG_csubst_cb_3.GeneTree.tsv ; rm all_OGG_csubst_b.GeneTree.tsv
rm all_OGG_csubst_cb_2.SpeciesTree.tsv ; rm all_OGG_csubst_cb_3.SpeciesTree.tsv ; rm all_OGG_csubst_b.SpeciesTree.tsv

for curr_OGG in `cat list_orthogroups.csubst.txt` ; do 

	if [ -f Csubst_results_GeneTrees/$curr_OGG.csubst_cb_2.tsv ] ; then sed "s/^/$curr_OGG	/g" Csubst_results_GeneTrees/$curr_OGG.csubst_cb_2.tsv >> all_OGG_csubst_cb_2.GeneTree.tsv ; fi
	if [ -f Csubst_results_GeneTrees/$curr_OGG.csubst_cb_3.tsv ] ; then sed "s/^/$curr_OGG	/g" Csubst_results_GeneTrees/$curr_OGG.csubst_cb_3.tsv >> all_OGG_csubst_cb_3.GeneTree.tsv ; fi
	sed "s/^/$curr_OGG	/g" Csubst_results_GeneTrees/$curr_OGG.csubst_b.tsv >> all_OGG_csubst_b.GeneTree.tsv

	if [ -f Csubst_results_SpeciesTrees/$curr_OGG.csubst_cb_2.tsv ] ; then sed "s/^/$curr_OGG	/g" Csubst_results_SpeciesTrees/$curr_OGG.csubst_cb_2.tsv >> all_OGG_csubst_cb_2.SpeciesTree.tsv ; fi
	if [ -f Csubst_results_SpeciesTrees/$curr_OGG.csubst_cb_3.tsv ] ; then sed "s/^/$curr_OGG	/g" Csubst_results_SpeciesTrees/$curr_OGG.csubst_cb_3.tsv >> all_OGG_csubst_cb_3.SpeciesTree.tsv ; fi
	sed "s/^/$curr_OGG	/g" Csubst_results_SpeciesTrees/$curr_OGG.csubst_b.tsv >> all_OGG_csubst_b.SpeciesTree.tsv

done


grep "branch_name" all_OGG_csubst_b.GeneTree.tsv | head -1 | sed 's/N1.*	branch_name/OGG	branch_name/g' > header_1
grep -v "branch_name" all_OGG_csubst_b.GeneTree.tsv > temp ; cat header_1 > all_OGG_csubst_b.GeneTree.tsv ; cat temp >> all_OGG_csubst_b.GeneTree.tsv ; rm temp ; rm header_1
grep "branch_name" all_OGG_csubst_b.SpeciesTree.tsv | head -1 | sed 's/N1.*	branch_name/OGG	branch_name/g' > header_1
grep -v "branch_name" all_OGG_csubst_b.SpeciesTree.tsv > temp ; cat header_1 > all_OGG_csubst_b.SpeciesTree.tsv ; cat temp >> all_OGG_csubst_b.SpeciesTree.tsv ; rm temp ; rm header_1

head -1 Csubst_results_GeneTrees/N1.HOG0000044.csubst_cb_2.tsv | sed 's/^branch_id_1/OGG	branch_id_1/g' > header_1
grep -v "branch_id_1" all_OGG_csubst_cb_2.GeneTree.tsv > temp ; cat header_1 > all_OGG_csubst_cb_2.GeneTree.tsv ; cat temp >> all_OGG_csubst_cb_2.GeneTree.tsv ; rm temp ; rm header_1
grep "branch_id_1" all_OGG_csubst_cb_3.GeneTree.tsv | head -1 | sed 's/N1.*	branch_id_1/OGG	branch_id_1/g' > header_1
grep -v "branch_id_1" all_OGG_csubst_cb_3.GeneTree.tsv > temp ; cat header_1 > all_OGG_csubst_cb_3.GeneTree.tsv ; cat temp >> all_OGG_csubst_cb_3.GeneTree.tsv ; rm temp ; rm header_1

head -1 Csubst_results_SpeciesTrees/N1.HOG0000044.csubst_cb_2.tsv | sed 's/^branch_id_1/OGG	branch_id_1/g' > header_1
grep -v "branch_id_1" all_OGG_csubst_cb_2.SpeciesTree.tsv > temp ; cat header_1 > all_OGG_csubst_cb_2.SpeciesTree.tsv ; cat temp >> all_OGG_csubst_cb_2.SpeciesTree.tsv ; rm temp ; rm header_1
grep "branch_id_1" all_OGG_csubst_cb_3.SpeciesTree.tsv | head -1 | sed 's/N1.*	branch_id_1/OGG	branch_id_1/g' > header_1
grep -v "branch_id_1" all_OGG_csubst_cb_3.SpeciesTree.tsv > temp ; cat header_1 > all_OGG_csubst_cb_3.SpeciesTree.tsv ; cat temp >> all_OGG_csubst_cb_3.SpeciesTree.tsv ; rm temp ; rm header_1

grep "Node" all_OGG_csubst_cb_2.GeneTree.tsv


###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###=============================== Launch Mutpred2 on every genes  ==================================================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================

#Prepare files for mutpred2 and run mutpred2
for curr_OGG in `cat list_orthogroups.toroot.txt` ; do sbatch -qos=1day -c 4 --mem=8G Mutpred2_run.sh $curr_OGG ; done


## Extract mutpred2 results

rm All_Mutpred2_scores.csv
for curr_OGG in `cat list_orthogroups.toroot.txt` ; do 
	if [ -f Coding_sequences_alignments/$curr_OGG.mutpred.results ] ; then
		cut -f1,2,3 -d "," Coding_sequences_alignments/$curr_OGG.mutpred.results | grep -v ",#" | sed "s/^/$curr_OGG,/g" >> All_Mutpred2_scores.csv
	fi
done
grep -v "ID,Substitution,MutPred2" All_Mutpred2_scores.csv > temp ; mv temp All_Mutpred2_scores.csv


###==================================================================================================================================
###==================================================================================================================================
###=================================  Compute Grantham distances   ==================================================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================

for curr_OGG in `cat list_orthogroups.toroot.txt` ; do 
	grep ">" Coding_sequences_alignments/$curr_OGG.mutpred.input | sed 's/>//g' | awk '{for (i=2; i<=NF; i++) print $1, $i}' | sed "s/^/$curr_OGG /g" >> Table_all_substitutions.txt
done

split -l 10000 --numeric-suffixes Table_all_substitutions.txt splitted_Table_all_substitutions

for file in splitted_Table_all_substitutions* ; do
	output_name=`echo "$file" | sed 's/$/.grantham/g'`
	sbatch --qos=30min -c 4 --mem=20G compute_grantham_dist.sh $file $output_name
done


cat *.grantham > Table_all_substitutions.grantham.txt


