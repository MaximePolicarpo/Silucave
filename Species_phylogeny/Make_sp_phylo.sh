#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Extract BUSCO gene sequences ####################################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

#Lets extract all busco genes present in at-least half of the species

mkdir Common_busco_genes
for species in `cat species_list.txt` ; do
	ls -l BUSCO_$species/run_actinopterygii_odb10/busco_sequences/single_copy_busco_sequences/ | sed 's/.* //g' | grep ".faa$" > Common_busco_genes/$species.singlecopy.list
done

#Extract busco genes present in at-least 4 species
cat Common_busco_genes/*.list | grep ".faa" | sort | uniq -c | sed 's/^ *//g' | sed 's/ /,/g' | awk -F ',' '{ if ($1 >= 4) print $2 }'  > BUSCO_list_half_sp.txt 

#Align BUSCO genes

mkdir Alignements
for i in `cat BUSCO_list_half_sp.txt` ; do sbatch --qos=30min -c 1 --mem=2G concat_gene.sh $i ; done

cd Alignements/

for i in *.faa ; do sbatch --qos=30min -c 5 --mem=20G --job-name=muscle_aln --wrap="muscle5.1 -align $i -output $i.aln" ; done


#Trim all alignments
for i in *.aln ; do sbatch --qos=30min -c 4 --mem=6G --job-name=trimal_aln lauch_trimal.sh $i ; done
sed -i 's/?/-/g' *.trimal

#Make a maximum likelihood tree with each gene and with the best model found by ModelFinder

for gene_alignment in *.trimal ; do 
	sbatch --qos=6hours -c 10 --mem=30G iqtree_launch.sh $gene_alignment
done


#Collapse every nodes that have a low bootstrap value (10%)

for gene_tree in *.treefile ; do 
	sbatch --qos=30min -c 2 --mem=4G collapse_lowsupport.sh $gene_tree $gene_tree.collapsed
done


#Concatenate all trimed alignments

python3 AMAS/amas/AMAS.py concat -f fasta -d aa -i Alignements/*.trimal --part-format nexus

sed -i 's/?/-/g' concatenated.out
mv concatenated.out AMAS_concatenated_alignment_BUSCO.fa
mv partitions.txt AMAS_concatenated_alignment_BUSCO.partition.nexus


#Create a file with three calibration dates
#nano Actino_calibration.txt : 
#Danio_rerio,Ictalurus_punctatus	-142
#Ictalurus_punctatus,Ancistrus_triradiatus	-96
#Ictalurus_punctatus,Clarias_gariepinus	-80


## run ASTRAL

cat Alignements/*.collapsed > All_genes_tree.concatnwk 
sbatch --qos=6hours -c 5 --mem=50G launch_astral.sh  


#Unroot astral tree
>R 
>library('ape')
>library("phytools")
>mytree <- read.tree("Astral_nolength.nwk")
>unrooted_tree <- unroot(mytree)
>write.tree(unrooted_tree, "Astral_unrooted_nolength.nwk")

#replace null branch length by 0.01
sed -i 's/NaN/0.01/g' Astral_unrooted_nolength.nwk

#date using the least-squared method

sbatch --qos=6hours -c 30 --mem=80G date_tree.sh

## Add node labels to the ASTRAL tree

>R
>library("ape")
>mytree <- read.tree("AMAS_concatenated_alignment_BUSCO.fa.timetree.nwk")
>mytree_nodelabel <- makeNodeLabel(mytree, method = "number", prefix = "Node")
>write.tree(mytree_nodelabel, file="AMAS_concatenated_alignment_BUSCO.fa.timetree.nwk.nodelabel")



#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################

# Test if the topology is the same usign iqtree and a a concatenaed BUSCO alignment

mkdir IQTREE_topology ; cd IQTREE_topology

cp ../AMAS_concatenated_alignment_BUSCO.fa ./

sbatch --qos=1day -c 10 --mem=40G iqtree_concat.sh 

