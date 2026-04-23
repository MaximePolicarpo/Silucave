# Silucave

## I - Reads mapping and processing

The master script provided [Reads_processing_mapping.sh](Reads_mapping/Reads_processing_mapping.sh) allows to:

- Clean raw illumina reads using fastp
- Map reads to reference genomes using bwa
- Call variants from reads mapped
- Generate consensus genomes and coding sequences from reads mapped
- Download catfishes genomes available on Genbank/Refseq
- Assess genome completness using BUSCO
- Annotate genomes with BRAKER
- Estimate GC percentage and genome size from raw reads with GenomeScope 

All the scripts called in the master script can be retrieved in the Reads_mapping folder
Raw reads can be retrieved from the SRA database (Bioproject ID : PRJNA1348312)

## II - Siluriformes phylogeny

The master script provided [Make_sp_phylo.sh](Species_phylogeny/Make_sp_phylo.sh) allows to:

- Align BUSCO genes
- Compute maximum likelihood phylogenies of BUSCO genes
- Make a species tree using ASTRAL
- Make a species tree using a concatenated alignment of BUSCO genes, with IQ-TREE

All the scripts called in the master script can be retrieved in the Species_phylogeny folder


## III - Orthogroups analysis

The master script provided [Orthogroups_analysis.sh](Orthogroups/Orthogroups_analysis.sh) allows to:

- Run orthofinder to assign genes to orthogroups
- Align all orthogroups using MAFFT and MUSCLE
- Compute maximum likelihood phylogenies of orthogroups using IQ-TREE
- Compute the dN/dS per branch using HyPhy
- Estimate the relative evolutionary rate of each gene using phangorn and the species phylogeny
- Assign GOterms to orthogroups
- Launch RERconverge
- Launch ABSREl, RELAX, CSUBST and MUTPRED2 for each orthogroup.

Proteome fasta file for each species is available in this Github (splitted in three files due to GitHub size limits. To obtain the fasta :
cat Proteomes_part_* > Proteomes_recombined.zip ; unzip Proteomes_recombined.zip)

## IV - Data analysis and figures

The R script provided [Silucave_mainscript_analysis.R](Data_analysis/Silucave_mainscript_analysis.R) allows to:

- Analyse assembly statistics and GenomeScope results
- Plot species phylogeny and mitochondrial phylogenies
- Analyse vision genes/pseudogenes and associated loss-of-function mutations
- Perform datations of cave colonization using vision genes dN/dS and pseudogenes proportion
- Analyse genome-wide dN/dS values
- Analyse RERconverge results
- Analyse RELAX results
- Analyse CSUBST results (computed with both the gene trees or the species tree)
- Analyse Mutpred2 (+ grantham distances) results

All files necessary to run this script are available on Zenodo (10.5281/zenodo.17432431 ; https://zenodo.org/records/17432431) and in the folder "Data_analysis"


