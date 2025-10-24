#!/bin/bash


#SBATCH --job-name=Consensus_Exons  # Job name


eval "$(conda shell.bash hook)"
conda activate miniprot

LMOD_DISABLE_SAME_NAME_AUTOSWAP=no


gene=$1 
ref_genome=Trichomycterus_rosablanca/GCF_030014385.1_fTriRos1.hap1_genomic.fna


##Generate a consensus of exon sequences for each gene and each species, and merge exons

module purge ; module load SAMtools ; module load BCFtools

#extract the gene strand
strand=`grep "$gene" Trichomycterus_rosablanca.gff.exons | cut -f7 | head -1`



rm -f Gene_LogFiles/slurm.$gene.out
rm -f Gene_LogFiles/error.$gene.out

for ID in RHH2 CHM6 RSS1 CSV83 RUI2 CUL4 ; do 

	if grep -q ">$ID" Consensus_Gene_Sequences.Trosa/$gene.consensus.final.fa ; then

		echo "$ID already in file"

	else 

		#Extract current exons consensus in the current individual
		xargs samtools faidx $ref_genome < Trichomycterus_rosablanca.GenesInfos/$gene.pos | bcftools consensus $ID.norm.flt-indels.bcf  >> Consensus_Gene_Sequences.Trosa/$gene.$ID.exons.fa
		
		#Reformat the fasta file for samtools
		sed -i 's/:/-/g'  Consensus_Gene_Sequences.Trosa/$gene.$ID.exons.fa
		fasta_formatter -i Consensus_Gene_Sequences.Trosa/$gene.$ID.exons.fa -o Consensus_Gene_Sequences.Trosa/$gene.$ID.exons.reformat.fa -w 60 ; rm Consensus_Gene_Sequences.Trosa/$gene.$ID.exons.fa
	
	
		if [ $strand == "+" ] ; then
			xargs samtools faidx Consensus_Gene_Sequences.Trosa/$gene.$ID.exons.reformat.fa < Trichomycterus_rosablanca.GenesInfos/$gene.tiret.pos | grep -v ">" | sed "1i\>$ID"  >> Consensus_Gene_Sequences.Trosa/$gene.consensus.final.fa
		else
	
			for curr_exon in `cat Trichomycterus_rosablanca.GenesInfos/$gene.tiret.pos.rev` ; do 
	
				samtools faidx Consensus_Gene_Sequences.Trosa/$gene.$ID.exons.reformat.fa $curr_exon | grep -v ">" > Consensus_Gene_Sequences.Trosa/$gene.$ID.temp.fa
				revseq Consensus_Gene_Sequences.Trosa/$gene.$ID.temp.fa Consensus_Gene_Sequences.Trosa/$gene.$ID.temp.fa.rev
				grep -v ">" Consensus_Gene_Sequences.Trosa/$gene.$ID.temp.fa.rev >> Consensus_Gene_Sequences.Trosa/$gene.$ID.rev
	
			done
	
			sed "1i\>$ID" Consensus_Gene_Sequences.Trosa/$gene.$ID.rev >> Consensus_Gene_Sequences.Trosa/$gene.consensus.final.fa
	
			rm Consensus_Gene_Sequences.Trosa/$gene.$ID.rev
			rm Consensus_Gene_Sequences.Trosa/$gene.$ID.temp.fa.rev
			rm Consensus_Gene_Sequences.Trosa/$gene.$ID.temp.fa
		fi
	
	
		rm Consensus_Gene_Sequences.Trosa/$gene.$ID.exons.reformat.fa*


	fi

done



