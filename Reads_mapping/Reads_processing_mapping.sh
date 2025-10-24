###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###=============================== Clean reads  ======================================================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


module purge ; module load SRA-Toolkit ; module load fastp

fastp -i RHP1_AACGTTCC-AGTACTCC-BHV5HYDSXY_L001_R1.fastq.gz -I RHP1_AACGTTCC-AGTACTCC-BHV5HYDSXY_L001_R2.fastq.gz -o RHP1.R1.cleaned.fastq.gz -O RHP1.R2.cleaned.fastq.gz -h RHP1.html -j RHP1.json
fastp -i CUL9_GCAGAATT-TGGCCGGT-BHV5HYDSXY_L001_R1.fastq.gz -I CUL9_GCAGAATT-TGGCCGGT-BHV5HYDSXY_L001_R2.fastq.gz -o CUL9.R1.cleaned.fastq.gz -O CUL9.R2.cleaned.fastq.gz -h CUL9.html -j CUL9.json
fastp -i CHM6_TAATACAG-GTGAATAT-BHV5HYDSXY_L001_R1.fastq.gz  -I CHM6_TAATACAG-GTGAATAT-BHV5HYDSXY_L001_R2.fastq.gz -o CHM6.R1.cleaned.fastq.gz -O CHM6.R2.cleaned.fastq.gz -h CHM6.html -j CHM6.json
fastp -i CSV83_ATGTAAGT-CATAGAGT-BHV5HYDSXY_L001_R1.fastq.gz -I CSV83_ATGTAAGT-CATAGAGT-BHV5HYDSXY_L001_R2.fastq.gz -o CSV83.R1.cleaned.fastq.gz -O CSV83.R2.cleaned.fastq.gz -h CSV83.html -j CSV83.json
fastp -i RHH2_AATCCGGA-AACTGTAG-BHV5HYDSXY_L001_R1.fastq.gz  -I RHH2_AATCCGGA-AACTGTAG-BHV5HYDSXY_L001_R2.fastq.gz  -o RHH2.R1.cleaned.fastq.gz -O RHH2.R2.cleaned.fastq.gz -h RHH2.html -j RHH2.json
fastp -i RSS1_CGGCGTGA-ACAGGCGC-BHV5HYDSXY_L001_R1.fastq.gz  -I RSS1_CGGCGTGA-ACAGGCGC-BHV5HYDSXY_L001_R2.fastq.gz  -o RSS1.R1.cleaned.fastq.gz -O RSS1.R2.cleaned.fastq.gz -h RSS1.html -j RSS1.json
fastp -i RUI2_GCACGGAC-TGCGAGAC-BHV5HYDSXY_L001_R1.fastq.gz  -I RUI2_GCACGGAC-TGCGAGAC-BHV5HYDSXY_L001_R2.fastq.gz  -o RUI2.R1.cleaned.fastq.gz -O RUI2.R2.cleaned.fastq.gz -h RUI2.html -j RUI2.json
fastp -i CUL4_GGTACCTT-GACGTCTT-BHV5HYDSXY_L001_R1.fastq.gz  -I CUL4_GGTACCTT-GACGTCTT-BHV5HYDSXY_L001_R2.fastq.gz  -o CUL4.R1.cleaned.fastq.gz -O CUL4.R2.cleaned.fastq.gz -h CUL4.html -j CUL4.json


###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###=============================== Map reads to references genomes  =================================================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================



#First index genomes
bwa-mem2 index GCA_025084515.1_ASM2508451v1_genomic.fna #Ancistrus_triradiatus genome
bwa-mem2 index GCF_001660625.3_Coco_2.0_genomic.fna #Ictalurus_punctatus genome
bwa-mem2 index GCF_030014385.1_fTriRos1.hap1_genomic.fna #Trichomycterus_rosablanca genome


#Now lets map to the closest species genomes
for line in `cat Table_mapping.tsv` ; do 
	ref_genome=`echo "$line" | cut -f2 -d ","`
	sample_name=`echo "$line" | cut -f1 -d ","`

	sbatch --qos=1day -c 40 --mem=50G -e error.$sample_name.out -o slurm.$sample_name.out launch_mapping.sh $ref_genome $sample_name
done


###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###=============================== Call variants  ==================================================================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


mkdir Log_BCF

for line in `cat Table_mapping.tsv` ; do 
	ref_genome=`echo "$line" | cut -f2 -d ","`
	sample_name=`echo "$line" | cut -f1 -d ","`

	sbatch --qos=1day -c 20 --mem=50G -e Log_BCF/error.$sample_name.out -o Log_BCF/slurm.$sample_name.out call_variants.sh $ref_genome $sample_name
done


###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###=============================== Generate consensus genomes and run BUSCO  ========================================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


mkdir Log_consensus

for line in `cat Table_mapping.tsv` ; do 
	ref_genome=`echo "$line" | cut -f2 -d ","`
	sample_name=`echo "$line" | cut -f1 -d ","`

	sbatch --qos=1day -c 20 --mem=50G -e Log_consensus/error.$sample_name.out -o Log_consensus/slurm.$sample_name.out Generate_consensus_genome.sh $ref_genome $sample_name
done


###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###=============================== Download all Siluriformes genomes  ===============================================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================



git clone https://github.com/pirovc/genome_updater.git
./genome_updater/genome_updater.sh -d "refseq,genbank" -T "7995" -f "assembly_report.txt" -t 1 -m -A 1 -o Siluriformes_assembly_reports -L curl #7995 is the taxonomic ID of Siluriformes
cp Siluriformes_assembly_reports/assembly_summary.txt ./assembly_summary.txt

#Remove partial assemblies 
awk -v FS="\t" '$14=="Full"' assembly_summary.txt > full_assembly_summary.txt

## Now lets download all the genomes, and the annotation (if it exists)

module purge ; module load SAMtools

IFS=$'\n'
for line in `cat full_assembly_summary.txt` ; do

	species=`echo "$line" | cut -f8 | sed "s/'//g" | sed 's/\.//g' | sed 's/ /_/g'`
	
	echo "$line" > line.txt
	
	rm -r $species ; mkdir $species
	
	if grep -q "GCF_" line.txt ; then 
	
		GCA_name=`echo "$line" | cut -f1`
		GCF_name=`echo "$line" | cut -f18`
		assembly_link=`echo "$line" | cut -f20 | sed 's/$/\//g' | sed "s/$GCA_name/$GCF_name/g" | sed 's/\/GCA\//\/GCF\//g'`
		assembly_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.fna.gz/g' | sed "s/$GCA_name/$GCF_name/g"`
		GFF_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.gff.gz/g' | sed "s/$GCA_name/$GCF_name/g"`
		report_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_assembly_report.txt/g' | sed "s/$GCA_name/$GCF_name/g"`
		
		report_full_link=`echo "$assembly_link$report_name"`
		full_fasta_link=`echo "$assembly_link$assembly_name"`
		full_gff_link=`echo "$assembly_link$GFF_name"`
		
		
		cd $species
		
		wget $full_gff_link ; wget $full_fasta_link ; wget $report_full_link
		gzip -d $assembly_name ; gzip -d $GFF_name 



		if grep -q "alt-scaffold" $report_name ; then

			grep -v "alt-scaffold" $report_name | cut -f7 | grep -v "^na$" | grep -v "#" | grep -v "RefSeq-Accn" | sed '/^[[:space:]]*$/d' > good_scaffold_list.txt
			grep  "alt-scaffold" $report_name | grep -v "^na$" > infos_removed_genome_parts.txt

			trimmed_genome_name=`echo "$assembly_name" | sed 's/.fna.gz/.trimmed.fna/g'`
			genome_name=`echo "$assembly_name" | sed 's/.fna.gz/.fna/g'`
			xargs samtools faidx $genome_name < good_scaffold_list.txt > $trimmed_genome_name
			gzip $genome_name


			echo "$genome_name" >> ../List_of_modified_genomes.txt

		fi


		
		cd ../
	
	
	else
	
		assembly_link=`echo "$line" | cut -f20 | sed 's/$/\//g'`
		assembly_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.fna.gz/g'`
		GFF_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.gff.gz/g'`
		report_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_assembly_report.txt/g'`

		
		full_fasta_link=`echo "$assembly_link$assembly_name"`
		full_gff_link=`echo "$assembly_link$GFF_name"`
		report_full_link=`echo "$assembly_link$report_name"`

		
		cd $species
		

		wget $full_gff_link ; wget $full_fasta_link ; wget $report_full_link
		gzip -d $assembly_name ; gzip -d $GFF_name 
		

		if grep -q "alt-scaffold" $report_name ; then

	
			grep -v "alt-scaffold" $report_name | cut -f5 | grep -v "^na$" | grep -v "#" | grep -v "GenBank-Accn" | sed '/^[[:space:]]*$/d' > good_scaffold_list.txt
			grep  "alt-scaffold" $report_name | grep -v "^na$" > infos_removed_genome_parts.txt

			trimmed_genome_name=`echo "$assembly_name" | sed 's/.fna.gz/.trimmed.fna/g'`
			genome_name=`echo "$assembly_name" | sed 's/.fna.gz/.fna/g'`
			xargs samtools faidx $genome_name < good_scaffold_list.txt > $trimmed_genome_name
			gzip $genome_name

			echo "$genome_name" >> ../List_of_modified_genomes.txt



		fi


		cd ../
	
	fi


	reduced_assembly_name=`echo "$assembly_name" | sed 's/_genomic.fna.gz//g'`
	if ls -l $species/ | grep -q ".gff" ; then 
		echo "$species,$reduced_assembly_name" >> Annotated_genome_list.txt
	else 
		echo "$species,$reduced_assembly_name" >> Non_annotated_genome_list.txt
	fi 


done


#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Run BUSCO on genome assemblises ##################################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


mkdir BUSCO_logs

IFS=$'\n'
for line in `cat full_assembly_summary.txt` ; do

	species=`echo "$line" | cut -f8 | sed "s/'//g" | sed 's/\.//g' | sed 's/ /_/g'`
	
	echo "$line" > curr_line.txt
		
	GCA_name=`echo "$line" | cut -f1`
	GCF_name=`echo "$line" | cut -f18`

	if grep -q "GCF_" curr_line.txt ; then 

		assembly_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.fna/g' | sed "s/$GCA_name/$GCF_name/g"`
		if grep -q "$assembly_name" List_of_modified_genomes.txt ; then assembly_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.trimmed.fna/g' | sed "s/$GCA_name/$GCF_name/g"` ; fi
	
	
	else
	
		assembly_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.fna/g'`
		if grep -q "$assembly_name" List_of_modified_genomes.txt ; then assembly_name=`echo "$line" | cut -f20 | sed 's/.*\///g' | sed 's/$/_genomic.trimmed.fna/g'` ; fi

		
	fi

	sbatch --qos=6hours -c 40 --mem=50G -e BUSCO_logs/error_$species.out -o BUSCO_logs/slurm_$species.out launch_busco.sh $species $assembly_name
done



#Extract BUSCO results


IFS=$'\n'
for line in `cat full_assembly_summary.txt` ; do

	species=`echo "$line" | cut -f8 | sed "s/'//g" | sed 's/\.//g' | sed 's/ /_/g'`


	Complete=`grep "Complete BUSCOs" BUSCO_$species/short_summary.specific.actinopterygii_odb10.BUSCO_$species.txt | sed 's/^	*//g' | sed 's/	.*//g'`
	Complete_single=`grep "Complete and single-copy" BUSCO_$species/short_summary.specific.actinopterygii_odb10.BUSCO_$species.txt | sed 's/^	*//g' | sed 's/	.*//g'`
	Complete_dup=`grep "Complete and duplicated" BUSCO_$species/short_summary.specific.actinopterygii_odb10.BUSCO_$species.txt | sed 's/^	*//g' | sed 's/	.*//g'`
	Fragmented=`grep "Fragmented" BUSCO_$species/short_summary.specific.actinopterygii_odb10.BUSCO_$species.txt | sed 's/^	*//g' | sed 's/	.*//g'`
	Missing=`grep "Missing" BUSCO_$species/short_summary.specific.actinopterygii_odb10.BUSCO_$species.txt | sed 's/^	*//g' | sed 's/	.*//g'`
	Total=`grep "Total BUSCO groups" BUSCO_$species/short_summary.specific.actinopterygii_odb10.BUSCO_$species.txt | sed 's/^	*//g' | sed 's/	.*//g'`

	Scaffold_nb=`grep "Number of scaffolds" BUSCO_$species/short_summary.specific.actinopterygii_odb10.BUSCO_$species.txt | sed 's/^	*//g' | sed 's/	.*//g'`
	Contig_nb=`grep "Number of contigs" BUSCO_$species/short_summary.specific.actinopterygii_odb10.BUSCO_$species.txt | sed 's/^	*//g' | sed 's/	.*//g'`
	Genome_size=`grep "Total length" BUSCO_$species/short_summary.specific.actinopterygii_odb10.BUSCO_$species.txt | sed 's/^	*//g' | sed 's/	.*//g'`
	Gaps_perc=`grep "Percent gaps" BUSCO_$species/short_summary.specific.actinopterygii_odb10.BUSCO_$species.txt | sed 's/^	*//g' | sed 's/	.*//g' | sed 's/%//g'`
	SN50=`grep "Scaffold N50" BUSCO_$species/short_summary.specific.actinopterygii_odb10.BUSCO_$species.txt | sed 's/^	*//g' | sed 's/	.*//g' | sed 's/ MB/000000/g' | sed 's/ KB/000/g'`
	CN50=`grep "Contigs N50" BUSCO_$species/short_summary.specific.actinopterygii_odb10.BUSCO_$species.txt | sed 's/^	*//g' | sed 's/	.*//g' | sed 's/ MB/000000/g' | sed 's/ KB/000/g'`


	echo "$species,$Complete,$Complete_single,$Complete_dup,$Fragmented,$Missing,$Total,$Scaffold_nb,$Contig_nb,$Genome_size,$Gaps_perc,$SN50,$CN50"  


done > BUSCO.results.csv




# Do the same for consensus genomes generated from reads mapping

for assembly in RHP1 CUL9 CHM6 CSV83 RHH2 RSS1 RUI2 CUL4 ; do
	assembly_file=$assembly.consensus_genome.fna
	sbatch --qos=6hours -c 40 --mem=50G launch_busco.new.sh $assembly_file $assembly
done


for assembly in RHP1 CUL9 CHM6 CSV83 RHH2 RSS1 RUI2 CUL4 ; do

	Complete=`grep "Complete BUSCOs" BUSCO_$assembly/short_summary.specific.actinopterygii_odb10.BUSCO_$assembly.txt | sed 's/^	*//g' | sed 's/	.*//g'`
	Complete_single=`grep "Complete and single-copy" BUSCO_$assembly/short_summary.specific.actinopterygii_odb10.BUSCO_$assembly.txt | sed 's/^	*//g' | sed 's/	.*//g'`
	Complete_dup=`grep "Complete and duplicated" BUSCO_$assembly/short_summary.specific.actinopterygii_odb10.BUSCO_$assembly.txt | sed 's/^	*//g' | sed 's/	.*//g'`
	Fragmented=`grep "Fragmented" BUSCO_$assembly/short_summary.specific.actinopterygii_odb10.BUSCO_$assembly.txt | sed 's/^	*//g' | sed 's/	.*//g'`
	Missing=`grep "Missing" BUSCO_$assembly/short_summary.specific.actinopterygii_odb10.BUSCO_$assembly.txt | sed 's/^	*//g' | sed 's/	.*//g'`
	Total=`grep "Total BUSCO groups" BUSCO_$assembly/short_summary.specific.actinopterygii_odb10.BUSCO_$assembly.txt | sed 's/^	*//g' | sed 's/	.*//g'`

	echo "$assembly,$Complete,$Complete_single,$Complete_dup,$Fragmented,$Missing,$Total" | tr ',' '\t'
done




##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################# Now Lets Launch BRAKER3 for non-annotated species ##################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

cut -f1 -d "," Non_annotated_genome_list.txt > species_to_annotate

echo "Ancistrus_triradiatus" > species_to_annotate
echo "Channallabes_apus" >> species_to_annotate
echo "Corydoras_maculifer" >> species_to_annotate
echo "Cranoglanis_bouderius" >> species_to_annotate
echo "Doras_micropoeus" >> species_to_annotate
echo "Duopalatinus_emarginatus" >> species_to_annotate
echo "Gogo_arcuatus" >> species_to_annotate
echo "Microcambeva_barbata" >> species_to_annotate
echo "Notoglanidium_macrostoma" >> species_to_annotate
echo "Plotosus_lineatus" >> species_to_annotate
echo "Prietella_phreatophila" >> species_to_annotate
echo "Scorpiodoras_heckelii" >> species_to_annotate
echo "Synodontis_membranacea" >> species_to_annotate



#SRA_data_filt.txt:
#Ancistrus_triradiatus,SRR14702161
#Corydoras_maculifer,SRR13849942
#Plotosus_lineatus,SRR20662526,SRR20662527,SRR18748376,SRR18748390,SRR18748391
#Cranoglanis_bouderius,SRR19216508,SRR19216380,SRR19216030,SRR19215954,SRR19214211


#Download SRA data for species for which there are some available


module load SRA-Toolkit

./download_SRA.single.sh Ancistrus_triradiatus
./download_SRA.sh Corydoras_maculifer
./download_SRA.sh Plotosus_lineatus
./download_SRA.sh Cranoglanis_bouderius


mv Cranoglanis_bouderius_R1.fastq Cranoglanis_bouderius_1.fastq
mv Cranoglanis_bouderius.R2.fastq Cranoglanis_bouderius_2.fastq

mv Plotosus_lineatus.R1.fastq Plotosus_lineatus_1.fastq
mv Plotosus_lineatus.R2.fastq Plotosus_lineatus_2.fastq


#Launch BRAKER

for species in `cat species_to_annotate` ; do 

	genome_file=`ls -l $species | grep ".fna$" | sed 's/.* //g'`
	masked_Genome=`echo "$genome_file"`
	working_dir=`echo "$species.BRAKER"`
	rm -r $working_dir ; mkdir $working_dir
	sed 's/ .*//g'$masked_Genome > $masked_Genome.renamed

	if grep -q "$species" SRA_data_species.txt ; then 

		sbatch --qos=1week -c 30 --mem=200G --job-name=RNA.$species launch_Braker_wSRA_ubuntu.sh /scicore/home/salzburg/polica0000/SiluCave/Genomes/$species/$masked_Genome.renamed /scicore/home/salzburg/polica0000/Horizontal_transfer_project/Homemade_annotation/BRAKER_Results/Custom_ProteinDB.renamed.fa $species $working_dir

	else
		
		sbatch --qos=1week -c 25 --mem=200G --job-name=noRNA.$species launch_Braker_woSRA_ubuntu.sh /scicore/home/salzburg/polica0000/SiluCave/Genomes/$species/$masked_Genome.renamed /scicore/home/salzburg/polica0000/Horizontal_transfer_project/Homemade_annotation/BRAKER_Results/Custom_ProteinDB.renamed.fa $species $working_dir

	fi

done


#Count the number of coding sequences annotated per species (verification that BRAKER worked)


for species in `cat species_to_annotate` ; do echo "$species" ; grep -c ">" $species.BRAKER/braker.codingseq ; done


#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Extract coding sequences from GFF3 files ######################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

#Extract coding sequences from GFF3 annotations download from RefSeq or Genbank

for line in `cat Annotated_genome_list.txt` ; do 
	species=`echo "$line" | cut -f1 -d ","`
	sbatch --qos=6hours -c 8 --mem=30G agat_process.sh $species
done



#Extract coding sequences from GFF3 annotations generated with BRAKER


for line in `cat species_to_annotate` ; do 
	species=`echo "$line" | cut -f1 -d ","`
	sbatch --qos=6hours -c 8 --mem=30G Parse_braker_results.sh $species
done


for species in `cat species_to_annotate` ; do 
	sed -i "s/>/>$species---/g" $species.BRAKER/$species.codingseq.LI.nostop
	sed -i "s/>/>$species---/g" $species.BRAKER/$species.prot.LI.nostop

done


#Put all coding sequences and protein fasta file into the same folders


rm -rf Coding_sequences ; mkdir Coding_sequences
rm -rf Proteins ; mkdir Proteins

IFS=$'\n'
for line in `cat Annotated_genome_list.txt` ; do 
	species=`echo "$line" | cut -f1 -d ","`
	grep ">" $species/$species.prot.nostop | sed 's/>//g' > tokeep.id
	xargs samtools faidx $species/$species.renamed.cds < tokeep.id  > Coding_sequences/$species.cds
	cat $species/$species.prot.nostop > Proteins/$species.prot
done
#Do the same steps above for the zebrafish .. 
cat Danio_rerio/Danio_rerio.prot.nostop > Proteins/Danio_rerio.prot
cat Danio_rerio/Danio_rerio.cds.nostop > Coding_sequences/Danio_rerio.cds


for species in `cat species_to_annotate` ; do 
	cat $species.BRAKER/$species.codingseq.LI.nostop > Coding_sequences/$species.cds
	cat $species.BRAKER/$species.prot.LI.nostop > Proteins/$species.prot
done


###### Now we do a final filter step to remove transposable elements.

ls -l Coding_sequences/ | sed 's/.* //g' | grep "cds" | sed 's/.cds//g'  > all_sp.txt

mkdir Coding_sequences_woTE
mkdir Proteins_woTE
for species in `cat all_sp.txt`  ; do 
	sbatch --qos=6hours -c 10 --mem=20G Remove_TE_annotations.sh $species
done



#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Extract coding sequences from consensus genomes generated from reads mapping on Trichomycterus_rosablanca  ########
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

###### In a final step, extract coding sequences from mapped genomes (CUL9, etc .. )


grep ">" Coding_sequences_woTE/Trichomycterus_rosablanca.cds  | sed 's/>//g' #Extract coding sequence from the reference genome of Trichomycterus_rosablanca


#Prepare exons of T. rosablanca

mkdir Trichomycterus_rosablanca.GenesInfos

grep "	CDS	" Trichomycterus_rosablanca/GCF_030014385.1_fTriRos1.hap1_genomic.gff.LI > Trichomycterus_rosablanca.gff.exons
cut -f1 Trichomycterus_rosablanca.gff.exons > Trichomycterus_rosablanca.exons_scaffold_name.txt
cut -f4,5 Trichomycterus_rosablanca.gff.exons | tr "\t" "-" > Trichomycterus_rosablanca.exons_locations.txt
paste -d ":" Trichomycterus_rosablanca.exons_scaffold_name.txt Trichomycterus_rosablanca.exons_locations.txt > Trichomycterus_rosablanca.exons_scaffolds_locations.txt
cut -f9 Trichomycterus_rosablanca.gff.exons | sed 's/.*Parent=rna-//g' | sed 's/.*Parent=id-//g' | sed 's/;.*//g' > Trichomycterus_rosablanca.exons_XM_id.txt
cut -f9 Trichomycterus_rosablanca.gff.exons | sed 's/.*ID=cds-//g' | sed 's/;.*//g' > Trichomycterus_rosablanca.exons_XP_id.txt
cut -f7 Trichomycterus_rosablanca.gff.exons > Trichomycterus_rosablanca.strand_exons

paste -d "," Trichomycterus_rosablanca.exons_scaffolds_locations.txt Trichomycterus_rosablanca.exons_XM_id.txt Trichomycterus_rosablanca.exons_XP_id.txt > Trichomycterus_rosablanca.exons_positions_infos.csv
paste -d "," Trichomycterus_rosablanca.exons_positions_infos.csv Trichomycterus_rosablanca.strand_exons > temp ; mv temp Trichomycterus_rosablanca.exons_positions_infos.csv
awk -F',' 'BEGIN {OFS=","} {print $0, NR}' Trichomycterus_rosablanca.exons_positions_infos.csv > temp ; mv temp Trichomycterus_rosablanca.exons_positions_infos.csv

cut -f2,3 -d "," Trichomycterus_rosablanca.exons_positions_infos.csv | sed 's/,/-/g' | sort | uniq > Trichomycterus_rosablanca.XP-XM.txt
cut -f2 -d "," Trichomycterus_rosablanca.exons_positions_infos.csv | sort | uniq > Trichomycterus_rosablanca.XM.txt

cut -f2,4,5 -d "," Trichomycterus_rosablanca.exons_positions_infos.csv | sed 's/,/_/g' > Trichomycterus_rosablanca.new_names.txt
paste -d "\t" Trichomycterus_rosablanca.exons_scaffolds_locations.txt Trichomycterus_rosablanca.new_names.txt  > Trichomycterus_rosablanca.names_mapping.tsv
sed -i 's/:/-/g' Trichomycterus_rosablanca.names_mapping.tsv

grep ">" Coding_sequences_woTE/Trichomycterus_rosablanca.cds  | sed 's/.*rna-//g' | sed 's/.*LOC/LOC/g' >  Trichomycterus_rosablanca.XM.noTE.txt


for gene in `cat Trichomycterus_rosablanca.XM.noTE.txt` ; do
	grep ",$gene," Trichomycterus_rosablanca.exons_positions_infos.csv | cut -f1 -d "," > Trichomycterus_rosablanca.GenesInfos/$gene.pos
	sed 's/:/-/g' Trichomycterus_rosablanca.GenesInfos/$gene.pos > Trichomycterus_rosablanca.GenesInfos/$gene.tiret.pos
	tac Trichomycterus_rosablanca.GenesInfos/$gene.tiret.pos > Trichomycterus_rosablanca.GenesInfos/$gene.tiret.pos.rev
done 


rm -rf Gene_LogFiles ; mkdir Gene_LogFiles
rm -rf Consensus_Gene_Sequences.Trosa ; mkdir Consensus_Gene_Sequences.Trosa

sbatch launcharray.trosa.sh



====== launcharray.trosa.sh

#!/bin/bash

#SBATCH --job-name=Consensus_Exons 
#SBATCH --qos=6hours
#SBATCH -c 2
#SBATCH --mem=5G
#SBATCH --error=Gene_LogFiles/error.%A_%a.out
#SBATCH --output=Gene_LogFiles/slurm.%A_%a.out
#SBATCH --array=1-22874%1000

#Remove existing log files
rm -f Gene_LogFiles/error.*.out
rm -f Gene_LogFiles/slurm.*.out

# Extract the gene corresponding to the array task ID
gene=$(sed -n "${SLURM_ARRAY_TASK_ID}p" Trichomycterus_rosablanca.XM.noTE.txt)

# Call the script to process each gene
./Generate_consensus_exons.Trosa.sh $gene


======================================================


for gene in `cat Trichomycterus_rosablanca.XM.noTE.txt` ; do x=`grep -c ">" Consensus_Gene_Sequences.Trosa/$gene.consensus.final.fa` ; echo "$gene,$x" ; done > list_fin


for gene in `cat Trichomycterus_rosablanca.XM.noTE.txt` ; do 
	fasta_formatter -i Consensus_Gene_Sequences.Trosa/$gene.consensus.final.fa -o Consensus_Gene_Sequences.Trosa/$gene.consensus.final.reformat.fa -w 60 
	mv Consensus_Gene_Sequences.Trosa/$gene.consensus.final.reformat.fa Consensus_Gene_Sequences.Trosa/$gene.consensus.final.fa

	samtools faidx Consensus_Gene_Sequences.Trosa/$gene.consensus.final.fa RHH2 | sed "s/>.*/>RHH2---$gene/g">> Coding_sequences_woTE/RHH2.cds
	samtools faidx Consensus_Gene_Sequences.Trosa/$gene.consensus.final.fa CHM6 | sed "s/>.*/>CHM6---$gene/g">> Coding_sequences_woTE/CHM6.cds
	samtools faidx Consensus_Gene_Sequences.Trosa/$gene.consensus.final.fa RSS1 | sed "s/>.*/>RSS1---$gene/g">> Coding_sequences_woTE/RSS1.cds
	samtools faidx Consensus_Gene_Sequences.Trosa/$gene.consensus.final.fa CSV83 | sed "s/>.*/>CSV83---$gene/g">> Coding_sequences_woTE/CSV83.cds
	samtools faidx Consensus_Gene_Sequences.Trosa/$gene.consensus.final.fa RUI2 | sed "s/>.*/>RUI2---$gene/g">> Coding_sequences_woTE/RUI2.cds
	samtools faidx Consensus_Gene_Sequences.Trosa/$gene.consensus.final.fa CUL4 | sed "s/>.*/>CUL4---$gene/g">> Coding_sequences_woTE/CUL4.cds

	rm Consensus_Gene_Sequences.Trosa/$gene.consensus.final.fa.fai
done


for curr_ID in RHH2 CHM6 RSS1 CSV83 RUI2 CUL4 ; do 

	transeq Coding_sequences_woTE/$curr_ID.cds Proteins_woTE/$curr_ID.prot ; sed -i 's/_1$//g' Proteins_woTE/$curr_ID.prot

	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Proteins_woTE/$curr_ID.prot | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\*/)' | tr "\t" "\n" > Proteins_woTE/$curr_ID.prot.nostop
	grep ">" Proteins_woTE/$curr_ID.prot.nostop | sed 's/>//g' | sort | uniq > good_id.txt
	xargs samtools faidx Coding_sequences_woTE/$curr_ID.cds < good_id.txt > Coding_sequences_woTE/$curr_ID.cds.nostop ; rm good_id.txt
	mv Coding_sequences_woTE/$curr_ID.cds.nostop Coding_sequences_woTE/$curr_ID.cds ; mv Proteins_woTE/$curr_ID.prot.nostop Proteins_woTE/$curr_ID.prot

done

rm Coding_sequences_woTE/*.fai ; rm Proteins_woTE/*.fai




#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Extract coding sequences from consensus genomes generated from reads mapping on Ancistrus_triradiatus  ############
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================



#Do the same as above

mkdir Ancistrus_triradiatus.GenesInfos

grep "	CDS	" Ancistrus_triradiatus.BRAKER/braker.gff3 > Ancistrus_triradiatus.gff.exons
cut -f1 Ancistrus_triradiatus.gff.exons > Ancistrus_triradiatus.exons_scaffold_name.txt
cut -f4,5 Ancistrus_triradiatus.gff.exons | tr "\t" "-" > Ancistrus_triradiatus.exons_locations.txt
paste -d ":" Ancistrus_triradiatus.exons_scaffold_name.txt Ancistrus_triradiatus.exons_locations.txt > Ancistrus_triradiatus.exons_scaffolds_locations.txt
cut -f9 Ancistrus_triradiatus.gff.exons | sed 's/.*Parent=//g' | sed 's/;.*//g' > Ancistrus_triradiatus.exons_XM_id.txt
cut -f9 Ancistrus_triradiatus.gff.exons | sed 's/.*Parent=//g' | sed 's/;.*//g' > Ancistrus_triradiatus.exons_XP_id.txt

cut -f7 Ancistrus_triradiatus.gff.exons > Ancistrus_triradiatus.strand_exons

paste -d "," Ancistrus_triradiatus.exons_scaffolds_locations.txt Ancistrus_triradiatus.exons_XM_id.txt Ancistrus_triradiatus.exons_XP_id.txt > Ancistrus_triradiatus.exons_positions_infos.csv
paste -d "," Ancistrus_triradiatus.exons_positions_infos.csv Ancistrus_triradiatus.strand_exons > temp ; mv temp Ancistrus_triradiatus.exons_positions_infos.csv
awk -F',' 'BEGIN {OFS=","} {print $0, NR}' Ancistrus_triradiatus.exons_positions_infos.csv > temp ; mv temp Ancistrus_triradiatus.exons_positions_infos.csv

cut -f2,3 -d "," Ancistrus_triradiatus.exons_positions_infos.csv | sed 's/,/-/g' | sort | uniq > Ancistrus_triradiatus.XP-XM.txt
cut -f2 -d "," Ancistrus_triradiatus.exons_positions_infos.csv | sort | uniq > Ancistrus_triradiatus.XM.txt

cut -f2,4,5 -d "," Ancistrus_triradiatus.exons_positions_infos.csv | sed 's/,/_/g' > Ancistrus_triradiatus.new_names.txt
paste -d "\t" Ancistrus_triradiatus.exons_scaffolds_locations.txt Ancistrus_triradiatus.new_names.txt  > Ancistrus_triradiatus.names_mapping.tsv
sed -i 's/:/-/g' Ancistrus_triradiatus.names_mapping.tsv

grep ">" Coding_sequences_woTE/Ancistrus_triradiatus.cds | sed 's/.*---//g' >  Ancistrus_triradiatus.XM.noTE.txt


for gene in `cat Ancistrus_triradiatus.XM.noTE.txt` ; do
	grep ",$gene," Ancistrus_triradiatus.exons_positions_infos.csv | cut -f1 -d "," > Ancistrus_triradiatus.GenesInfos/$gene.pos
	sed 's/:/-/g' Ancistrus_triradiatus.GenesInfos/$gene.pos > Ancistrus_triradiatus.GenesInfos/$gene.tiret.pos
	tac Ancistrus_triradiatus.GenesInfos/$gene.tiret.pos > Ancistrus_triradiatus.GenesInfos/$gene.tiret.pos.rev
done 


rm -rf Gene_LogFiles ; mkdir Gene_LogFiles
rm -rf Consensus_Gene_Sequences.Atri ; mkdir Consensus_Gene_Sequences.Atri

sbatch launcharray.atri.sh



====== launcharray.atri.sh

#!/bin/bash

#SBATCH --job-name=Consensus_Exons 
#SBATCH --qos=6hours
#SBATCH -c 2
#SBATCH --mem=5G
#SBATCH --error=Gene_LogFiles/error.%A_%a.out
#SBATCH --output=Gene_LogFiles/slurm.%A_%a.out
#SBATCH --array=1-17417%1000

#Remove existing log files
rm -f Gene_LogFiles/error.*.out
rm -f Gene_LogFiles/slurm.*.out

# Extract the gene corresponding to the array task ID
gene=$(sed -n "${SLURM_ARRAY_TASK_ID}p" Ancistrus_triradiatus.XM.noTE.txt)

# Call the script to process each gene
./Generate_consensus_exons.Atri.sh $gene

======================================================


for gene in `cat Ancistrus_triradiatus.XM.noTE.txt` ; do grep -c ">" Consensus_Gene_Sequences.Atri/$gene.consensus.final.fa ; done > list_fin


for gene in `cat Ancistrus_triradiatus.XM.noTE.txt` ; do 

	echo "$gene" 
	fasta_formatter -i Consensus_Gene_Sequences.Atri/$gene.consensus.final.fa -o Consensus_Gene_Sequences.Atri/$gene.consensus.final.reformat.fa -w 60 
	mv Consensus_Gene_Sequences.Atri/$gene.consensus.final.reformat.fa Consensus_Gene_Sequences.Atri/$gene.consensus.final.fa

	samtools faidx Consensus_Gene_Sequences.Atri/$gene.consensus.final.fa RHP1 | sed "s/>.*/>RHP1---$gene/g">> Coding_sequences_woTE/RHP1.cds
	samtools faidx Consensus_Gene_Sequences.Atri/$gene.consensus.final.fa CUL9 | sed "s/>.*/>CUL9---$gene/g">> Coding_sequences_woTE/CUL9.cds

	rm Consensus_Gene_Sequences.Atri/$gene.consensus.final.fa.fai
done

transeq Coding_sequences_woTE/RHP1.cds Proteins_woTE/RHP1.prot ; sed -i 's/_1$//g' Proteins_woTE/RHP1.prot
transeq Coding_sequences_woTE/CUL9.cds Proteins_woTE/CUL9.prot ; sed -i 's/_1$//g' Proteins_woTE/CUL9.prot


awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Proteins_woTE/RHP1.prot | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\*/)' | tr "\t" "\n" > Proteins_woTE/RHP1.prot.nostop
grep ">" Proteins_woTE/RHP1.prot.nostop | sed 's/>//g' | sort | uniq > good_id.txt
xargs samtools faidx Coding_sequences_woTE/RHP1.cds < good_id.txt > Coding_sequences_woTE/RHP1.cds.nostop ; rm good_id.txt
mv Coding_sequences_woTE/RHP1.cds.nostop Coding_sequences_woTE/RHP1.cds ; mv Proteins_woTE/RHP1.prot.nostop Proteins_woTE/RHP1.prot
rm Coding_sequences_woTE/*.fai ; rm Proteins_woTE/*.fai


awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Proteins_woTE/CUL9.prot | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\*/)' | tr "\t" "\n" > Proteins_woTE/CUL9.prot.nostop
grep ">" Proteins_woTE/CUL9.prot.nostop | sed 's/>//g' | sort | uniq > good_id.txt
xargs samtools faidx Coding_sequences_woTE/CUL9.cds < good_id.txt > Coding_sequences_woTE/CUL9.cds.nostop ; rm good_id.txt
mv Coding_sequences_woTE/CUL9.cds.nostop Coding_sequences_woTE/CUL9.cds ; mv Proteins_woTE/CUL9.prot.nostop Proteins_woTE/CUL9.prot
rm Coding_sequences_woTE/*.fai ; rm Proteins_woTE/*.fai



#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Extract genome assembly sizes #####################################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


for species in `cat all_sp.txt` ; do 
	genome_assembly=`ls -l /$species/ | grep ".fna$" | sed 's/.* //g'`
	size=`cut -f2 $species/$genome_assembly.fai | awk '{sum+=$1;} END{print sum;}'`
	echo "$species,$size" >> assembly_sizes.csv
done


#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## JellyFish Kmer count #############################################################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================

conda create -n GenomeScope
conda activate GenomeScope
conda install bioconda::kmc


ls -l | grep ".cleaned.fastq" | sed 's/.* //g' | grep -v "SRR" | sed 's/\..*//g' | sort | uniq > list_sp.txt


#Run genomescope on our species reads
for species_ID in `cat list_sp.txt` ; do sbatch --qos=6hours -c 10 --mem=64G kmer_genomescope.sh $species_ID 21 ; done


### Now download reads the rest of the species found in genbank or refseq

module purge ; module load SRA-Toolkit

IFS=$'\n'
for line in `cat SRA_DNA_reads.tsv` ; do 

	if echo "$line" | grep -q "No reads available" ; then
		echo "No reads available"
	else 

		SRA_accession=`echo "$line" | cut -f2`
		species=`echo "$line" | cut -f1`

		fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip $SRA_accession

		mv $SRA_accession\_pass_1.fastq.gz $species.R1.fastq.gz
		mv $SRA_accession\_pass_2.fastq.gz $species.R2.fastq.gz
done


grep -v "No reads" SRA_DNA_reads.tsv | cut -f1 > species_SRA.txt


for species in `cat species_SRA.txt` ; do ls -lh | grep "$species.*.fastq.gz" | wc -l ; done

#Clean reads
for species in `cat species_SRA.txt` ; do sbatch --qos=1day -c 8 --mem=50G clean_reads_fastp.sh $species ; done


#Run genomescope on cleaned reads
for species_ID in `cat species_SRA.txt` ; do sbatch --qos=6hours -c 10 --mem=64G kmer_genomescope.SRA.sh $species_ID ; done



sbatch --qos=6hours -c 10 --mem=64G kmer_genomescope.SRA.sh Microcambeva_barbata
sbatch --qos=6hours -c 10 --mem=64G kmer_genomescope.SRA.sh Trichomycterus_rosablanca
sbatch --qos=6hours -c 10 --mem=64G kmer_genomescope.SRA.sh Clarias_gariepinus
sbatch --qos=6hours -c 10 --mem=64G kmer_genomescope.SRA.sh Plotosus_lineatus
sbatch --qos=6hours -c 10 --mem=64G kmer_genomescope.SRA.sh Silurus_meridionalis
sbatch --qos=6hours -c 10 --mem=64G kmer_genomescope.SRA.sh Synodontis_membranacea
sbatch --qos=6hours -c 10 --mem=64G kmer_genomescope.SRA.sh Cranoglanis_bouderius
sbatch --qos=6hours -c 10 --mem=64G kmer_genomescope.SRA.sh Neoarius_graeffei
sbatch --qos=6hours -c 10 --mem=64G kmer_genomescope.SRA.sh Ameiurus_melas


for species_ID in `cat species_SRA.txt` ; do wc -l $species_ID.GenomeScope_rslt/summary.txt ; done

#Extract results from genomescopre

cat list_sp.txt species_SRA.txt > all_species_list.txt

for species_ID in `cat all_species_list.txt` ; do

	homo_min=`grep "Homozygous" $species_ID.GenomeScope_rslt/summary.txt | sed 's/  */	/g' | cut -f3 | sed 's/%//g'`
	homo_max=`grep "Homozygous" $species_ID.GenomeScope_rslt/summary.txt | sed 's/  */	/g' | cut -f4 | sed 's/%//g'`
	hetero_min=`grep "Heterozygous" $species_ID.GenomeScope_rslt/summary.txt | sed 's/  */	/g' | cut -f3 | sed 's/%//g'`
	hetero_max=`grep "Heterozygous" $species_ID.GenomeScope_rslt/summary.txt | sed 's/  */	/g' | cut -f4 | sed 's/%//g'`
	genome_size_min=`grep "Genome Haploid Length" $species_ID.GenomeScope_rslt/summary.txt | sed 's/  */	/g' | cut -f4 | sed 's/,//g'`
	genome_size_max=`grep "Genome Haploid Length" $species_ID.GenomeScope_rslt/summary.txt | sed 's/  */	/g' | cut -f6 | sed 's/,//g'`
	model_fit_min=`grep "Model Fit" $species_ID.GenomeScope_rslt/summary.txt | sed 's/  */	/g' | cut -f3 |  sed 's/%//g'`
	model_fit_max=`grep "Model Fit" $species_ID.GenomeScope_rslt/summary.txt | sed 's/  */	/g' | cut -f4 |  sed 's/%//g'`

	echo "$species_ID,$homo_min,$homo_max,$hetero_min,$hetero_max,$genome_size_min,$genome_size_max,$model_fit_min,$model_fit_max" 
done > GenomeScope_summary.csv




