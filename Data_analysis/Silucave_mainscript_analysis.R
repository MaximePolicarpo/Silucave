##### Libraries  ---------------------------------

rm(list=ls())

set.seed(2712)
options(scipen=0, digits=7)
options(bitmapType = "cairo")

library("caper")
library("patchwork")
library("ape")
library(dplyr)
library(ggplot2)
library(tidyverse)
library("lattice")
library(reshape2)
library(adephylo)
library(phylobase)
library(data.table)
library(phytools)
library("corrplot")
library("geiger")
library("ggpubr")
library(RColorBrewer)
library(gt)
library(scales)
library(phangorn)
library(castor)
library(VennDiagram)
library(venn)
library(ggtree)
library(ggtreeExtra)
library("viridis")
library(ggnewscale)
library("topGO")
library(GOstats)
library(rstatix)

split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}





##### Data load - Species trees   ---------------------------------

species_tree <- 
  read.tree("AMAS_concatenated_alignment_BUSCO.fa.timetree.nwk")

species_tree <- drop.tip(species_tree, "Channallabes_apus")

species_tree <- makeNodeLabel(species_tree, method="number", prefix="Node")
species_list <- species_tree$tip.label

species_df <- as.data.frame(species_list)
colnames(species_df) <- "species"

#Add habitat to species df

cave_species <- 
  c("Prietella_phreatophila",
    "CHM6",
    "CSV83",
    "CUL4",
    "CUL9",
    "Trichomycterus_rosablanca"
  )

species_df <- 
  species_df %>%
  mutate(Habitat = if_else(
    species %in% cave_species,
    "Cave", 
    "Surface"
  ))


species_tree_len <- species_tree
species_tree_len$edge.length <- species_tree_len$edge.length * 1000
##### Import OGG - gene names   ---------------------------------

OGG_names_df <- 
  read.table("OGG_Names_Product.tsv",
             header=FALSE,
             sep="\t")
colnames(OGG_names_df) <- c("OGG", "gene_ID", "gene_name", "evalue_name")


##### Mitochondrial phylogeny   ---------------------------------

cox1_nt_tree <- read.tree("cox1.nuc.aln.trimal.treefile")
cox1_nt_tree_rooted <- midpoint_root(cox1_nt_tree)

cytb_nt_tree <- read.tree("cytb.nuc.aln.trimal.treefile")
cytb_nt_tree_rooted <- midpoint_root(cytb_nt_tree)

#### BUSCO analysis  ---------------------------------

BUSCO_df <- 
  read.table("BUSCO_miniprot.results.csv",
             sep="\t",
             header=FALSE)

colnames(BUSCO_df) <- c("species", "B_Complete", "B_Complete_single", "B_Complete_duplicated",
                        "B_Fragmented", "B_Missing", "B_Pseudogene")


BUSCO_df <- left_join(BUSCO_df, species_df)
BUSCO_df <- BUSCO_df %>% mutate(B_Complete_single_ps = B_Complete_single + B_Pseudogene)


BUSCO_df_long <- 
  as.data.frame(BUSCO_df %>%
                  dplyr::select(species, B_Complete_single, B_Complete_duplicated, 
                                B_Fragmented, B_Missing, B_Pseudogene) %>%
                  pivot_longer(!species, names_to = "category", values_to = "count"))



BUSCO_df_long$category <-
  factor(BUSCO_df_long$category,
         levels=rev(c("B_Complete_single", "B_Complete_duplicated", 
                      "B_Fragmented", "B_Pseudogene","B_Missing")))


BUSCO_df_long %>%
  ggplot(., aes(x=species, y=count, fill=category)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  scale_fill_manual(values =
                      c("B_Complete_single"="#56B4E9",
                        "B_Complete_duplicated"="#009E73", 
                        "B_Fragmented"="#D55E00",
                        "B_Missing"="gray",
                        "B_Pseudogene" = "black")) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


##### Comparison iqtree vs astralIII   ---------------------------------

astral_phylo <- read.tree("AMAS_concatenated_alignment_BUSCO.fa.treefile.rooted")
astral_phylo <- drop.tip(astral_phylo, "Channallabes_apus")

iqtree_phylo <- read.tree("IQTREE_topology/AMAS_concatenated_alignment_BUSCO.fa.treefile.rooted")
iqtree_phylo <- drop.tip(iqtree_phylo, "Channallabes_apus")


astral_iqtree_df <- as.data.frame(iqtree_phylo$tip.label)
colnames(astral_iqtree_df) <- c("astral")
astral_iqtree_df <- astral_iqtree_df %>% mutate(iqtree = astral)


p1 <- ggtree(astral_phylo) 
p2 <- ggtree(iqtree_phylo)

d1 <- p1$data
d2 <- p2$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1

pp <- p1 + geom_tree(data=d2) +      
  ggnewscale::new_scale_fill() 

dd <- bind_rows(d1, d2) %>% 
  filter(!is.na(label))


pp + geom_line(aes(x, y, group=label), data=dd, color='grey') +
  geom_tiplab() +
  geom_tiplab(data = d2, hjust=1)



#### GenomeScope results   ---------------------------------

genomescope_df <- 
  read.table("GenomeScope_summary.csv",
             sep=",",
             header=FALSE)

colnames(genomescope_df) <- 
  c("species", "homo_min", "homo_max", "hetero_min", "hetero_max", "genome_size_min",
    "genome_size_max", "model_fit_min", "model_fit_max")

genomescope_df <- left_join(genomescope_df, species_df, by="species")

genomescope_df$genome_size_max <- as.numeric(genomescope_df$genome_size_max)
genomescope_df$hetero_max <- as.numeric(genomescope_df$hetero_max)
genomescope_df$hetero_min <- as.numeric(genomescope_df$hetero_min)


caper_genomeP <- 
  comparative.data(phy = species_tree, 
                   data = genomescope_df,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)


Heterozyogisity_vs_Habitat <-
  pgls(hetero_max ~ Habitat, 
       data = caper_genomeP, 
       lambda = 1)
summary(Heterozyogisity_vs_Habitat)

GenomeSize_vs_Habitat <-
  pgls(genome_size_max ~ Habitat, 
       data = caper_genomeP, 
       lambda = "ML")
summary(GenomeSize_vs_Habitat)


genomescope_df %>%
  group_by(Habitat) %>%
  summarise(mean_H = mean(hetero_max))


#Compare GenomeScope genome sizes and assembly sizes

assembly_sizes_df <- 
  read.table("assembly_sizes.csv", sep=",", header=FALSE)
colnames(assembly_sizes_df) <- c("species", "assembly_size")

scope_vs_assem <- 
  left_join(genomescope_df, assembly_sizes_df, by="species") %>%
  filter(! is.na(assembly_size))

scope_vs_assem_cor <- cor.test(scope_vs_assem$genome_size_max, scope_vs_assem$assembly_size,
                               method="pearson")


scope_vs_assem %>%
  ggplot(., aes(x=genome_size_max, y=assembly_size, color=Habitat)) +
  geom_point(size=2) +
  geom_abline(intercept = 0, slope=1, color="gray", linetype="dashed") +
  theme_classic() +
  scale_color_manual(values = c("Cave" = "#CC79A7", "Surface" = "#56B4E9")) + 
  geom_smooth(method = "lm", color = "black") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Genome size - GenomeScope2") +
  ylab("Genome assembly size")

#Hetero vs genome size
                                             
cor.test(genomescope_df$hetero_max, genomescope_df$genome_size_max, method="pearson")
genomescope_df_noout <- genomescope_df %>% filter(species != "Neoarius_graeffei")
cor.test(genomescope_df_noout$hetero_max, genomescope_df_noout$genome_size_max, method="pearson")

genomescope_df %>%
  ggplot(., aes(y=hetero_max, x=genome_size_max, color=Habitat)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = c("Cave" = "#CC79A7", "Surface" = "#56B4E9")) + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  ylab("Heterozygosity") +
  xlab("Genome size (bp)")
                                             
#Compare GenomeScope genome sizes and animal genome size database


genomedb_sizes_df <- 
  read.table("genomesize_animaldb.tsv", sep="\t", header=TRUE)

genomedb_sizes_df <- 
  genomedb_sizes_df %>%
  rowwise() %>%
  mutate(Mb = pg * 0.978 * 10^9)


scope_vs_animaldb <- 
  left_join(genomescope_df, genomedb_sizes_df, by="species") %>%
  filter(! is.na(Mb))

scope_vs_animaldb %>% filter(species_db == species)
scope_vs_assem_cor <- cor.test(scope_vs_animaldb$genome_size_max, scope_vs_animaldb$Mb,
                               method="pearson")


scope_vs_animaldb %>%
  ggplot(., aes(x=genome_size_max, y=Mb, color=Habitat)) +
  geom_point(size=2) +
  geom_abline(intercept = 0, slope=1, color="gray", linetype="dashed") +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = c("Cave" = "#CC79A7", "Surface" = "#56B4E9")) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Genome size - GenomeScope2") +
  ylab("Genome size - Animal GSD")


#Add species with no reads based on their genome assembly size

Ancistrus_triradiatus_GS <-
  as.data.frame(cbind("Ancistrus_triradiatus", NA, NA, NA, NA, 993800000, 993800000, NA, NA, "Surface"))

colnames(Ancistrus_triradiatus_GS) <- colnames(genomescope_df)

Corydoras_maculifer_GS <- 
  as.data.frame(cbind("Corydoras_maculifer", NA, NA, NA, NA, 636400000, 636400000, NA, NA, "Surface"))
colnames(Corydoras_maculifer_GS) <- colnames(genomescope_df)

Pangasius_djambal_GS <- 
  as.data.frame(cbind("Pangasius_djambal", NA, NA, NA, NA, 839800000, 839800000, NA, NA, "Surface"))
colnames(Pangasius_djambal_GS) <- colnames(genomescope_df)



genomescope_df <-
  rbind(genomescope_df, Ancistrus_triradiatus_GS, Corydoras_maculifer_GS, Pangasius_djambal_GS)

genomescope_df$genome_size_max <- as.numeric(genomescope_df$genome_size_max)
genomescope_df$hetero_max <- as.numeric(genomescope_df$hetero_max)

#Draw a tree with the results

p <- 
  ggtree(species_tree, size=2) %<+% 
  ace_values_labels_all +
  aes(color=as.numeric(Cave)) +
  scale_color_viridis() + 
  ggnewscale::new_scale_colour() + 
  geom_tiplab(size = 3, offset = 0.01) +
  #theme(legend.position="none") +
  new_scale_fill() +
  new_scale_color()



p + 
  geom_facet(
    panel="Genome_Size",
    data= genomescope_df,
    geom = geom_point,
    mapping = aes(x=1, color=hetero_max, size=genome_size_max),
    width = 0.5
  ) +
  scale_color_distiller(palette = "YlOrBr", direction = 1) +
  theme(legend.position = "none") +
  xlim_tree(0.5)
  

genomescope_df %>% dplyr::select(species, genome_size_max) %>% arrange(genome_size_max)

#Plot the heterozygosity vs genome size

genomescope_df %>%
  ggplot(., aes(x=genome_size_max, y=hetero_max, color=Habitat)) +
  geom_point(size = 3) + 
  scale_color_manual(values = c("Surface" = "black", "Cave" = "#FFC107")) +
  xlab("Genome size") +
  ylab("Heterozygosity") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")


#### Load Vision genes tables  ---------------------------------

#Load the sequences table
vision_sequences_df <- 
  read.table("Table_vision_seqs.tsv",
             sep="\t",
             header=FALSE)
colnames(vision_sequences_df) <- 
  c("gene_clade","species","gene_name","genomic_position","CDS_length","gene_type","note",
    "exon_Count","sequence")
vision_sequences_df <- vision_sequences_df %>% filter(species != "") %>% filter(gene_type != "Not found")

#Load the LoF tables, with the read count frequencies
stop_df <- read.table("Stop_table.tsv", sep="\t",header=FALSE)
colnames(stop_df) <- c("gene_name", "species", "lof_type", "aa_pos")
stop_df <- stop_df %>% mutate(cds_pos = aa_pos*3) %>% mutate(nb_nuc = "") %>% dplyr::select(c("gene_name", "species", "lof_type", "cds_pos", "nb_nuc", "aa_pos"))

frameshift_df <- read.table("Frameshift_table.tsv", sep="\t",header=FALSE)
colnames(frameshift_df) <- c("gene_name", "species", "lof_type", "cds_pos", "nb_nuc")

other_lof_df <- read.table("Other_LoF_table.tsv", sep="\t",header=FALSE)
colnames(other_lof_df) <- c("gene_name", "species", "lof_type", "cds_pos")
other_lof_df <- other_lof_df %>% mutate(nb_nuc = "") 

read_counts_df <- read.table("Read_counts_table.tsv", sep="\t",header=FALSE)
colnames(read_counts_df) <- c("gene_name", "species", "lof_type", "aa_pos", "position", "support_reads", "non_support_reads")

#Merge tables
stop_df <- left_join(stop_df, read_counts_df, by=c("gene_name", "species", "lof_type", "aa_pos"))
stop_df <- stop_df %>% dplyr::select(-aa_pos)
colnames(read_counts_df) <- c("gene_name", "species", "lof_type", "cds_pos", "position", "support_reads", "non_support_reads")
frameshift_df <- left_join(frameshift_df, read_counts_df, by=c("gene_name", "species", "lof_type", "cds_pos"))
other_lof_df <- left_join(other_lof_df, read_counts_df, by=c("gene_name", "species", "lof_type", "cds_pos"))

vision_lof_df <- rbind(stop_df, frameshift_df, other_lof_df)
vision_lof_df <- 
  vision_lof_df %>% 
  mutate(lof_name_temp = paste(gene_name,lof_type, sep="_")) %>%
  mutate(lof_name = paste(lof_name_temp,cds_pos, sep="_")) %>%
  dplyr::select(-lof_name_temp)

#For Ancistrus_triradiatus, no reads are available. All LoF will be considered as homozygous. 

vision_lof_df <- 
  vision_lof_df %>% 
  mutate(support_reads = ifelse(is.na(support_reads), 10, support_reads)) %>%
  mutate(non_support_reads = ifelse(is.na(non_support_reads), 0, non_support_reads))
  

#Verify that the LoF table and Sequence table are similar in term of pseudogene nb


count_seq_table <- 
  as.data.frame(
  vision_sequences_df %>%
    filter(gene_type == "Pseudogene") %>%
    group_by(species) %>%
    summarise(countP1 = n())
)


count_lof_table <-
  as.data.frame(
    vision_lof_df %>%
    dplyr::select(species, gene_name) %>%
    distinct() %>%
    group_by(species) %>%
    summarise(countP2 = n())
)


left_join(count_seq_table, count_lof_table, by="species") %>% filter(countP2 != countP1)



#Verify that the LoF table and Read count table are similar in term of pseudogene nb


count_lof_table_second <-
  as.data.frame(
    read_counts_df %>%
      dplyr::select(species, gene_name) %>%
      distinct() %>%
      group_by(species) %>%
      summarise(countP3 = n())
  )


left_join(count_lof_table_second, count_lof_table, by="species") %>% filter(countP2 != countP3)


#Compute the genotype for each LoF mutation using a binomial law. 

vision_lof_df <- 
  vision_lof_df %>% 
  rowwise() %>%
  mutate(total_reads = non_support_reads + support_reads) %>%
  mutate(non_support_reads_prop = non_support_reads/total_reads) %>%
  mutate(support_reads_prop = support_reads/total_reads) %>% 
  mutate(minimum_read_number = min(support_reads, non_support_reads)) %>%
  mutate(binom_proba = pbinom(minimum_read_number, total_reads , p = 0.5)) %>%
  mutate(Genotype = case_when(
    binom_proba >= 0.05 & binom_proba <= 0.95 & non_support_reads != 0 & support_reads != 0 ~ "Heterozygous",
    binom_proba >= 0.05 & binom_proba <= 0.95 & non_support_reads == 0  & support_reads != 0 ~ "LoF_Homozygous",
    binom_proba >= 0.05 & binom_proba <= 0.95 & non_support_reads != 0  & support_reads == 0 ~ "no_LoF_Homozygous",
    binom_proba < 0.05 & (support_reads > non_support_reads) ~ "LoF_Homozygous",
    binom_proba > 0.95 & (support_reads > non_support_reads) ~ "LoF_Homozygous",
    binom_proba < 0.05 & (support_reads < non_support_reads) ~ "no_LoF_Homozygous",
    binom_proba > 0.95 & (support_reads < non_support_reads) ~ "no_LoF_Homozygous",
  )) %>%
  mutate(LoF_genotype = case_when(
    Genotype == "LoF_Homozygous" ~ 2,
    Genotype == "Heterozygous"  ~ 1,
    Genotype == "no_LoF_Homozygous" ~ 0
  )) %>%
  mutate(non_LoF_genotype = case_when(
    Genotype == "LoF_Homozygous" ~ 0,
    Genotype == "Heterozygous"  ~ 1,
    Genotype == "no_LoF_Homozygous" ~ 2
  ))  %>% 
  mutate(total_allele_nb = LoF_genotype + non_LoF_genotype) 

vision_lof_df <- as.data.frame(vision_lof_df)

vision_lof_df %>% filter(Genotype == "no_LoF_Homozygous")


#Count the number of LoF mutations per species

vision_lof_df %>%
  group_by(species, Genotype) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

vision_lof_df %>%
  pull(Genotype) %>% unique()

vision_lof_df_observed <- 
  vision_lof_df %>%
  filter(Genotype != "no_LoF_Homozygous")



vision_lof_df_observed_caves <- 
  vision_lof_df_observed %>%
  filter(species %in% c("Prietella_phreatophila",
                        "Trichomycterus_rosablanca",
                        "CUL9",
                        "CSV83",
                        "CUL4", 
                        "CHM6"))

vision_lof_df_observed_caves$orthogroup <- vision_lof_df_observed_caves$gene_name
vision_lof_df_observed_caves$individual <- vision_lof_df_observed_caves$species
vision_lof_df_observed_caves <- vision_lof_df_observed_caves %>% mutate(prot_pos = round(cds_pos/3))

vision_lof_df_observed_caves_export <- 
  vision_lof_df_observed_caves %>%
  dplyr::select(orthogroup, gene_name, species, individual, lof_type, prot_pos, nb_nuc, Genotype)



#Merge the LoF table with the Sequence table. Consider a gene as pseudogene
#if there is at-least one LoF present at heterozygous/homozygous state.
#Consider a pseudogene as heterozygous if there is no homozygous LoF 
#and homozygous if there is at-least one homozygous mutation

vision_sequences_df_final <- as.data.frame(NULL)
for (curr_line in 1:nrow(vision_sequences_df)){
  curr_species <- vision_sequences_df[curr_line,]$species
  curr_gene_name <- vision_sequences_df[curr_line,]$gene_name
  curr_gene_type <- vision_sequences_df[curr_line,]$gene_type
  
  nb_homozygous_lof <- 
    nrow(vision_lof_df %>% 
    filter(species == curr_species) %>% 
    filter(gene_name == curr_gene_name) %>%
    filter(Genotype == "LoF_Homozygous"))
  
  nb_heterozygous_lof <- 
    nrow(vision_lof_df %>% 
    filter(species == curr_species) %>% 
    filter(gene_name == curr_gene_name) %>%
    filter(Genotype == "Heterozygous"))
  
  
  if(nb_homozygous_lof > 0){
    curr_line_mut <- vision_sequences_df[curr_line,] %>% mutate(corrected_gene_type = "Pseudogene_Homozygous") 
  } 
  
  if(nb_homozygous_lof == 0 & nb_heterozygous_lof > 0){
    curr_line_mut <- vision_sequences_df[curr_line,] %>% mutate(corrected_gene_type = "Pseudogene_Heterozygous") 
  }   
  
  if(nb_homozygous_lof == 0 & nb_heterozygous_lof == 0 & curr_gene_type == "Pseudogene"){
    curr_line_mut <- vision_sequences_df[curr_line,] %>% mutate(corrected_gene_type = "Complete") 
  }   
  
  if(nb_homozygous_lof == 0 & nb_heterozygous_lof == 0 & curr_gene_type == "Incomplete"){
    curr_line_mut <- vision_sequences_df[curr_line,] %>% mutate(corrected_gene_type = "Incomplete") 
  }  
  
  if(nb_homozygous_lof == 0 & nb_heterozygous_lof == 0 & curr_gene_type == "Complete"){
    curr_line_mut <- vision_sequences_df[curr_line,] %>% mutate(corrected_gene_type = "Complete") 
  }  
  
  vision_sequences_df_final <- rbind(vision_sequences_df_final, curr_line_mut)
  
}


#Check which gene_type changed

nrow(vision_sequences_df_final)
nrow(vision_sequences_df)


vision_sequences_df_final %>%
  filter(corrected_gene_type != gene_type) %>%
  dplyr::select(species, gene_name, gene_type, corrected_gene_type) 


vision_sequences_df_final_export <- 
  vision_sequences_df_final %>%
  filter(species %in% c("Prietella_phreatophila",
                        "Trichomycterus_rosablanca",
                        "CUL9",
                        "CSV83",
                        "CUL4", 
                        "CHM6")) %>%
  mutate(individual = species) %>%
  dplyr::select(gene_clade, species, individual, gene_name, genomic_position, corrected_gene_type)


#### Vision genes analysis  ---------------------------------

#Summary of gene types 
summary_vision_df <- 
  as.data.frame(
    vision_sequences_df_final %>%
  group_by(species, corrected_gene_type) %>%
  summarise(count = n())
  )

summary_vision_df_wide <- 
  as.data.frame(
    summary_vision_df %>%
      pivot_wider(names_from = corrected_gene_type, values_from = count, values_fill = 0)
  ) %>%
  mutate(Total = Complete + Incomplete + Pseudogene_Homozygous +  Pseudogene_Heterozygous ) %>%
  mutate(Pseudogene = Pseudogene_Homozygous +  Pseudogene_Heterozygous ) %>%
  mutate(prop_pseudo = Pseudogene / Total)
  
summary_vision_df_wide <- left_join(summary_vision_df_wide, species_df, by="species")

#Compute the total number of LoF mutations

nb_lof_df <- 
  as.data.frame(vision_lof_df %>%
  filter(Genotype %in% c("LoF_Homozygous", "Heterozygous")) %>%
  group_by(species) %>%
  summarise(lof_nb = n()))

#Compute the total number of vision nucleotides

nb_nuc_df <- 
  as.data.frame(vision_sequences_df_final %>%
  group_by(species) %>%
  summarise(nt_nb = sum(CDS_length)))

#Combine tables

summary_vision_df_wide <- left_join(summary_vision_df_wide, nb_lof_df, by="species")
summary_vision_df_wide <- left_join(summary_vision_df_wide, nb_nuc_df, by="species")
summary_vision_df_wide <- summary_vision_df_wide %>% mutate(lof_nb = ifelse(is.na(lof_nb), 0, lof_nb))
summary_vision_df_wide <- summary_vision_df_wide %>% mutate(lof_per_nt = lof_nb/nt_nb)


#Print the table

gt(summary_vision_df_wide %>%
     dplyr::select(species, Complete, Incomplete, Pseudogene, Total, prop_pseudo, nt_nb, lof_nb, lof_per_nt) %>%
     arrange(prop_pseudo))


#### Compute stats on observed LoF mutations ---------------------------------


### What are the number and type of observed LoF ? 


vision_lof_df %>%
  filter(Genotype != "no_LoF_Homozygous") %>%
  dplyr::select(gene_name, lof_type, cds_pos) %>%
  distinct() %>%
  group_by(lof_type) %>%
  summarise(count = n()) %>%
  ggplot(., aes(x="", y=count, fill=lof_type)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  scale_fill_manual(values = c(
    "Ins" = "#000000",
    "Del" = "gray",
    "STOP" = "#D55E00",
    "Splice_site" = "#CC79A7",
    "Start_loss" = "#009E73",
    "Stop_loss" = "#0072B2"
  )) + 
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x=element_blank(),
    plot.title=element_text(size=14, face="bold")
  ) +
  theme(legend.position="none") 



### What is the distribution of indel sizes ? 


vision_lof_df %>%
  filter(Genotype != "no_LoF_Homozygous") %>%
  filter(lof_type %in% c("Del", "Ins")) %>%
  dplyr::select(gene_name, lof_type, cds_pos, nb_nuc) %>%
  distinct() %>%
  group_by(lof_type, nb_nuc) %>%
  summarise(count = n()) %>%
  ggplot(., aes(x=as.numeric(nb_nuc), y=count, fill=lof_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c(
    "Ins" = "#000000",
    "Del" = "gray"
  )) + 
  theme_classic() +
  theme(legend.position="none") +
  xlab("Number of nucleotides") +
  ylab("count") +
  scale_x_continuous(breaks=c(0,1,2,4,5,8,10, 13,14, 16, 20,43))


#### pGLS pseudogene vs Habitat -- BUSCO + Vision ---------------------------------


BUSCO_df_long <- 
  as.data.frame(BUSCO_df %>%
                  dplyr::select(species, B_Complete, B_Fragmented, B_Missing, B_Pseudogene) %>%
                  pivot_longer(!species, names_to = "category", values_to = "count")
                )


summary_vision_df_wide_long <- 
  as.data.frame(summary_vision_df_wide %>%
                  dplyr::select(species, Complete, Pseudogene_Homozygous, Pseudogene_Heterozygous, Incomplete) %>%
                  pivot_longer(!species, names_to = "category", values_to = "count")
  )


BUSCO_df_long$category <-
  factor(BUSCO_df_long$category,
         levels=(c("B_Complete", "B_Fragmented", "B_Pseudogene","B_Missing")))
summary_vision_df_wide_long$category <-
  factor(summary_vision_df_wide_long$category,
         levels=(c("Complete", "Incomplete", "Pseudogene_Homozygous","Pseudogene_Heterozygous")))



species_tree_len_woDr <- drop.tip(species_tree_len, "Danio_rerio") 


ace_values_labels_all_g <- 
  ace_values_labels_all %>% mutate(Habitat = if_else(
    Cave > 0,
    "Cave",
    "Surface"
  ))


ggtree(species_tree_len_woDr, size=2) %<+% 
  ace_values_labels_all_g +
  aes(color=Habitat) +
  scale_color_manual(values = c("Cave" = "#CC79A7", "Surface" = "#0072B2")) + 
  #scale_color_viridis() + 
  ggnewscale::new_scale_colour() + 
  geom_tiplab(size = 3, offset = 0.01) + 
  theme(legend.position='none')  +
  geom_facet(panel = 'BUSCO', 
             data = BUSCO_df_long, geom = geom_bar, 
             mapping = aes(x = as.numeric(count), fill = as.factor(category)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values =
                      c("B_Complete"="#009E73",
                        "B_Fragmented"="#56B4E9",
                        "B_Missing"="gray",
                        "B_Pseudogene" = "#D55E00"))  +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = 'Vision', 
             data = summary_vision_df_wide_long, geom = geom_bar, 
             mapping = aes(x = as.numeric(count), fill = as.factor(category)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values =
                      c("Complete"="#009E73",
                        "Incomplete"="#56B4E9",
                        "Pseudogene_Homozygous"="#D55E00",
                        "Pseudogene_Heterozygous" = "#E69F00"))  +
  theme_tree2(legend.position = 'none') +
  xlim_tree(500)



summary_vision_df_wide <- left_join(summary_vision_df_wide, BUSCO_df, by=c("species", "Habitat"))
summary_vision_df_wide <- summary_vision_df_wide %>% mutate(B_prop_pseudo = B_Pseudogene / (B_Complete + B_Fragmented + B_Pseudogene))




caper_busco_vision <- 
  comparative.data(phy = species_tree, 
                   data = summary_vision_df_wide,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)

Vision_propPseudo_vs_Habitat <-
  pgls(prop_pseudo ~ Habitat, 
       data = caper_busco_vision, 
       lambda = "ML")
summary(Vision_propPseudo_vs_Habitat)


BUSCO_propPseudo_vs_Habitat <-
  pgls(B_prop_pseudo ~ Habitat, 
       data = caper_busco_vision, 
       lambda = "ML")
summary(BUSCO_propPseudo_vs_Habitat)



summary_vision_df_wide %>% 
  filter(Habitat == "Surface") %>% 
  filter(species != "Danio_rerio") %>%
  pull(B_prop_pseudo) %>% max()


summary_vision_df_wide %>% 
  filter(Habitat == "Cave") %>% 
  filter(species != "Danio_rerio") %>%
  pull(B_prop_pseudo) %>% min()


#### Tile plot of Vision pseudogenes ---------------------------------



uniq_genes <- 
  vision_sequences_df_final %>%
  filter(species != "Danio_rerio") %>%
  group_by(gene_name) %>%
  summarise(count = n()) %>%
  filter(count > 0) %>%
  pull(gene_name) %>%
  unique()

uniq_ind <- 
  vision_sequences_df_final %>%
  pull(species) %>%
  unique()

all_combinations <- expand.grid(species = uniq_ind, gene_name = uniq_genes)


vision_sequences_df_final_expanded <- 
  all_combinations %>%
  left_join(vision_sequences_df_final, by = c("species", "gene_name")) %>%
  mutate(corrected_gene_type = ifelse(is.na(corrected_gene_type), "Not_found", corrected_gene_type))


order_genes_plot <- 
  c("rh1.1","rh1.2","rh2.1","lws1","exorh","opn3","opn4m2","opn4m3",
    "opn4x1","opn5","opn6a","opn6b","opn7a","opn7b","opn7c","opn8a",
    "parapinopsin-1","parietopsin","rgr1","rgr2","rrh","tmt1a","tmt3a","tmt3b",
    "va1","va2","cryaa","cryba1b","cryba1l1","cryba2a","cryba4","crybb1","crybb1l1",
    "crybb1l2","crybgx","crygm5","crygn2","rpe65a","arr3a","arr3b","saga","sagb","gcap1",
    "gcap2","gcap3","gcap4","gcap7.1","gcap7.2","grk1a","grk1b","grk7a","rcv1a","rcv1b",
    "rcv2b","gc2.1","gc2.2","gucy2f","gc3","pde6a","pde6b","pde6c","pde6ga","pde6gb",
    "pde6ha","pde6hb","gnat1","gnat2","gnb3a","gnb3b","gngt1","gngt2","gja8b")


vision_sequences_df_final_expanded$gene_name <- 
  factor(vision_sequences_df_final_expanded$gene_name,
         levels = order_genes_plot)

vision_sequences_df_final_expanded$corrected_gene_type <- 
  factor(vision_sequences_df_final_expanded$corrected_gene_type,
         levels = c("Complete", "Incomplete", "Pseudogene_Homozygous", "Pseudogene_Heterozygous", "Not_found")) 


species_tree_len_woDr <- drop.tip(species_tree_len, "Danio_rerio") 
species_tree_len_woDr_notir <- species_tree_len_woDr
species_tree_len_woDr_notir$tip.label <- gsub("_", " ", species_tree_len_woDr_notir$tip.label)
ace_values_labels_all_notir <- ace_values_labels_all
ace_values_labels_all_notir$label <- gsub("_", " ", ace_values_labels_all_notir$label)
ace_values_labels_all_notir <- ace_values_labels_all_notir %>% dplyr::select(label, Cave)
genomescope_df_notir <- genomescope_df
genomescope_df_notir$species <- gsub("_", " ", genomescope_df_notir$species)
vision_sequences_df_final_expanded_notir <- vision_sequences_df_final_expanded
vision_sequences_df_final_expanded_notir$species <- gsub("_", " ", vision_sequences_df_final_expanded_notir$species)

ace_values_labels_all_notir <- 
  ace_values_labels_all_notir %>% mutate(Habitat = if_else(
  Cave > 0,
  "Cave",
  "Surface"
))

p1 <- 
  ggtree(species_tree_len_woDr_notir, size=2) %<+% 
  ace_values_labels_all_notir +
  #aes(color=as.numeric(Cave)) +
  aes(color=Habitat) +
  scale_color_manual(values = c("Cave" = "#CC79A7", "Surface" = "#0072B2")) + 
  #scale_color_viridis() + 
  ggnewscale::new_scale_colour() + 
  geom_tiplab(size = 3, offset = 0.01, fontface = 3) + 
  theme(legend.position='none') +
  geom_facet(
    panel="Genome_Size",
    data= genomescope_df_notir %>% dplyr::select(-Habitat),
    geom = geom_point,
    mapping = aes(x=1, color=hetero_max, size=genome_size_max),
    width = 0.5
  ) +
  scale_color_distiller(palette = "YlOrBr", direction = 1) +
  ggnewscale::new_scale_colour() + 
  geom_facet(panel = "Gene", 
             data = vision_sequences_df_final_expanded_notir, 
             geom = geom_tile, 
             aes(x = as.integer(gene_name), 
                 fill = corrected_gene_type), 
             width = 1, 
             color="black") +
  scale_fill_manual(values =
                      c("Complete"="#009E73",
                        "Incomplete"="#56B4E9",
                        "Pseudogene_Homozygous"="#D55E00",
                        "Pseudogene_Heterozygous" = "#E69F00",
                        "Not_found" = "gray"))  +
  theme(legend.position = "none",
        strip.text = element_blank()) +
  xlim_tree(300) 



p2 <- 
  ggtree(species_tree_len_woDr_notir, size=2) %<+% 
  ace_values_labels_all_notir +
  #aes(color=as.numeric(Cave)) +
  #scale_color_viridis() + 
  aes(color=Habitat) +
  scale_color_manual(values = c("Cave" = "#CC79A7", "Surface" = "#0072B2")) + 
  ggnewscale::new_scale_colour() + 
  geom_tiplab(size = 3, offset = 0.01, fontface = 3) + 
  geom_facet(
    panel="Genome_Size",
    data= genomescope_df_notir %>% dplyr::select(-Habitat),
    geom = geom_point,
    mapping = aes(x=1, color=hetero_max, size=genome_size_max),
    width = 0.5
  ) +
  scale_color_distiller(palette = "YlOrBr", direction = 1) +
  ggnewscale::new_scale_colour() + 
  geom_facet(panel = "Gene", 
             data = vision_sequences_df_final_expanded_notir, 
             geom = geom_tile, 
             aes(x = as.integer(gene_name), 
                 fill = corrected_gene_type), 
             width = 1, 
             color="black") +
  scale_fill_manual(values =
                      c("Complete"="#009E73",
                        "Incomplete"="#56B4E9",
                        "Pseudogene_Homozygous"="#D55E00",
                        "Pseudogene_Heterozygous" = "#E69F00",
                        "Not_found" = "gray"))  +
  theme(
        strip.text = element_blank()) +
  xlim_tree(300) 


#Check which genes have been pseudogenized in internal branches of the tree (common LOF)


lof_mult_species_df <- 
  as.data.frame(
    vision_lof_df %>%
  filter(Genotype != "no_LoF_Homozygous") %>%
  group_by(gene_name, lof_type, cds_pos) %>%
  summarise(count = n()) %>%
  filter(count >= 2)
  )
 
lof_mult_species_df_filt <- as.data.frame(NULL)
for(curr_line in 1:nrow(lof_mult_species_df)){
  curr_gene <- lof_mult_species_df[curr_line,]$gene_name
  curr_type <- lof_mult_species_df[curr_line,]$lof_type
  curr_pos <- lof_mult_species_df[curr_line,]$cds_pos
  
  curr_species <- vision_lof_df %>% 
    filter(gene_name == curr_gene & lof_type == curr_type & cds_pos == curr_pos) %>% pull(species)

  curr_species_str <- paste(curr_species, collapse = ',')
  
  curr_df <- lof_mult_species_df[curr_line,] %>% mutate(species = curr_species_str)
  lof_mult_species_df_filt <- rbind(lof_mult_species_df_filt, curr_df)
  
}

lof_mult_species_df_filt %>%
  dplyr::select(gene_name, species) %>%
  distinct() %>% arrange(species)

#### Tile plot of Vision LoF ---------------------------------

vision_lof_df_sim <-
  vision_lof_df %>%
  filter(Genotype != "no_LoF_Homozygous") %>%
  dplyr::select(species, lof_name, Genotype)


vision_lof_df_sim %>% filter(Genotype != "LoF_Homozygous")

uniq_species <- species_tree$tip.label
uniq_LoF <- vision_lof_df_sim$lof_name %>% unique()

all_combinations <- expand.grid(species = uniq_species, lof_name = uniq_LoF)


vision_lof_df_sim_expanded <- all_combinations %>%
  left_join(vision_lof_df_sim, by = c("species", "lof_name")) %>%
  mutate(Genotype = ifelse(is.na(Genotype), "no_LoF_Homozygous", Genotype))

vision_lof_df_sim_expanded$gene_name <- gsub("_.*", "", vision_lof_df_sim_expanded$lof_name)


vision_lof_df_sim_expanded$lof_name <- factor(vision_lof_df_sim_expanded$lof_name,
                                        levels = vision_lof_df_sim_expanded$lof_name %>% unique()) 
vision_lof_df_sim_expanded$Genotype <- factor(vision_lof_df_sim_expanded$Genotype,
                                   levels = c("no_LoF_Homozygous", "LoF_Homozygous", "Heterozygous")) 



order_genes_plot <- 
  c("rh1.1","rh1.2","rh2.1","rh2.2","exorh","lws1","opn3","opn4m2","opn4m3",
    "opn4x1","opn5","opn6","opn6a","opn6b","opn7a","opn7b","opn7c","opn8a",
    "parapinopsin-1","parietopsin","rgr1","rgr2","rrh","tmt1a","tmt3a","tmt3b",
    "va1","va2","cryaa","cryba1b","cryba1l1","cryba2a","cryba4","crybb1","crybb1l1",
    "crybb1l2","crybgx","crygm5","crygn2","rpe65a","arr3a","arr3b","saga","sagb","gcap1",
    "gcap2","gcap3","gcap4","gcap7.1","gcap7.2","grk1a","grk1b","grk7a","rcv1a","rcv1b",
    "rcv2b","gc2.1","gc2.2","gucy2f","gc3","pde6a","pde6b","pde6c","pde6ga","pde6gb",
    "pde6ha","pde6hb","gnat1","gnat2","gnb3a","gnb3b","gngt1","gngt2","gja8b")

vision_lof_df_sim_expanded_ordered <- as.data.frame(NULL)
for(curr_gene in order_genes_plot){
  curr_df <- vision_lof_df_sim_expanded %>% filter(gene_name == curr_gene)
  vision_lof_df_sim_expanded_ordered <- rbind(vision_lof_df_sim_expanded_ordered, curr_df)
}



species_tree_len_woDr <- drop.tip(species_tree_len, "Danio_rerio") 


p2 <-
  ggtree(species_tree_len_woDr, size=2) %<+% 
  ace_values_labels_all_g +
  aes(color=Habitat) +
  scale_color_manual(values = c("Cave" = "#CC79A7", "Surface" = "#0072B2")) +
  #scale_color_viridis() + 
  ggnewscale::new_scale_colour() + 
  geom_tiplab(size = 3, offset = 0.01) + 
  theme(legend.position='none') +
  geom_facet(panel = "Gene", 
             data = vision_lof_df_sim_expanded_ordered, 
             geom = geom_tile, 
             aes(x = as.integer(lof_name), 
                 fill = Genotype), 
             width = 1, 
             color="black") +
  scale_fill_manual(values = c("no_LoF_Homozygous" = "#009E73", 
                               "LoF_Homozygous" = "#D55E00",
                               "Heterozygous" = "#E69F00")) +
  theme(legend.position = "none",
        strip.text = element_blank()) +
  xlim_tree(300)  


lof_mult_species_df_filt %>%
  group_by(species) %>%
  summarise(count = n())
  
lof_mult_species_df_filt %>%
  filter(species %in% c("Ancistrus_triradiatus,RHP1,CUL9"))


lof_mult_species_df_filt %>%
  filter(species %in% c("RHP1,CUL9", "CUL9,RHP1"))

#### Compute a LoF rate based on observed LoF mutations ---------------------------------

#Transition/transversion ratio computed with PAML on a concatenated alignment of vision genes

R_ratio = 2.11911 

#Stop codon probability

#extract sequences of I. punctatus (=reference)
concatenated_vision_seq <- 
  vision_sequences_df_final %>%
  filter(species == "Ictalurus_punctatus") %>%
  pull(sequence) %>%
  paste(., collapse = '') %>%
  as.character(.)


n=0
x=0

for(codon_stop in seq(3, nchar(concatenated_vision_seq), 3)){
  codon_start <- codon_stop-2
  codon_seq <- substr(concatenated_vision_seq, codon_start, codon_stop)
  
  if(codon_seq == "CAA" | codon_seq == "CGA" | codon_seq == "CAG"){
    x <- x + R_ratio/((3*R_ratio)+6)
    n <- n + 3
  } else if(codon_seq == "GAA" | codon_seq == "AAA" | codon_seq == "AGA" | codon_seq == "TGT" | codon_seq == "TGC" | codon_seq == "AAG" | codon_seq == "GAG" | codon_seq == "TTG" | codon_seq == "TCG" | codon_seq == "GGA"){
    x <- x + 1.0/((3*R_ratio)+6)
    n <- n + 3
    
  } else if(codon_seq == "TTA" | codon_seq == "TCA" | codon_seq == "TAC" | codon_seq == "TAT" | codon_seq == "TTA" | codon_seq == "TAC"){
    x <- x + 2.0/((3*R_ratio)+6)
    n <- n + 3
    
  } else if(codon_seq == "TGG"){
    x <= x + (2.0*R_ratio)/((3*R_ratio)+6)
    n <- n + 3
    
  } else {
    
    x <- x + 0.0
    n <- n + 3
    
  }
  
  
  
}

proba_stop_gain <- x/(n/3)



vision_lof_df %>%
  filter(Genotype != "no_LoF_Homozygous") 


#Probability of Frameshift
nb_stop <- 
  length(
    vision_lof_df %>%
      filter(Genotype != "no_LoF_Homozygous") %>%
      filter(lof_type == "STOP") %>%
      pull(lof_name) %>%
      unique()
  )


nb_fs <- 
  length(
    vision_lof_df %>%
      filter(Genotype != "no_LoF_Homozygous") %>%
      filter(lof_type %in% c("Ins", "Del")) %>%
      pull(lof_name) %>%
      unique()
  )




proba_frameshift <- proba_stop_gain * (nb_fs/nb_stop)


#Proba of splice site mutation (take the mean number of exons over species)

  
  
intron_number <- 
  vision_sequences_df_final %>%
  mutate(intron_count = as.numeric(exon_Count) - 1) %>%
  group_by(species) %>%
  summarise(total_intron = sum(intron_count)) %>%
  pull(total_intron) %>%
  sum()


base_number <- 
  vision_sequences_df_final %>%
  group_by(species) %>%
  summarise(sum_CDS_length = sum(CDS_length)) %>%
  pull(sum_CDS_length) %>%
  sum()

proba_splice <- 4 * (intron_number / base_number)

#Proba of Start and Stop loss

proba_start_loss <- 
  3 * (
    vision_sequences_df_final %>% group_by(species) %>% summarise(count = n()) %>% pull(count) %>% sum()
    / base_number
  )
proba_stop_loss <- 0.852 * proba_start_loss


# Sum of LoF probas


LoF_proba <- proba_stop_loss + proba_start_loss + proba_splice + proba_frameshift + proba_stop_gain


#Make a nice table

LoF_rates_df <- 
  as.data.frame(
    cbind(proba_stop_gain, proba_frameshift, 
          proba_start_loss, proba_stop_loss, proba_splice, LoF_proba))

colnames(LoF_rates_df) <- 
  c("Stop_gain", "Frameshift", "Start_loss", "Stop_loss", "Splice", "Total")


#### Get the observed distribution of LoF per gene per species  ---------------------------------

lof_obs_per_sp_df <- 
  as.data.frame(
    vision_lof_df %>%
      filter(Genotype != "no_LoF_Homozygous")  %>%
      group_by(species, gene_name) %>%
      summarise(lof_nb = n())
  )

list_species <- lof_obs_per_sp_df$species %>% unique()

noLoF_obs_df <- as.data.frame(NULL)
for(curr_species in list_species){
  
  pseudo_prez <- 
    lof_obs_per_sp_df %>%
    filter(species == curr_species) %>%
    pull(gene_name)
  
  curr_df <- 
    vision_sequences_df_final %>%
    filter(species == curr_species) %>%
    filter(! gene_name %in% pseudo_prez) %>%
    mutate(lof_nb = 0) %>%
    dplyr::select(species, gene_name,lof_nb)
  
  noLoF_obs_df <- rbind(noLoF_obs_df, curr_df)
}



LoF_distrib_obs_df <- rbind(lof_obs_per_sp_df, noLoF_obs_df)


LoF_distrib_obs_df_summary <-
  as.data.frame(LoF_distrib_obs_df %>%
                  group_by(species, lof_nb) %>%
                  summarise(count = n()))


#Add the proportion

LoF_distrib_obs_df_summary_total <- 
  as.data.frame(
    LoF_distrib_obs_df_summary %>%
      group_by(species) %>%
      mutate(total = sum(count))
  )

LoF_distrib_obs_df_summary_total <- 
  LoF_distrib_obs_df_summary_total %>%
  mutate(prop_count = count / total)

#Finally, extend for each category to be present in each species (from 0 to max lof_nb)

LoF_distrib_obs_df_summary_total %>% pull(lof_nb) %>% max()
  
list_sp <- LoF_distrib_obs_df_summary_total$species %>% unique
LoF_categories <- LoF_distrib_obs_df_summary_total$lof_nb %>% unique()

all_combinations <- expand.grid(species = list_sp, lof_nb = LoF_categories)


LoF_distrib_obs_df_summary_total_exp <- 
  all_combinations %>%
  left_join(LoF_distrib_obs_df_summary_total, by = c("species", "lof_nb")) %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)),
         across(where(is.character), ~ replace_na(., "0")))


#lets plot

list_sp <- LoF_distrib_obs_df_summary_total_exp$species %>% unique



for(curr_sp in list_sp){
  
  
  p1 <- 
    LoF_distrib_obs_df_summary_total_exp %>%
    filter(species == curr_sp) %>%
    ggplot(., aes(x=as.numeric(lof_nb), y=count)) +
    geom_bar(stat="identity", position = "dodge") +
    theme_classic() +
    scale_x_continuous("# of LoF", 
                       c(0,1, 2, 3, 4, 5, 6, 7, 8)) +
    #scale_x_continuous("# of LoF", 
    #                   x_axis_ticks) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          axis.title.x=element_blank()) +
    theme(legend.position="none")  +
    xlab("Number of LoF mutation") +
    ylab("Number of genes") 
    
  
  print(p1)
  
}



list_cave_sp <- c("Prietella_phreatophila", "CUL9","Trichomycterus_rosablanca", "CHM6",
                  "CSV83", "CUL4")
for(curr_sp in list_cave_sp){
  
  
  p1 <- 
    LoF_distrib_obs_df_summary_total_exp %>%
    filter(species == curr_sp) %>%
    ggplot(., aes(x=as.numeric(lof_nb), y=prop_count)) +
    geom_bar(stat="identity", position = "dodge") +
    theme_classic() +
    scale_x_continuous("# of LoF", 
                       c(0,1, 2, 3, 4, 5, 6, 7, 8)) +
    ylim(0, 1) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          axis.title.x=element_blank()) +
    theme(legend.position="none") +
    xlab("Number of LoF mutation") +
    ylab("Proportion of genes") 
  
  print(p1)
  
}



LoF_distrib_obs_df_summary_total_exp %>%
  filter(species == "Prietella_phreatophila") %>%
  ggplot(., aes(x=as.numeric(lof_nb), y=prop_count)) +
  geom_bar(stat="identity", position = "dodge", fill="gray", color="black") +
  theme_classic() +
  scale_x_continuous("# of LoF", 
                     c(0,1, 2, 3, 4, 5, 6, 7, 8)) +
  ylim(0, 1) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  xlab("Number of LoF mutation") +
  ylab("Proportion of genes") 


LoF_distrib_obs_df_summary_total_exp %>%
  filter(species == "CUL9") %>%
  ggplot(., aes(x=as.numeric(lof_nb), y=prop_count)) +
  geom_bar(stat="identity", position = "dodge", fill="gray", color="black") +
  theme_classic() +
  scale_x_continuous("# of LoF", 
                     c(0,1, 2, 3, 4, 5, 6, 7, 8)) +
  ylim(0, 1) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  xlab("Number of LoF mutation") +
  ylab("Proportion of genes") 


LoF_distrib_obs_df_summary_total_exp %>%
  filter(species %in% c("CUL4", "CSV83", "CHM6", "Trichomycterus_rosablanca")) %>%
  ggplot(., aes(x=as.numeric(lof_nb), y=prop_count, fill=species)) +
  geom_bar(stat="identity", position = "dodge", color="black") +
  scale_fill_manual(values = c("Trichomycterus_rosablanca" = "#CC79A7",
                               "CUL4" = "#0072B2",
                               "CHM6" = "#E69F00",
                               "CSV83" = "#009E73")) +
  theme_classic() +
  scale_x_continuous("# of LoF", 
                     c(0,1, 2, 3, 4, 5, 6, 7, 8)) +
  ylim(0, 1) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  xlab("Number of LoF mutation") +
  ylab("Proportion of genes") 


#### Compare LoF distribution with neutral gene proportion - P. p  ---------------------------------

LoF_distrib_obs_df_summary_total_exp$species %>% unique()
LoF_distrib_obs_df_summary_total_exp %>%
  filter(species == "Prietella_phreatophila") %>%
  filter(lof_nb == 1) %>% pull(count)
  
LoF_distrib_obs_df_summary_total_exp %>%
  filter(species == "Prietella_phreatophila") %>%
  filter(lof_nb > 1) %>% pull(count) %>% sum()


LoF_distrib_obs_df_summary_total_exp %>%
  filter(species == "Prietella_phreatophila") %>%
  arrange(lof_nb)


#Total number of genes in P. P
Total_Number_Genes <- 
  nrow(vision_sequences_df_final %>% 
         filter(species == "Prietella_phreatophila")) 

#Number of pseudogenes
Total_Number_Pseudo <- 
  nrow(vision_sequences_df_final %>% 
         filter(species == "Prietella_phreatophila") %>%
         filter(corrected_gene_type == "Pseudogene_Homozygous")) #no heterozygous pseudo 

#Number of Loss-of-function mutations
Number_Lof <- 
  LoF_distrib_obs_df_summary_total_exp %>%
  filter(species == "Prietella_phreatophila") %>%
  mutate(tot_mut = lof_nb * count) %>%
  pull(tot_mut) %>%
  sum()

#Number of splice site mutations

Nb_lof_intron <- 
  nrow(vision_lof_df %>%
  filter(species == "Prietella_phreatophila") %>%
  filter(Genotype %in% c("LoF_Homozygous", "Heterozygous")) %>%
  filter(lof_type == "Splice_site"))


#Number of LoF other than splice site mutations
Nb_lof_cds <- Number_Lof - Nb_lof_intron

#Lets start the simulation

#65, 59, 52, 46, 39 and 33 genes correspond to 1, 0.9, 0.8, 0.7, 0.6 and 0.5 in gene proportion
P_p_vision_sequences_df <- 
  vision_sequences_df_final %>%
  filter(species == "Prietella_phreatophila")

Table_rslt_simu_final <- as.data.frame(NULL)
nb_simu <- 10000
#for(NB_neutral_gene in c(65 ,59, 52, 46, 39, 33)){
#  
#  nb_genes_0_lof <- c()
#  nb_genes_1_lof <- c()
#  nb_genes_2_lof <- c()
#  nb_genes_3_lof <- c()
#  nb_genes_4_lof <- c()
#  nb_genes_5_lof <- c()
#  nb_genes_6_lof <- c()
#  
#  for(i in rep(1, nb_simu)){
#    current_table <- sample_n(P_p_vision_sequences_df, NB_neutral_gene) #pick X random genes evolving as neutral
#    
#    simulation_rslt <- bind_rows(
#      sample_n(current_table, size=Nb_lof_cds, weight=CDS_length, replace = TRUE),
#      sample_n(current_table, size=Nb_lof_intron, weight=exon_Count, replace = TRUE))
#    
#    
#    table_Lof_nb <- simulation_rslt %>% group_by(gene_name) %>% summarise(LoF_nb = n()) %>%
#      group_by(LoF_nb) %>% summarise(Count = n())
#    
#    
#    nb_genes_0_lof <- c(nb_genes_0_lof, Total_Number_Genes - length(unique(simulation_rslt %>% pull(gene_name))))
#    nb_genes_1_lof <- c(nb_genes_1_lof, table_Lof_nb %>% filter(LoF_nb == 1) %>% pull(Count))
#    nb_genes_2_lof <- c(nb_genes_2_lof, table_Lof_nb %>% filter(LoF_nb == 2) %>% pull(Count))
#    nb_genes_3_lof <- c(nb_genes_3_lof, table_Lof_nb %>% filter(LoF_nb == 3) %>% pull(Count))
#    nb_genes_4_lof <- c(nb_genes_4_lof, table_Lof_nb %>% filter(LoF_nb == 4) %>% pull(Count))
#    nb_genes_5_lof <- c(nb_genes_5_lof, sum(table_Lof_nb %>% filter(LoF_nb == 5) %>% pull(Count)))
#    nb_genes_6_lof <- c(nb_genes_6_lof, sum(table_Lof_nb %>% filter(LoF_nb >= 6) %>% pull(Count)))
#    
#  }
#  
#  mean_0_lof <- (sum(nb_genes_0_lof)/nb_simu) * 100 / Total_Number_Genes
#  mean_1_lof <- (sum(nb_genes_1_lof)/nb_simu) * 100 / Total_Number_Genes
#  mean_2_lof <- (sum(nb_genes_2_lof)/nb_simu) * 100 / Total_Number_Genes
#  mean_3_lof <- (sum(nb_genes_3_lof)/nb_simu) * 100 / Total_Number_Genes
#  mean_4_lof <- (sum(nb_genes_4_lof)/nb_simu) * 100 / Total_Number_Genes
#  mean_5_lof <- (sum(nb_genes_5_lof)/nb_simu) * 100 / Total_Number_Genes
#  mean_6_lof <- (sum(nb_genes_6_lof)/nb_simu) * 100 / Total_Number_Genes
#  
#  
#  LoF_count_vector <- c(mean_0_lof, mean_1_lof, mean_2_lof, mean_3_lof, mean_4_lof, mean_5_lof, mean_6_lof)
#  LoF_nb_vector <- c(0, 1, 2, 3, 4, 5, 6)
#  Nb_neutral_gene_vector <- rep(NB_neutral_gene, 7)
#  
#  Table_rslt_simu_ongoing <- as.data.frame(cbind(Nb_neutral_gene_vector, LoF_nb_vector, LoF_count_vector))
#  
#  Table_rslt_simu_final <- rbind(Table_rslt_simu_final, Table_rslt_simu_ongoing)
#  
#}
#
#write.table(Table_rslt_simu_final, "Table_rslt_simu_final.tsv",
#            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

Table_rslt_simu_final <- 
  read.table("Table_rslt_simu_final.tsv",
             sep="\t",
             header=TRUE)



Table_rslt_simu_final %>% filter(Nb_neutral_gene_vector == 65)

## Observed numbers of LoF ##

obs_table <-
  LoF_distrib_obs_df_summary_total_exp %>%
  filter(species == "Prietella_phreatophila") %>%
  arrange(lof_nb)

obs_1_lof <- (obs_table %>% filter(lof_nb == 1) %>% pull(count)) * 100 / Total_Number_Genes
obs_2_lof <- (obs_table %>% filter(lof_nb == 2) %>% pull(count)) * 100 / Total_Number_Genes
obs_3_lof <- (obs_table %>% filter(lof_nb == 3) %>% pull(count)) * 100 / Total_Number_Genes
obs_4_lof <- (obs_table %>% filter(lof_nb == 4) %>% pull(count)) * 100 / Total_Number_Genes
obs_5_lof <- (sum(obs_table %>% filter(lof_nb == 5) %>% pull(count))) * 100 / Total_Number_Genes
obs_6_lof <- (sum(obs_table %>% filter(lof_nb >= 6) %>% pull(count))) * 100 / Total_Number_Genes
obs_0_lof <- (Total_Number_Genes - Total_Number_Pseudo) * 100 / Total_Number_Genes


order_lof_nb_plot <- 
  c(33, 39, 46, 52, 59, 65)


Table_rslt_simu_final$Nb_neutral_gene_vector <- 
  factor(Table_rslt_simu_final$Nb_neutral_gene_vector,
         levels = rev(order_lof_nb_plot))


#Draw the graphic with T. subt only


Table_rslt_simu_final %>% 
  ggplot(., aes(x=LoF_nb_vector, y=LoF_count_vector, group=as.factor(Nb_neutral_gene_vector))) +
  geom_bar(stat="identity", position="dodge", aes(fill= as.factor(Nb_neutral_gene_vector))) +
  theme_minimal() +
  scale_fill_manual(values= c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"),
                    name = "Neutral gene number", labels = c("65" ,"59", "52", "46", "39", "33")) +
  geom_segment(aes(x = -0.45, y = obs_0_lof, xend = 0.5, yend = obs_0_lof, color = "red")) +
  geom_segment(aes(x = 0.55, y = obs_1_lof, xend = 1.4, yend = obs_1_lof, color = "red")) +
  geom_segment(aes(x = 1.55, y = obs_2_lof, xend = 2.4, yend = obs_2_lof, color = "red")) +
  geom_segment(aes(x = 2.55, y = obs_3_lof, xend = 3.4, yend = obs_3_lof, color = "red")) +
  geom_segment(aes(x = 3.55, y = obs_4_lof, xend = 4.4, yend = obs_4_lof, color = "red")) +
  geom_segment(aes(x = 4.55, y = obs_5_lof, xend = 5.4, yend = obs_5_lof, color = "red")) +
  geom_segment(aes(x = 5.55, y = obs_6_lof, xend = 6.4, yend = obs_6_lof, color = "red")) +
  xlab("Number of LoF") +
  ylab("Percentage of genes") 


#### Shared pseudogenes number - L.d x P.p  ---------------------------------

#Extract gene names from L. dentata and Prietella

Dentata_Table <- read.table("Dentata_Table.tsv", header=FALSE, sep="\t")
Ld_genes <- Dentata_Table %>% pull(V2)
Pp_genes <- P_p_vision_sequences_df %>% pull(gene_name)

P_p_vision_sequences_df_ren <- P_p_vision_sequences_df
# Remove clade specific duplications and rename genes to match both species gene names

Pp_genes[Pp_genes == "rh1.1"] <- "rh1"
Pp_genes <- Pp_genes[Pp_genes != "rh1.2"]
Pp_genes[Pp_genes == "rh2.1"] <- "rh2"
Pp_genes[Pp_genes == "gcap7.1"] <- "gcap7"
Pp_genes <- Pp_genes[Pp_genes != "gcap7.2"]
Pp_genes[Pp_genes == "gc2.1"] <- "gc2"
Pp_genes <- Pp_genes[Pp_genes != "gc2.2"]
Ld_genes[Ld_genes == "gngt2b"] <- "gngt2"

Dentata_Table[Dentata_Table$V2 == "gngt2b", "V2"] <- "gngt2"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "rh1.1", "gene_name"] <- "rh1"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "rh2.1", "gene_name"] <- "rh2"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "gcap7.1", "gene_name"] <- "gcap7"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "gc2.1", "gene_name"] <- "gc2"
P_p_vision_sequences_df_ren[55, "gene_name"] <- "gc2"

#Compute the number of genes in common

nb_genes_common <- length(intersect(Ld_genes, Pp_genes))
nb_genes_Ld <- length(Ld_genes)
nb_genes_Pp <- length(Pp_genes)


#Compute the number of pseudogenes in common

Pp_pseudogenes <- 
  P_p_vision_sequences_df %>% 
  filter(corrected_gene_type == "Pseudogene_Homozygous") %>%
  pull(gene_name)

Ld_pseudogenes <- Dentata_Table %>% filter(V6 == "Pseudogene") %>% pull(V2)
  

Pp_pseudogenes[Pp_pseudogenes == "rh1.1"] <- "rh1"
Pp_pseudogenes <- Pp_pseudogenes[Pp_pseudogenes != "rh1.2"]
Pp_pseudogenes[Pp_pseudogenes == "rh2.1"] <- "rh2"
Pp_pseudogenes[Pp_pseudogenes == "gcap7.1"] <- "gcap7"
Pp_pseudogenes <- Pp_pseudogenes[Pp_pseudogenes != "gcap7.2"]
Pp_pseudogenes[Pp_pseudogenes == "gc2.1"] <- "gc2"
Pp_pseudogenes <- Pp_pseudogenes[Pp_pseudogenes != "gc2.2"]
Ld_pseudogenes[Ld_pseudogenes == "gngt2b"] <- "gngt2"

nb_pseudo_common <- length(intersect(Ld_pseudogenes, Pp_pseudogenes))
pseudo_in_common <- intersect(Ld_pseudogenes, Pp_pseudogenes)

  
  
nb_pseudo_Ld <- length(Ld_pseudogenes)
nb_pseudo_Pp <- length(Pp_pseudogenes)


#Expected number of pseudogenes in common:
(nb_pseudo_Pp/nb_genes_Pp) * (nb_pseudo_Ld/nb_genes_Ld) * nb_genes_common


#List of common genes
genes_in_common <- intersect(Ld_genes, Pp_genes)

#Number of pseudogenes among common genes in P.p 
nb_subset_pseudo_Pp <- 
  nrow(P_p_vision_sequences_df_ren %>% filter(gene_name %in% genes_in_common) %>% 
         filter(corrected_gene_type == "Pseudogene_Homozygous"))

#Number of pseudogenes among common genes in L.dent
nb_subset_pseudo_Ld <- 
  nrow(Dentata_Table %>% filter(V2 %in% genes_in_common) %>% filter(V6 == "Pseudogene"))

### Lets start simulations taking into account genes lengths and the number of introns

#Compute mean length of genes


t1 <- P_p_vision_sequences_df_ren %>% 
  filter(gene_name %in% genes_in_common) %>% 
  dplyr::select(gene_name, CDS_length, exon_Count)

t2 <- Dentata_Table %>% 
  filter(V2 %in% genes_in_common) %>% 
  dplyr::select(V2, V5, V8) 
colnames(t2) <- colnames(t1)

t1$CDS_length <- as.numeric(t1$CDS_length)
t2$CDS_length <- as.numeric(t2$CDS_length)
t1$exon_Count <- as.numeric(t1$exon_Count)
t2$exon_Count <- as.numeric(t2$exon_Count)

t1 <- t1 %>% mutate(intron_Count = exon_Count-1)
t2 <- t2 %>% mutate(intron_Count = exon_Count-1)



mean_CDS_length_df <- 
  as.data.frame(left_join(t1, t2, by="gene_name") %>%
  rowwise() %>%
  mutate(mean_CDS_length = mean(c(CDS_length.x, CDS_length.y))) %>%
    mutate(mean_intron_nb = mean(c(intron_Count.x, intron_Count.y))) %>%
  dplyr::select(gene_name, mean_CDS_length, mean_intron_nb))




#Normalize cds length and intron number by the LoF proba rate in exon or in intron 


mean_CDS_length_df <- mean_CDS_length_df %>%
  mutate(exon_weight = (proba_stop_gain + proba_frameshift + proba_stop_loss + proba_start_loss) * mean_CDS_length) %>%
  mutate(intron_weight = proba_splice * mean_intron_nb) %>%
  mutate(Proba_mut = intron_weight+exon_weight)


#Lets start the simulations
rslt_simu_inter <- c()
neutal_vector_sim <- c()
#for(neutral_portion in seq(nb_subset_pseudo_Pp, nrow(mean_CDS_length_df), 1)){
#  for(i in rep(1, 10000)){
#    
#    subset_table_picked <- sample_n(mean_CDS_length_df, neutral_portion)
#    
#    nb_intersect <- 
#      length(
#        intersect(
#          sample_n(subset_table_picked, nb_subset_pseudo_Ld, weight=Proba_mut) %>% pull(gene_name),
#          sample_n(subset_table_picked, nb_subset_pseudo_Pp, weight=Proba_mut) %>% pull(gene_name)
#          )
#        )
#    
#    rslt_simu_inter <- c(rslt_simu_inter, nb_intersect)
#    neutal_vector_sim <- c(neutal_vector_sim, neutral_portion)
#    
#  }
#}
#
#Simu_df_exon_intron <- as.data.frame(cbind(rslt_simu_inter, neutal_vector_sim))
#
#write.table(Simu_df_exon_intron, "Simu_df_exon_intron.tsv",
#            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

Simu_df_exon_intron <- 
  read.table("Simu_df_exon_intron.tsv", header=TRUE, sep="\t")

test_exon_intron <- Simu_df_exon_intron %>% group_by(neutal_vector_sim) %>%
  summarise(n_10_obs = (sum(rslt_simu_inter == nb_pseudo_common)),
            n_simu = n(),
            p_value = n_10_obs / 10000) 



as.data.frame(
  Simu_df_exon_intron %>%
  group_by(neutal_vector_sim) %>%
  summarise(mean_common = mean(rslt_simu_inter))
)



Simu_df_exon_intron %>% filter(neutal_vector_sim == 30) %>% filter(rslt_simu_inter == 9)
Simu_df_exon_intron %>% filter(neutal_vector_sim == 48) %>% filter(rslt_simu_inter == 9)

Simu_df_exon_intron %>% filter(neutal_vector_sim == nb_genes_common) %>%
  group_by(rslt_simu_inter) %>%
  summarise(nb_obs = n()) %>%
  ggplot(., aes(x=rslt_simu_inter, y=nb_obs)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_x_continuous("Number of common pseudogenes", c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)) +
  ylab("Number of observations") 


Simu_df_exon_intron %>% filter(neutal_vector_sim == nb_genes_common) %>%
  group_by(rslt_simu_inter) %>%
  summarise(nb_obs = n()) %>%
  mutate(nb_obs_prop = nb_obs/10000) %>%
  ggplot(., aes(x=rslt_simu_inter, y=nb_obs_prop)) +
  geom_bar(stat="identity", color="black", fill="gray") +
  theme_minimal() +
  scale_x_continuous("Number of common pseudogenes", c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)) +
  ylab("Proportion of observations") #+
  #geom_hline(yintercept = 0.05, col="red", linetype="dashed")


Simu_df_exon_intron %>% filter(neutal_vector_sim == nb_genes_common) %>% pull(rslt_simu_inter) %>% mean(.)


test_exon_intron %>% 
  ggplot(., aes(x=neutal_vector_sim, y=p_value)) +
  geom_point() + 
  theme_minimal() +
  xlab("Number of neutral genes") +
  ylab("Probability of finding 9 common pseudogenes") #+



#### Shared pseudogenes number - A.r x P.p  ---------------------------------

#Extract gene names from A. rosae and Prietella


Perco_vision_df <- read.table("Perco_Table.tsv", header=FALSE, sep="\t")
colnames(Perco_vision_df) <- 
  c("Gene_clade","Species","gene_name","Genomic_position","CDS_length","Gene_Type","Note",
    "Exon_Count","Sequence")

Tr_vision_df <- Perco_vision_df %>% filter(Species == "Amblyopsis_rosae")
 
Tr_genes <- Tr_vision_df %>% pull(gene_name)
Pp_genes <- P_p_vision_sequences_df %>% pull(gene_name)

P_p_vision_sequences_df_ren <- P_p_vision_sequences_df




# Remove clade specific duplications and rename genes to match both species gene names

P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "exorh", "gene_name"] <- "Exorh"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "lws1", "gene_name"] <- "Lws1"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "opn3", "gene_name"] <- "Opn3"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "opn4m3", "gene_name"] <- "Opn4m3"
P_p_vision_sequences_df_ren[55, "gene_name"] <- "gc2"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "opn7c", "gene_name"] <- "Opn7c"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "parapinopsin-1", "gene_name"] <- "Parapinopsin-1"
P_p_vision_sequences_df_ren$gene_name <- gsub("cry", "Cry", P_p_vision_sequences_df_ren$gene_name)
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "rrh", "gene_name"] <- "Rrh"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "rh1.1", "gene_name"] <- "Rh1.1"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "rh2.1", "gene_name"] <- "Rh2.1"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "tmt1b", "gene_name"] <- "Tmt1b"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "tmt3a", "gene_name"] <- "Tmt3a"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "va2", "gene_name"] <- "Va2"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "rpe65a", "gene_name"] <- "Rpe65a"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "gngt2", "gene_name"] <- "gngt2b"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "parietopsin", "gene_name"] <- "Parietopsin"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "rgr1", "gene_name"] <- "Rgr1"
P_p_vision_sequences_df_ren[P_p_vision_sequences_df_ren$gene_name == "rgr2", "gene_name"] <- "Rgr2"

Pp_genes <- P_p_vision_sequences_df_ren %>% pull(gene_name)



#Compute the number of genes in common

nb_genes_common <- length(intersect(Tr_genes, Pp_genes))
nb_genes_Tr <- length(Tr_genes)
nb_genes_Pp <- length(Pp_genes)


#Compute the number of pseudogenes in common

Pp_pseudogenes <- 
  P_p_vision_sequences_df_ren %>% 
  filter(corrected_gene_type == "Pseudogene_Homozygous") %>%
  pull(gene_name)

Tr_pseudogenes <- Tr_vision_df %>% filter(Gene_Type == "Pseudogene") %>% pull(gene_name)


nb_pseudo_common <- length(intersect(Tr_pseudogenes, Pp_pseudogenes))
pseudo_in_common <- intersect(Tr_pseudogenes, Pp_pseudogenes)



nb_pseudo_Tr <- length(Tr_pseudogenes)
nb_pseudo_Pp <- length(Pp_pseudogenes)


#Expected number of pseudogenes in common:
(nb_pseudo_Pp/nb_genes_Pp) * (nb_pseudo_Tr/nb_genes_Tr) * nb_genes_common


#List of common genes
genes_in_common <- intersect(Tr_genes, Pp_genes)

#Number of pseudogenes among common genes in P.p 
nb_subset_pseudo_Pp <- 
  nrow(P_p_vision_sequences_df_ren %>% filter(gene_name %in% genes_in_common) %>% 
         filter(corrected_gene_type == "Pseudogene_Homozygous"))

#Number of pseudogenes among common genes in L.dent
nb_subset_pseudo_Tr <- 
  nrow(Tr_vision_df %>% filter(gene_name %in% genes_in_common) %>% filter(Gene_Type == "Pseudogene"))

### Lets start simulations taking into account genes lengths and the number of introns

#Compute mean length of genes

t1 <- P_p_vision_sequences_df_ren %>% 
  filter(gene_name %in% genes_in_common) %>% 
  dplyr::select(gene_name, CDS_length, exon_Count)

t2 <- Tr_vision_df %>% 
  filter(gene_name %in% genes_in_common) %>%
  dplyr::select(gene_name, CDS_length, Exon_Count) 
colnames(t2) <- colnames(t1)

t1$CDS_length <- as.numeric(t1$CDS_length)
t2$CDS_length <- as.numeric(t2$CDS_length)
t1$exon_Count <- as.numeric(t1$exon_Count)
t2$exon_Count <- as.numeric(t2$exon_Count)

t1 <- t1 %>% mutate(intron_Count = exon_Count-1)
t2 <- t2 %>% mutate(intron_Count = exon_Count-1)



mean_CDS_length_df <- 
  as.data.frame(left_join(t1, t2, by="gene_name") %>%
                  rowwise() %>%
                  mutate(mean_CDS_length = mean(c(CDS_length.x, CDS_length.y))) %>%
                  mutate(mean_intron_nb = mean(c(intron_Count.x, intron_Count.y))) %>%
                  dplyr::select(gene_name, mean_CDS_length, mean_intron_nb))




#Normalize cds length and intron number by the LoF proba rate in exon or in intron 


mean_CDS_length_df <- mean_CDS_length_df %>%
  mutate(exon_weight = (proba_stop_gain + proba_frameshift + proba_stop_loss + proba_start_loss) * mean_CDS_length) %>%
  mutate(intron_weight = proba_splice * mean_intron_nb) %>%
  mutate(Proba_mut = intron_weight+exon_weight)

#
##Lets start the simulations
#rslt_simu_inter <- c()
#neutal_vector_sim <- c()
#for(neutral_portion in seq(nb_subset_pseudo_Tr, nrow(mean_CDS_length_df), 1)){
#  for(i in rep(1, 10000)){
#    
#    subset_table_picked <- sample_n(mean_CDS_length_df, neutral_portion)
#    
#    nb_intersect <- 
#      length(
#        intersect(
#          sample_n(subset_table_picked, nb_subset_pseudo_Tr, weight=Proba_mut) %>% pull(gene_name),
#          sample_n(subset_table_picked, nb_subset_pseudo_Pp, weight=Proba_mut) %>% pull(gene_name)
#          )
#        )
#    
#    rslt_simu_inter <- c(rslt_simu_inter, nb_intersect)
#    neutal_vector_sim <- c(neutal_vector_sim, neutral_portion)
#    
#  }
#}
#
#Simu_df_exon_intron <- as.data.frame(cbind(rslt_simu_inter, neutal_vector_sim))
#
#write.table(Simu_df_exon_intron, "Simu_df_exon_intron.rosae.tsv",
#            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

Simu_df_exon_intron <- 
  read.table("Simu_df_exon_intron.rosae.tsv", header=TRUE, sep="\t")

test_exon_intron <- Simu_df_exon_intron %>% group_by(neutal_vector_sim) %>%
  summarise(n_10_obs = (sum(rslt_simu_inter == nb_pseudo_common)),
            n_simu = n(),
            p_value = n_10_obs / 10000) 



as.data.frame(
  Simu_df_exon_intron %>%
    group_by(neutal_vector_sim) %>%
    summarise(mean_common = mean(rslt_simu_inter))
)



Simu_df_exon_intron %>% filter(neutal_vector_sim == nb_genes_common) %>%
  group_by(rslt_simu_inter) %>%
  summarise(nb_obs = n()) %>%
  ggplot(., aes(x=rslt_simu_inter, y=nb_obs)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_x_continuous("Number of common pseudogenes", c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)) +
  ylab("Number of observations") 


Simu_df_exon_intron %>% filter(neutal_vector_sim == nb_genes_common) %>%
  group_by(rslt_simu_inter) %>%
  summarise(nb_obs = n()) %>%
  mutate(nb_obs_prop = nb_obs/10000) %>%
  ggplot(., aes(x=rslt_simu_inter, y=nb_obs_prop)) +
  geom_bar(stat="identity", color="black", fill="gray") +
  theme_minimal() +
  scale_x_continuous("Number of common pseudogenes", c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)) +
  ylab("Proportion of observations") #+


Simu_df_exon_intron %>% filter(neutal_vector_sim == nb_genes_common) %>% pull(rslt_simu_inter) %>% mean(.)



test_exon_intron %>% 
  ggplot(., aes(x=neutal_vector_sim, y=p_value)) +
  geom_point() + 
  theme_minimal() +
  xlab("Number of neutral genes") +
  ylab("Probability of finding 7 common pseudogenes") #+


Simu_df_exon_intron %>% 
  filter(neutal_vector_sim == 28) %>% 
  pull(rslt_simu_inter) %>% mean(.)


#### Load concatenated dN/dS  -- Meredith - noClean  ---------------------------------


neutral.dNdS <- 1

## Prietella_phreatophila

Prietella_dNdS_bs <- 
  read.table("Meredith_method/Prietella_phreatophila/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(Prietella_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

DivergenceTime_P_A <- 12.9362 * 1000000


Prietella_dNdS_bs <- 
  Prietella_dNdS_bs %>%
  mutate(Tp = 
           DivergenceTime_P_A - ((DivergenceTime_P_A * (Cave_dNdS - neutral.dNdS)) / 
                                   (Surface_dNdS - neutral.dNdS))) 
  

Prietella_dNdS_bs <- 
  Prietella_dNdS_bs %>% 
  mutate(Gp = Tp/3)


## CUL9

CUL9_dNdS_bs <- 
  read.table("Meredith_method/CUL9/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(CUL9_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

DivergenceTime_CUL9_RHP1 <- 8.188 * 1000000

CUL9_dNdS_bs <- 
  CUL9_dNdS_bs %>%
  mutate(Tp = 
           DivergenceTime_CUL9_RHP1 - ((DivergenceTime_CUL9_RHP1 * (Cave_dNdS - neutral.dNdS)) / 
                                   (Surface_dNdS - neutral.dNdS))) 
CUL9_dNdS_bs <- 
  CUL9_dNdS_bs %>% 
  mutate(Gp = Tp/3)
  


## CSV83

CSV83_dNdS_bs <- 
  read.table("Meredith_method/CSV83/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(CSV83_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

DivergenceTime_CSVandCUL4_RUI2 <- 14.6389 * 1000000

CSV83_dNdS_bs <- 
  CSV83_dNdS_bs %>%
  mutate(Tp = 
           DivergenceTime_CSVandCUL4_RUI2 - ((DivergenceTime_CSVandCUL4_RUI2 * (Cave_dNdS - neutral.dNdS)) / 
                                               (Surface_dNdS - neutral.dNdS))) 
CSV83_dNdS_bs <- 
  CSV83_dNdS_bs %>% 
  mutate(Gp = Tp/3)


## CUL4

CUL4_dNdS_bs <- 
  read.table("Meredith_method/CUL4/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(CUL4_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

DivergenceTime_CSVandCUL4_RUI2 <- 14.6389 * 1000000

CUL4_dNdS_bs <- 
  CUL4_dNdS_bs %>%
  mutate(Tp = 
           DivergenceTime_CSVandCUL4_RUI2 - ((DivergenceTime_CSVandCUL4_RUI2 * (Cave_dNdS - neutral.dNdS)) / 
                                               (Surface_dNdS - neutral.dNdS))) 
CUL4_dNdS_bs <- 
  CUL4_dNdS_bs %>% 
  mutate(Gp = Tp/3)



## CHM6

CHM6_dNdS_bs <- 
  read.table("Meredith_method/CHM6/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(CHM6_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

DivergenceTime_CHM6_RSS1 <- 10.5187 * 1000000

CHM6_dNdS_bs <- 
  CHM6_dNdS_bs %>%
  mutate(Tp = 
           DivergenceTime_CHM6_RSS1 - ((DivergenceTime_CHM6_RSS1 * (Cave_dNdS - neutral.dNdS)) / 
                                         (Surface_dNdS - neutral.dNdS))) 
CHM6_dNdS_bs <- 
  CHM6_dNdS_bs %>% 
  mutate(Gp = Tp/3)



## Trichomycterus_rosablanca

Trichomycterus_rosablanca_dNdS_bs <- 
  read.table("Meredith_method/Trichomycterus_rosablanca/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(Trichomycterus_rosablanca_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

DivergenceTime_T_RSS1 <- 15.1434 * 1000000

Trichomycterus_rosablanca_dNdS_bs <- 
  Trichomycterus_rosablanca_dNdS_bs %>%
  mutate(Tp = 
           DivergenceTime_T_RSS1 - ((DivergenceTime_T_RSS1 * (Cave_dNdS - neutral.dNdS)) / 
                                      (Surface_dNdS - neutral.dNdS))) 
Trichomycterus_rosablanca_dNdS_bs <- 
  Trichomycterus_rosablanca_dNdS_bs %>% 
  mutate(Gp = Tp/3)






## Summary 

Summary_dNdS_datation_noclean <- 
  rbind(
    Prietella_dNdS_bs %>% mutate(species = "Prietella_phreatophila"),
    CUL9_dNdS_bs %>% mutate(species = "CUL9"),
    CSV83_dNdS_bs %>% mutate(species = "CSV83"),
    CUL4_dNdS_bs %>% mutate(species = "CUL4"),
    CHM6_dNdS_bs %>% mutate(species = "CHM6"),
    Trichomycterus_rosablanca_dNdS_bs %>% mutate(species = "Trichomycterus_rosablanca")
  )


Summary_dNdS_datation_noclean %>% filter(bootstrap_nb == "original_aln")


Summary_dNdS_datation_noclean %>%
  filter(bootstrap_nb != "original_aln") %>%
  filter(species == "Trichomycterus_rosablanca") %>%
  summarise(
    lower_95  = quantile(Tp, 0.025, na.rm = TRUE),
    upper_95  = quantile(Tp, 0.975, na.rm = TRUE)
  )


Summary_dNdS_datation_noclean %>%
  ggplot(., aes(x=Gp, color=species)) +
  geom_density() +
  theme_classic() +
  scale_color_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9"
  )) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  xlim(0, 1700000) +
  geom_vline(xintercept = (11700/3), color= "gray") +
  geom_vline(xintercept = (126000/3), color= "gray") +
  geom_vline(xintercept = (781000/3), color= "gray") +
  geom_vline(xintercept = (2580000/3), color= "gray") +
  geom_vline(xintercept = (5333000/3), color= "gray") 


# geom_vline(xintercept = Prietella_dNdS_bs %>% filter(bootstrap_nb == "original_aln") %>% pull(Gp)) +

dens <- density(Summary_dNdS_datation_noclean %>% filter(species == "Prietella_phreatophila") %>% pull(Gp))
max(dens$y)


Summary_dNdS_datation_noclean %>%
  group_by(species) %>%
  summarise(mean_Tp = mean(Tp))


#Look at quantile values

Summary_dNdS_datation_noclean %>%
  filter(species == "Prietella_phreatophila") %>%
  pull(Tp) %>% quantile(0.05)

Summary_dNdS_datation_noclean %>%
  filter(species == "Prietella_phreatophila") %>%
  pull(Tp) %>% quantile(0.95)

#### Load concatenated dN/dS  -- Meredith - Clean  ---------------------------------


neutral.dNdS <- 1

## Prietella_phreatophila

Prietella_dNdS_bs <- 
  read.table("Meredith_method_CleanData/Prietella_phreatophila/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(Prietella_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

DivergenceTime_P_A <- 12.9362 * 1000000


Prietella_dNdS_bs <- 
  Prietella_dNdS_bs %>%
  mutate(Tp = 
           DivergenceTime_P_A - ((DivergenceTime_P_A * (Cave_dNdS - neutral.dNdS)) / 
                                   (Surface_dNdS - neutral.dNdS))) 


Prietella_dNdS_bs <- 
  Prietella_dNdS_bs %>% 
  mutate(Gp = Tp/3)


## CUL9

CUL9_dNdS_bs <- 
  read.table("Meredith_method_CleanData/CUL9/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(CUL9_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

DivergenceTime_CUL9_RHP1 <- 8.188 * 1000000

CUL9_dNdS_bs <- 
  CUL9_dNdS_bs %>%
  mutate(Tp = 
           DivergenceTime_CUL9_RHP1 - ((DivergenceTime_CUL9_RHP1 * (Cave_dNdS - neutral.dNdS)) / 
                                         (Surface_dNdS - neutral.dNdS))) 
CUL9_dNdS_bs <- 
  CUL9_dNdS_bs %>% 
  mutate(Gp = Tp/3)



## CSV83

CSV83_dNdS_bs <- 
  read.table("Meredith_method_CleanData/CSV83/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(CSV83_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

DivergenceTime_CSVandCUL4_RUI2 <- 14.6389 * 1000000

CSV83_dNdS_bs <- 
  CSV83_dNdS_bs %>%
  mutate(Tp = 
           DivergenceTime_CSVandCUL4_RUI2 - ((DivergenceTime_CSVandCUL4_RUI2 * (Cave_dNdS - neutral.dNdS)) / 
                                               (Surface_dNdS - neutral.dNdS))) 
CSV83_dNdS_bs <- 
  CSV83_dNdS_bs %>% 
  mutate(Gp = Tp/3)


## CUL4

CUL4_dNdS_bs <- 
  read.table("Meredith_method_CleanData/CUL4/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(CUL4_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

DivergenceTime_CSVandCUL4_RUI2 <- 14.6389 * 1000000

CUL4_dNdS_bs <- 
  CUL4_dNdS_bs %>%
  mutate(Tp = 
           DivergenceTime_CSVandCUL4_RUI2 - ((DivergenceTime_CSVandCUL4_RUI2 * (Cave_dNdS - neutral.dNdS)) / 
                                               (Surface_dNdS - neutral.dNdS))) 
CUL4_dNdS_bs <- 
  CUL4_dNdS_bs %>% 
  mutate(Gp = Tp/3)



## CHM6

CHM6_dNdS_bs <- 
  read.table("Meredith_method_CleanData/CHM6/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(CHM6_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

DivergenceTime_CHM6_RSS1 <- 10.5187 * 1000000

CHM6_dNdS_bs <- 
  CHM6_dNdS_bs %>%
  mutate(Tp = 
           DivergenceTime_CHM6_RSS1 - ((DivergenceTime_CHM6_RSS1 * (Cave_dNdS - neutral.dNdS)) / 
                                         (Surface_dNdS - neutral.dNdS))) 
CHM6_dNdS_bs <- 
  CHM6_dNdS_bs %>% 
  mutate(Gp = Tp/3)



## Trichomycterus_rosablanca

Trichomycterus_rosablanca_dNdS_bs <- 
  read.table("Meredith_method_CleanData/Trichomycterus_rosablanca/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(Trichomycterus_rosablanca_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

DivergenceTime_T_RSS1 <- 15.1434 * 1000000

Trichomycterus_rosablanca_dNdS_bs <- 
  Trichomycterus_rosablanca_dNdS_bs %>%
  mutate(Tp = 
           DivergenceTime_T_RSS1 - ((DivergenceTime_T_RSS1 * (Cave_dNdS - neutral.dNdS)) / 
                                      (Surface_dNdS - neutral.dNdS))) 
Trichomycterus_rosablanca_dNdS_bs <- 
  Trichomycterus_rosablanca_dNdS_bs %>% 
  mutate(Gp = Tp/3)


## Summary 

Summary_dNdS_datation_clean <- 
  rbind(
    Prietella_dNdS_bs %>% mutate(species = "Prietella_phreatophila"),
    CUL9_dNdS_bs %>% mutate(species = "CUL9"),
    CSV83_dNdS_bs %>% mutate(species = "CSV83"),
    CUL4_dNdS_bs %>% mutate(species = "CUL4"),
    CHM6_dNdS_bs %>% mutate(species = "CHM6"),
    Trichomycterus_rosablanca_dNdS_bs %>% mutate(species = "Trichomycterus_rosablanca")
  )


Summary_dNdS_datation_clean %>% filter(bootstrap_nb == "original_aln")


#### Prepare values for datations using the number of pseudogenes  ---------------------------------

species_to_date <- 
  species_df %>%
  filter(Habitat == "Cave") %>%
  pull(species)


#Extract the mean length of genes per species


gene_length_df <- 
  as.data.frame(
    vision_sequences_df_final %>%
      filter(species %in% species_to_date) %>%
      group_by(species) %>%
      summarise(mean_CDS_length = mean(CDS_length))
  )


summary_vision_df_wide_date <-
  left_join(summary_vision_df_wide, gene_length_df, by="species") %>%
  filter(! is.na(mean_CDS_length)) %>%
  dplyr::select(species, Total, Complete, Incomplete, Pseudogene, mean_CDS_length)


#remove genes/pseudogenes found pseudogenised in a internal branch

summary_vision_df_wide_date[(summary_vision_df_wide_date$species == "CUL9"),"Total"] <- 64
summary_vision_df_wide_date[(summary_vision_df_wide_date$species == "CUL9"),"Pseudogene"] <- 7

#Low nb = nb pseudo - nb pseudo of the closest surface species
#Medium nb = nb pseudo
#High nb = nb pseudo + nb of genes missing compared to the closest surface species

summary_vision_df_wide_date <- 
  summary_vision_df_wide_date %>%
  mutate(ref_surface_sp = case_when(
    species == "Prietella_phreatophila" ~ "Ameiurus_melas",
    species %in% c("CUL4", "CSV83") ~ "RUI2",
    species == "CHM6" ~ "RHH2",
    species == "CUL9" ~ "RHP1",
    species == "Trichomycterus_rosablanca" ~ "RHH2"
  ))


#### Perform datation (pseudogene number) of cavefishes - mutation rate = 10E-8  ---------------------------------

species_to_date <- 
  summary_vision_df_wide_date$species

Datation_results_df_normal <- as.data.frame(NULL)
curr_LoF_rate <- LoF_rates_df %>% pull(Total)
for(curr_species in species_to_date){
  
  curr_gene_number_normal <- summary_vision_df_wide_date %>% filter(species == curr_species) %>% pull(Total)
  curr_Dn_normal <- summary_vision_df_wide_date %>% filter(species == curr_species) %>% pull(Pseudogene)  
  curr_meanCDSlength <- summary_vision_df_wide_date %>% filter(species == curr_species) %>% pull(mean_CDS_length)  
  
  
  
  curr_Datation_results_df <- as.data.frame(seq(1, 3000000, 10))
  colnames(curr_Datation_results_df) <- c("curr_time")
  
  if(curr_Dn_normal > 0){
    curr_Datation_results_df_normal <- 
      curr_Datation_results_df %>%
      mutate(gene_number = curr_gene_number_normal) %>%
      mutate(meanCDSlength = curr_meanCDSlength) %>%
      mutate(LoF_rate = as.numeric(curr_LoF_rate)) %>%
      mutate(Dn = curr_Dn_normal)
    
    curr_Datation_results_df_normal <- 
      as.data.frame(
        curr_Datation_results_df_normal %>%
          rowwise() %>%
          mutate(
            curr_proba=(((factorial(gene_number))/((factorial(Dn)*(factorial(gene_number-Dn))))))*((1-exp((-10^-8)*LoF_rate*curr_time*meanCDSlength))^Dn)*(exp((-10^-8)*LoF_rate*curr_time*meanCDSlength))^(gene_number-Dn))
      )
    
    curr_Datation_results_df_normal <- 
      curr_Datation_results_df_normal %>% mutate(species = curr_species)
    
    colnames(curr_Datation_results_df_normal) <-
      c("GenerationNb","GeneNb","mean_cds_length", "LoF_rate", "PseudoNb","probability", "species")
    
    Datation_results_df_normal <- 
      rbind(Datation_results_df_normal, curr_Datation_results_df_normal)
    
    
  }
  
}

Datation_results_df_normal$probability <- as.numeric(Datation_results_df_normal$probability)


#Draw the results


Datation_results_df_normal %>%
  ggplot(., aes(x=GenerationNb, y=probability, color=species)) +
  geom_line(size=1.2) +
  theme_classic() +
  scale_color_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9"
  )) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  xlim(0, 1300000) +
  ylim(0, 0.4) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color= "gray") +
  geom_vline(xintercept = (11700/3), color= "gray") +
  geom_vline(xintercept = (126000/3), color= "gray") +
  geom_vline(xintercept = (781000/3), color= "gray") +
  geom_vline(xintercept = (2580000/3), color= "gray") +
  geom_vline(xintercept = (5333000/3), color= "gray") 
  
  


#Get the max proba for each species

MaxProba_gen <- 
  as.data.frame(
    Datation_results_df_normal %>%
      group_by(species) %>%
      slice_max(probability, n=1) %>%
      arrange(species) %>%
      dplyr::select(GenerationNb, species, probability)
  )
colnames(MaxProba_gen) <- c("Best_GenNb", "species", "MaxProba")


#Get the 0.05 confidence intervals for each species

LowProba05_gen <- 
  as.data.frame(
    Datation_results_df_normal %>%
      filter(probability > 0.05) %>%
      group_by(species) %>%
      slice_min(GenerationNb, n = 1) %>%
      arrange(species) %>%
      dplyr::select(GenerationNb, species)
  )
colnames(LowProba05_gen) <- c("MinGen_p_05", "species")

HighProba05_gen <- 
  as.data.frame(
    Datation_results_df_normal %>%
      filter(probability > 0.05) %>%
      group_by(species) %>%
      slice_max(GenerationNb, n = 1) %>%
      arrange(species) %>%
      dplyr::select(GenerationNb, species)
  )
colnames(HighProba05_gen) <- c("MaxGen_p_05", "species")

Summary_datation <- left_join(MaxProba_gen, LowProba05_gen, by=c("species"))
Summary_datation <- left_join(Summary_datation, HighProba05_gen, by=c("species"))


gt(Summary_datation %>% arrange(Best_GenNb))


# Check how much Mya it would be (gen time between 1 and 5)

Summary_datation <- 
  Summary_datation %>% 
  mutate(MinDate_p_05_1yr = MinGen_p_05 * 1) %>%
  mutate(MinDate_p_05_5yr = MinGen_p_05 * 5) %>%
  mutate(MaxDate_p_05_1yr = MaxGen_p_05 * 1) %>%
  mutate(MaxDate_p_05_5yr = MaxGen_p_05 * 5) %>%
  mutate(Best_GenNb_1yr = Best_GenNb * 1)  %>%
  mutate(Best_GenNb_5yr = Best_GenNb * 5) 



#### Plot datations results on the same graphic  ---------------------------------

## Plot the pseudogene nb in generations


Datation_results_df_normal %>%
  filter(probability > 1e-5) %>%
  ggplot(aes(x = GenerationNb, y = probability, color = species)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9"
  )) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.subtitle = element_text(size = 16),
    legend.position = "none"
  ) +
  #xlim(0, 1500000) +
  ylim(0, 0.4) +
  #geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray") +
  #geom_vline(xintercept = c(11700, 126000, 781000, 2580000, 5333000) / 3, color = "gray") +
  ylab("Probability") +
  xlab("Generations") +
  scale_x_continuous(
    limits = c(0, 1500000),
    breaks = c(0, 100000, 200000, 300000, 400000, 500000, 750000, 1000000, 1250000, 1500000)
  ) 


## Plot the meredith data in years


Summary_dNdS_datation_noclean %>%
  ggplot(., aes(x = Tp, color = species)) +
  geom_density(size = 1.2)  + 
  scale_color_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9"
  )) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.subtitle = element_text(size = 16),
    legend.position = "none"
  ) +
  #xlim(-200000, 6000000) +
  geom_vline(xintercept = c(11700, 126000, 781000, 2580000, 5333000), color = "gray") +
  xlab("Years") +
  scale_x_continuous(
    limits = c(-200000, 6000000),
    breaks = c(0, 500000, 1000000, 1500000, 2000000, 2500000, 
               3000000, 3500000, 4000000, 4500000,
               5000000, 5500000, 6000000)
  ) 


## Combine the two

p1 <- Datation_results_df_normal %>%
  filter(probability > 1e-5) %>%
  ggplot(aes(x = GenerationNb, y = probability, color = species)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9"
  )) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.x = element_blank(),
    plot.subtitle = element_text(size = 16),
    legend.position = "none"
  ) +
  xlim(0, 1500000) +
  ylim(0, 0.4) +
  #geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(11700, 126000, 781000, 2580000, 5333000) / 3, color = "gray") +
  ylab("Probability")



# Combine both plots -- No clean

combined_plot_noclean  <- p1 +
  geom_density(
    data = Summary_dNdS_datation_noclean,
    aes(x = Gp, y = ..scaled.. / 5, color = species),
    inherit.aes = FALSE,
    linetype="dashed", size=0.8
  ) +
  scale_y_continuous(
    name = "Probability",
    sec.axis = sec_axis(~ . * 10, name = "Density")
  )



# Combine both plots -- Clean

combined_plot_clean <- p1 +
  geom_density(
    data = Summary_dNdS_datation_clean,
    aes(x = Gp, y = ..scaled.. / 5, color = species),
    inherit.aes = FALSE,
    linetype="dashed", size=0.8
  ) +
  scale_y_continuous(
    name = "Probability",
    sec.axis = sec_axis(~ . * 10, name = "Density")
  )




#### Draw nexus tree with datations ---------------------------------

library("deeptime")
library(phyloch)
library(EvoPhylo)

dNdS_datation_df <- Summary_dNdS_datation_noclean

beast_species_tree <- treeio::read.beast("AMAS_concatenated_alignment_BUSCO.fa.timetree.nex")


branch_node_prietella <- beast_species_tree@phylo$tip.label == "Prietella_phreatophila"
index_prietella <- Ntip(beast_species_tree) - which(branch_node_prietella)  + 2
dNdS_max_prietella <- as.numeric(dNdS_datation_df %>% filter(species == "Prietella_phreatophila") %>% pull(Tp) %>% quantile(0.95))
dNdS_min_prietella <- as.numeric(dNdS_datation_df %>% filter(species == "Prietella_phreatophila") %>% pull(Tp) %>% quantile(0.05))


branch_node_CHM6 <- beast_species_tree@phylo$tip.label == "CHM6"
index_CHM6 <- Ntip(beast_species_tree) - which(branch_node_CHM6)  + 3
dNdS_max_CHM6 <- as.numeric(dNdS_datation_df %>% filter(species == "CHM6") %>% pull(Tp) %>% quantile(0.95))
dNdS_min_CHM6 <- as.numeric(dNdS_datation_df %>% filter(species == "CHM6") %>% pull(Tp) %>% quantile(0.05))


branch_node_Trichomycterus <- beast_species_tree@phylo$tip.label == "Trichomycterus_rosablanca"
index_Trichomycterus <- Ntip(beast_species_tree) - which(branch_node_Trichomycterus)  + 5
dNdS_max_Trichomycterus <- as.numeric(dNdS_datation_df %>% filter(species == "Trichomycterus_rosablanca") %>% pull(Tp) %>% quantile(0.95))
dNdS_min_Trichomycterus <- as.numeric(dNdS_datation_df %>% filter(species == "Trichomycterus_rosablanca") %>% pull(Tp) %>% quantile(0.05))

branch_node_CSV83 <- beast_species_tree@phylo$tip.label == "CSV83"
index_CSV83 <- Ntip(beast_species_tree) - which(branch_node_CSV83)  + 5
dNdS_max_CSV83 <- as.numeric(dNdS_datation_df %>% filter(species == "CSV83") %>% pull(Tp) %>% quantile(0.95))
dNdS_min_CSV83 <- as.numeric(dNdS_datation_df %>% filter(species == "CSV83") %>% pull(Tp) %>% quantile(0.05))


branch_node_CUL4 <- beast_species_tree@phylo$tip.label == "CUL4"
index_CUL4 <- Ntip(beast_species_tree) - which(branch_node_CUL4)  + 5
dNdS_max_CUL4 <- as.numeric(dNdS_datation_df %>% filter(species == "CUL4") %>% pull(Tp) %>% quantile(0.95))
dNdS_min_CUL4 <- as.numeric(dNdS_datation_df %>% filter(species == "CUL4") %>% pull(Tp) %>% quantile(0.05))

branch_node_CUL9 <- beast_species_tree@phylo$tip.label == "CUL9"
index_CUL9 <- Ntip(beast_species_tree) - which(branch_node_CUL9)  - 6
dNdS_max_CUL9 <- as.numeric(dNdS_datation_df %>% filter(species == "CUL9") %>% pull(Tp) %>% quantile(0.95))
dNdS_min_CUL9 <- as.numeric(dNdS_datation_df %>% filter(species == "CUL9") %>% pull(Tp) %>% quantile(0.05))


p <- 
  ggtree(beast_species_tree) + 
  geom_range("CI_height", color='#FE6100', size=2, alpha=.5) + 
  geom_segment(aes(x = (-912361 * 3 / 1000000), xend = (-595421 * 3 / 1000000), y = index_prietella + 0.2, yend = index_prietella + 0.2), color = "#0072B2", size = 2, alpha=.5) +
  geom_segment(aes(x = (-dNdS_max_prietella / 1000000), xend = (-dNdS_min_prietella / 1000000), y = index_prietella - 0.2, yend = index_prietella - 0.2), color = "#CC79A7", size = 2, alpha=.5) +
  
  geom_segment(aes(x = (-72461 * 3 / 1000000), xend = (-841 * 3 / 1000000), y = index_CHM6 + 0.2, yend = index_CHM6 + 0.2), color = "#0072B2", size = 2, alpha=.5) +
  geom_segment(aes(x = (-dNdS_max_CHM6 / 1000000), xend = (-dNdS_min_CHM6 / 1000000), y = index_CHM6 - 0.2, yend = index_CHM6 - 0.2), color = "#CC79A7", size = 2, alpha=.5) +
  
  geom_segment(aes(x = (-72511 * 3 / 1000000), xend = (-841 * 3 / 1000000), y = index_Trichomycterus + 0.2, yend = index_Trichomycterus + 0.2), color = "#0072B2", size = 2, alpha=.5) +
  geom_segment(aes(x = (-dNdS_max_Trichomycterus / 1000000), xend = (-dNdS_min_Trichomycterus / 1000000), y = index_Trichomycterus - 0.2, yend = index_Trichomycterus - 0.2), color = "#CC79A7", size = 2, alpha=.5) +
  
  geom_segment(aes(x = (-136971 * 3 / 1000000), xend = (-24971 * 3 / 1000000), y = index_CSV83 + 0.2, yend = index_CSV83 + 0.2), color = "#0072B2", size = 2, alpha=.5) +
  geom_segment(aes(x = (-dNdS_max_CSV83 / 1000000), xend = (-dNdS_min_CSV83 / 1000000), y = index_CSV83 - 0.2, yend = index_CSV83 - 0.2), color = "#CC79A7", size = 2, alpha=.5) +
  
  geom_segment(aes(x = (-157681 * 3 / 1000000), xend = (-36621 * 3 / 1000000), y = index_CUL4 + 0.2, yend = index_CUL4 + 0.2), color = "#0072B2", size = 2, alpha=.5) +
  geom_segment(aes(x = (-dNdS_max_CUL4 / 1000000), xend = (-dNdS_min_CUL4 / 1000000), y = index_CUL4 - 0.2, yend = index_CUL4 - 0.2), color = "#CC79A7", size = 2, alpha=.5) +
  
  geom_segment(aes(x = (-212321 * 3 / 1000000), xend = (-66741 * 3 / 1000000), y = index_CUL9 + 0.2, yend = index_CUL9 + 0.2), color = "#0072B2", size = 2, alpha=.5) +
  geom_segment(aes(x = (-dNdS_max_CUL9 / 1000000), xend = (-dNdS_min_CUL9 / 1000000), y = index_CUL9 - 0.2, yend = index_CUL9 - 0.2), color = "#CC79A7", size = 2, alpha=.5) +
  
  coord_geo(xlim = c(-100, 0), ylim = c(-2, Ntip(beast_species_tree)+1),
            neg = TRUE, abbrv = TRUE, dat="epochs") +
  geom_tiplab() +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank())

revts(p)


ggtree(beast_species_tree) + 
  geom_tiplab()



#### Pearson correlation between the two datations ---------------------------------

GLS_function <- function(gls_rslt) {
  GLS_cc <- coef(gls_rslt)
  GLS_f <- function(x) GLS_cc[1] + GLS_cc[2]*x
  return(GLS_f)
}

#First do a linear regression to find the genetation time (w and without Prietella)

Summary_datation_simp <- Summary_datation %>% dplyr::select(species, Best_GenNb) 
dNdS_datation_df_simp <- Summary_dNdS_datation_noclean %>% filter(bootstrap_nb == "original_aln") %>% dplyr::select(species, Tp)
Combined_datation_df <- left_join(Summary_datation_simp, dNdS_datation_df_simp, by="species")



lm_pseudo_dNdS <- 
  lm(data = Combined_datation_df, 
     Tp ~ Best_GenNb)
summary_lm_pseudo_dNdS <- summary(lm_pseudo_dNdS)
pseudo_dNdS_function <- GLS_function(lm_pseudo_dNdS)
pseudo_dNdS_r2 = formatC(summary_lm_pseudo_dNdS$r.squared, digits = 2)
pseudo_dNdS_pval = formatC(summary_lm_pseudo_dNdS$coefficients[8], digits = 3)

Combined_datation_df %>%
  ggplot(., aes(x=Best_GenNb, y=Tp)) +
  geom_point() +
  theme_classic() +
  stat_function(fun = pseudo_dNdS_function, color="black") +
  labs(subtitle = bquote(R^2 ~ "=" ~ .(pseudo_dNdS_r2) ~ ";" ~ pvalue ~ "=" ~.(pseudo_dNdS_pval))) +
  xlab("Datation - Peudogenes number (Generations)") +
  ylab("Datation - dN/dS (Years)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 




lm_pseudo_dNdS <- 
  lm(data = Combined_datation_df %>% filter(species != "Prietella_phreatophila"), 
     Tp ~ Best_GenNb)
summary_lm_pseudo_dNdS <- summary(lm_pseudo_dNdS)
pseudo_dNdS_function <- GLS_function(lm_pseudo_dNdS)
pseudo_dNdS_r2 = formatC(summary_lm_pseudo_dNdS$r.squared, digits = 2)
pseudo_dNdS_pval = formatC(summary_lm_pseudo_dNdS$coefficients[8], digits = 3)

Combined_datation_df %>%
  filter(species != "Prietella_phreatophila") %>%
  ggplot(., aes(x=Best_GenNb, y=Tp)) +
  geom_point() +
  theme_classic() +
  stat_function(fun = pseudo_dNdS_function, color="black") +
  labs(subtitle = bquote(R^2 ~ "=" ~ .(pseudo_dNdS_r2) ~ ";" ~ pvalue ~ "=" ~.(pseudo_dNdS_pval))) +
  xlab("Datation - Peudogenes number (Generations)") +
  ylab("Datation - dN/dS (Years)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 


#Now estimate the generation time based on the reconciliation of the two datations


generation_time_df <- as.data.frame(NULL)
for(curr_species in cave_species){
  curr_TP_df <- 
    Summary_dNdS_datation_noclean %>%
    filter(species == curr_species) %>%
    sample_n(size = 10000, replace=TRUE) %>%
    dplyr::select(Tp)
  
  
  curr_Gen_df <- 
    Datation_results_df_normal %>%
    filter(species == curr_species) %>%
    filter(probability > 0.05) %>%
    sample_n(size = 10000, replace = TRUE, weight = probability) %>%
    dplyr::select(GenerationNb)
  
  curr_df <- as.data.frame(cbind(curr_TP_df, curr_Gen_df))
  curr_df <- 
    curr_df %>% 
    mutate(species = curr_species) %>%
    mutate(generation_time = Tp/GenerationNb)
  
  generation_time_df <- rbind(generation_time_df, curr_df)
 
}



generation_time_df <- read.table("generation_time_df.tsv",
                                 sep="\t", header=TRUE)


generation_time_df %>%
  ggplot(., aes(x = generation_time, color = species)) +
  geom_density(size = 1.2)  + 
  scale_color_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9"
  )) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.subtitle = element_text(size = 16),
    legend.position = "none"
  ) +
  xlab("Generation time (years)") +
  xlim(0, 30) +
  geom_vline(xintercept = 8.5, color="darkblue", linetype="dashed") +
  geom_vline(xintercept = 4, color="darkblue", linetype="dashed") +
  geom_vline(xintercept = 3, color="darkblue", linetype="dashed") +
  geom_vline(xintercept = 10, color="darkblue", linetype="dashed") +
  geom_vline(xintercept = 5, color="darkblue", linetype="dashed")


peak_densities <- 
  generation_time_df %>%
  group_by(species) %>%
  summarize(
    peak_generation_time = {
      dens <- density(generation_time)  # Compute density
      dens$x[which.max(dens$y)]        # Find x-value at maximum density
    }
  )



#Generation time = 3 years, pearson corr


Summary_datation_simp <- Summary_datation_simp %>% mutate(Tp_pseudo = Best_GenNb*3)

cor_test_result <- 
  cor.test(
    Combined_datation_df %>% pull(Tp_pseudo),
    Combined_datation_df %>% pull(Tp),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)


Combined_datation_df %>%
  ggplot(., aes(x=Tp_pseudo, y=Tp, color=species)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype="dashed", color="gray") +
  theme_classic() +
  geom_smooth(method = "lm", color = "black") +
  scale_color_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9"
  )) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Datation - number of pseudogenes") +
  ylab("Datation - dN/dS")


##### Vision dN/dS   ---------------------------------

Vision_dNdS_df <- 
  as.data.frame(
    fread("dNdS_per_vision_gene.csv",
          header=FALSE,
          sep=",")
  )

colnames(Vision_dNdS_df) <- c("HOG","species_gene", "LB", "MLE", "UB", "dN", "dS")

#Remove saturated genes or genes with not enough mutations
Vision_dNdS_df <- 
  Vision_dNdS_df %>%  
  filter(dN < 1) %>% 
  filter(dS < 1) %>% 
  filter(dS > 0.01) %>% 
  filter(MLE < 5) 

#Add the species name
Vision_dNdS_df_filt <- as.data.frame(NULL)
for(curr_species in species_list){
  
  curr_species_df <- 
    Vision_dNdS_df %>%
    filter(grepl(curr_species, species_gene)) %>%
    mutate(species = curr_species)
  
  Vision_dNdS_df_filt <- rbind(Vision_dNdS_df_filt, curr_species_df)
}

#Compute a mean dN/dS per species

vision_meandNdS_per_sp_df <- 
  Vision_dNdS_df_filt %>%
  group_by(species) %>%
  summarise(mean_dNdS = mean(MLE),
            median_dNdS = median(MLE))
vision_meandNdS_per_sp_df <- as.data.frame(vision_meandNdS_per_sp_df)


#Perform a pGLS between the habitat and the mean dNdS
species_vision_dNdS_df <- 
  left_join(species_df, vision_meandNdS_per_sp_df, by="species")


caper_dNdS_vision <- 
  comparative.data(phy = species_tree, 
                   data = species_vision_dNdS_df,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)

viison_dNdS_vs_Habitat_mean <-
  pgls(mean_dNdS ~ Habitat, 
       data = caper_dNdS_vision, 
       lambda = "ML")
summary(viison_dNdS_vs_Habitat_mean)

vision_dNdS_vs_Habitat_median <-
  pgls(median_dNdS ~ Habitat, 
       data = caper_dNdS_vision, 
       lambda = "ML")
summary(vision_dNdS_vs_Habitat_median)



species_vision_dNdS_df %>%
  ggplot(., aes(x=Habitat, y=mean_dNdS)) +
  geom_boxplot()


Vision_dNdS_df_filt %>%
  ggplot(., aes(x=species, y=MLE)) +
  geom_boxplot() 

env_df <- ace_values_labels_all %>% filter(! grepl("Node", label)) %>% dplyr::select(Cave, label)
colnames(env_df) <- c("Cave", "species")
Vision_dNdS_df_filt <- left_join(Vision_dNdS_df_filt, env_df, by="species")


Vision_dNdS_df_filt %>%
  ggplot(., aes(x=species, y=MLE)) +
  geom_violin() +
  geom_jitter(position=position_jitter(0.2)) +
  
  ylab("dN/dS") +
  theme_classic() +
  theme(axis.text=element_text(size=5),
        axis.title=element_text(size=5),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")



Vision_dNdS_df_filt_plot <- Vision_dNdS_df_filt %>% dplyr::select(species, MLE, Cave)
colnames(Vision_dNdS_df_filt_plot) <- c("species", "MLE", "Habitat")
p <- 
  ggtree(species_tree, size=2) %<+% 
  ace_values_labels_all +
  aes(color=as.numeric(Cave)) +
  scale_color_viridis() + 
  ggnewscale::new_scale_colour() + 
  geom_tiplab(size = 3, offset = 0.01) +
  theme(legend.position="none") +
  new_scale_fill() +
  new_scale_color()

p +
  geom_fruit(
    data=Vision_dNdS_df_filt_plot,
    geom=geom_boxplot,
    mapping = aes(
      x=MLE,
      y=species,
      fill=as.factor(Habitat),
      color=as.factor(Habitat)
    ),
    size=1,
    pwidth=3,
    offset = 0.8,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3,
      limits = c(0, 2)    
    ),
    grid.params=list()
  ) +
  scale_fill_manual(values = c("0" = "gray", "1" = "#FFC107")) +
  scale_color_manual(values = c("0" = "gray", "1" = "#FFC107")) +
  xlim_tree(0.5)


Vision_dNdS_df_filt %>%
  mutate(sp_color = if_else(
    Cave == 1,
    species,
    "Surface"
  )) %>%
  ggplot(., aes(x=reorder(species, MLE, mean), y=MLE, fill=as.factor(sp_color))) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9",
    "Surface" = "lightgray"
  )) +
  ylab("dN/dS") +
  xlab("Species") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  ylim(0, 2.5)



#Make pairwise comparisons

pairwise_wilcox_dNdS_vision <- 
  pairwise.wilcox.test(
    x=Vision_dNdS_df_filt %>% pull(MLE),
    g=Vision_dNdS_df_filt %>% pull(species),
    p.adjust.method="bonferroni"
  )
my_pvals <- pairwise_wilcox_dNdS_vision$p.value
my_pvals[my_pvals > 0.05 & !is.na(my_pvals)] <- NA
my_transformed_pvals=-log10(my_pvals)

melted_pvals <- melt(my_transformed_pvals, na.rm = FALSE)  # Keep NAs for coloring
colnames(melted_pvals) <- c("Row", "Column", "Value")



ggplot(melted_pvals, aes(x = Row, y = Column, fill = Value)) +
  geom_tile(color = "black") +  # Add borders to tiles
  scale_fill_gradient2(
    low = "lightcoral",
    mid = "firebrick2",
    high = "firebrick4",
    midpoint = 0,
    na.value = "white"  # White for NA values
  ) +
  labs(x = "", y = "", fill = "Log P-value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(),
    panel.grid = element_blank()
  ) +
  coord_fixed()  


#### Load Other candidate genes tables  ---------------------------------

#Load the sequences table
candidate_sequences_df <- 
  read.table("Other_genes/table_other_genes.tsv",
             sep="\t",
             header=FALSE)
colnames(candidate_sequences_df) <- 
  c("gene_clade","species","gene_name","genomic_position","CDS_length","gene_type","note",
    "exon_Count","sequence")
candidate_sequences_df <- candidate_sequences_df %>% filter(species != "") %>% filter(gene_type != "Not found")

#Load the LoF tables, with the read count frequencies
cand_stop_df <- read.table("Other_genes/Stop_table.tsv", sep="\t",header=FALSE)
colnames(cand_stop_df) <- c("gene_name", "species", "lof_type", "aa_pos")
cand_stop_df <- cand_stop_df %>% mutate(cds_pos = aa_pos*3) %>% mutate(nb_nuc = "") %>% dplyr::select(c("gene_name", "species", "lof_type", "cds_pos", "nb_nuc", "aa_pos"))

cand_frameshift_df <- read.table("Other_genes/Frameshift_table.tsv", sep="\t",header=FALSE)
colnames(cand_frameshift_df) <- c("gene_name", "species", "lof_type", "cds_pos", "nb_nuc")

cand_read_counts_df <- read.table("Other_genes/Read_counts_table.tsv", sep="\t",header=FALSE)
colnames(cand_read_counts_df) <- c("gene_name", "species", "lof_type", "aa_pos", "position", "support_reads", "non_support_reads")

#Merge tables
cand_stop_df <- left_join(cand_stop_df, cand_read_counts_df, by=c("gene_name", "species", "lof_type", "aa_pos"))
cand_stop_df <- cand_stop_df %>% dplyr::select(-aa_pos)
colnames(cand_read_counts_df) <- c("gene_name", "species", "lof_type", "cds_pos", "position", "support_reads", "non_support_reads")
cand_frameshift_df <- left_join(cand_frameshift_df, cand_read_counts_df, by=c("gene_name", "species", "lof_type", "cds_pos"))

cand_lof_df <- rbind(cand_stop_df, cand_frameshift_df)
cand_lof_df <- 
  cand_lof_df %>% 
  mutate(lof_name_temp = paste(gene_name,lof_type, sep="_")) %>%
  mutate(lof_name = paste(lof_name_temp,cds_pos, sep="_")) %>%
  dplyr::select(-lof_name_temp)



#Verify that the LoF table and Sequence table are similar in term of pseudogene nb


count_seq_table <- 
  as.data.frame(
    candidate_sequences_df %>%
      filter(gene_type == "Pseudogene") %>%
      group_by(species) %>%
      summarise(countP1 = n())
  )


count_lof_table <-
  as.data.frame(
    cand_lof_df %>%
      dplyr::select(species, gene_name) %>%
      distinct() %>%
      group_by(species) %>%
      summarise(countP2 = n())
  )


left_join(count_seq_table, count_lof_table, by="species") %>% filter(countP2 != countP1)



#Verify that the LoF table and Read count table are similar in term of pseudogene nb


count_lof_table_second <-
  as.data.frame(
    cand_read_counts_df %>%
      dplyr::select(species, gene_name) %>%
      distinct() %>%
      group_by(species) %>%
      summarise(countP3 = n())
  )


left_join(count_lof_table_second, count_lof_table, by="species") %>% filter(countP2 != countP3)


#Compute the genotype for each LoF mutation using a binomial law. 

cand_lof_df <- 
  cand_lof_df %>% 
  rowwise() %>%
  mutate(total_reads = non_support_reads + support_reads) %>%
  mutate(non_support_reads_prop = non_support_reads/total_reads) %>%
  mutate(support_reads_prop = support_reads/total_reads) %>% 
  mutate(minimum_read_number = min(support_reads, non_support_reads)) %>%
  mutate(binom_proba = pbinom(minimum_read_number, total_reads , p = 0.5)) %>%
  mutate(Genotype = case_when(
    binom_proba >= 0.05 & binom_proba <= 0.95 & non_support_reads != 0 & support_reads != 0 ~ "Heterozygous",
    binom_proba >= 0.05 & binom_proba <= 0.95 & non_support_reads == 0  & support_reads != 0 ~ "LoF_Homozygous",
    binom_proba >= 0.05 & binom_proba <= 0.95 & non_support_reads != 0  & support_reads == 0 ~ "no_LoF_Homozygous",
    binom_proba < 0.05 & (support_reads > non_support_reads) ~ "LoF_Homozygous",
    binom_proba > 0.95 & (support_reads > non_support_reads) ~ "LoF_Homozygous",
    binom_proba < 0.05 & (support_reads < non_support_reads) ~ "no_LoF_Homozygous",
    binom_proba > 0.95 & (support_reads < non_support_reads) ~ "no_LoF_Homozygous",
  )) %>%
  mutate(LoF_genotype = case_when(
    Genotype == "LoF_Homozygous" ~ 2,
    Genotype == "Heterozygous"  ~ 1,
    Genotype == "no_LoF_Homozygous" ~ 0
  )) %>%
  mutate(non_LoF_genotype = case_when(
    Genotype == "LoF_Homozygous" ~ 0,
    Genotype == "Heterozygous"  ~ 1,
    Genotype == "no_LoF_Homozygous" ~ 2
  ))  %>% 
  mutate(total_allele_nb = LoF_genotype + non_LoF_genotype) 

cand_lof_df <- as.data.frame(cand_lof_df)

cand_lof_df %>% filter(Genotype == "no_LoF_Homozygous")


#Count the number of LoF mutations per species

cand_lof_df %>%
  group_by(species, Genotype) %>%
  summarise(count = n()) %>%
  arrange(desc(count))


#Merge the LoF table with the Sequence table. Consider a gene as pseudogene
#if there is at-least one LoF present at heterozygous/homozygous state.
#Consider a pseudogene as heterozygous if there is no homozygous LoF 
#and homozygous if there is at-least one homozygous mutation

candidate_sequences_df_final <- as.data.frame(NULL)
for (curr_line in 1:nrow(candidate_sequences_df)){
  curr_species <- candidate_sequences_df[curr_line,]$species
  curr_gene_name <- candidate_sequences_df[curr_line,]$gene_name
  curr_gene_type <- candidate_sequences_df[curr_line,]$gene_type
  
  nb_homozygous_lof <- 
    nrow(cand_lof_df %>% 
           filter(species == curr_species) %>% 
           filter(gene_name == curr_gene_name) %>%
           filter(Genotype == "LoF_Homozygous"))
  
  nb_heterozygous_lof <- 
    nrow(cand_lof_df %>% 
           filter(species == curr_species) %>% 
           filter(gene_name == curr_gene_name) %>%
           filter(Genotype == "Heterozygous"))
  
  
  if(nb_homozygous_lof > 0){
    curr_line_mut <- candidate_sequences_df[curr_line,] %>% mutate(corrected_gene_type = "Pseudogene_Homozygous") 
  } 
  
  if(nb_homozygous_lof == 0 & nb_heterozygous_lof > 0){
    curr_line_mut <- candidate_sequences_df[curr_line,] %>% mutate(corrected_gene_type = "Pseudogene_Heterozygous") 
  }   
  
  if(nb_homozygous_lof == 0 & nb_heterozygous_lof == 0 & curr_gene_type == "Pseudogene"){
    curr_line_mut <- candidate_sequences_df[curr_line,] %>% mutate(corrected_gene_type = "Complete") 
  }   
  
  if(nb_homozygous_lof == 0 & nb_heterozygous_lof == 0 & curr_gene_type == "Incomplete"){
    curr_line_mut <- candidate_sequences_df[curr_line,] %>% mutate(corrected_gene_type = "Incomplete") 
  }  
  
  if(nb_homozygous_lof == 0 & nb_heterozygous_lof == 0 & curr_gene_type == "Complete"){
    curr_line_mut <- candidate_sequences_df[curr_line,] %>% mutate(corrected_gene_type = "Complete") 
  }  
  
  candidate_sequences_df_final <- rbind(candidate_sequences_df_final, curr_line_mut)
  
}


#Check which gene_type changed

nrow(candidate_sequences_df_final)
nrow(candidate_sequences_df)


candidate_sequences_df_final %>%
  filter(corrected_gene_type != gene_type) %>%
  dplyr::select(species, gene_name, gene_type, corrected_gene_type) 

#### Candidate genes analysis  ---------------------------------

#Summary of gene types 
summary_candidates_df <- 
  as.data.frame(
    candidate_sequences_df_final %>%
      group_by(species, corrected_gene_type) %>%
      summarise(count = n())
  )

summary_candidates_df_wide <- 
  as.data.frame(
    summary_candidates_df %>%
      pivot_wider(names_from = corrected_gene_type, values_from = count, values_fill = 0)
  ) %>%
  mutate(Total = Complete + Incomplete + Pseudogene_Homozygous) %>%
  mutate(Pseudogene = Pseudogene_Homozygous) %>%
  mutate(prop_pseudo = Pseudogene / Total)

summary_candidates_df_wide <- left_join(summary_candidates_df_wide, species_df, by="species")

#Compute the total number of LoF mutations

cand_nb_lof_df <- 
  as.data.frame(cand_lof_df %>%
                  filter(Genotype %in% c("LoF_Homozygous", "Heterozygous")) %>%
                  group_by(species) %>%
                  summarise(lof_nb = n()))

#Compute the total number of nucleotides

cand_nb_nuc_df <- 
  as.data.frame(candidate_sequences_df_final %>%
                  group_by(species) %>%
                  summarise(nt_nb = sum(CDS_length)))

#Combine tables

summary_candidates_df_wide <- left_join(summary_candidates_df_wide, cand_nb_lof_df, by="species")
summary_candidates_df_wide <- left_join(summary_candidates_df_wide, cand_nb_nuc_df, by="species")
summary_candidates_df_wide <- summary_candidates_df_wide %>% mutate(lof_nb = ifelse(is.na(lof_nb), 0, lof_nb))
summary_candidates_df_wide <- summary_candidates_df_wide %>% mutate(lof_per_nt = lof_nb/nt_nb)


#Print the table

gt(summary_candidates_df_wide %>%
     dplyr::select(species, Complete, Incomplete, Pseudogene, Total, prop_pseudo, nt_nb, lof_nb, lof_per_nt) %>%
     arrange(prop_pseudo))

#### Tile plot of candidate genes pseudogenes ---------------------------------



uniq_genes <- 
  candidate_sequences_df_final %>%
  filter(species != "Danio_rerio") %>%
  group_by(gene_name) %>%
  summarise(count = n()) %>%
  filter(count > 0) %>%
  pull(gene_name) %>%
  unique()

uniq_ind <- 
  candidate_sequences_df_final %>%
  pull(species) %>%
  unique()

all_combinations <- expand.grid(species = uniq_ind, gene_name = uniq_genes)


candidate_sequences_df_final_expanded <- 
  all_combinations %>%
  left_join(candidate_sequences_df_final, by = c("species", "gene_name")) %>%
  mutate(corrected_gene_type = ifelse(is.na(corrected_gene_type), "Not_found", corrected_gene_type))


order_genes_plot <- 
  c("plaat1","cry1a", "cry1b", "cry2a", "cry3", "per1b", "per2", "adamts20", "fhl2b", "gch2", "mc1r",
    "myo7ab", "oca2", "pax7b", "pmelb", "slc2a11b", "tmem33", "trpm1a", "tspan10", "tyrp1a", 
    "tyrp1b")


candidate_sequences_df_final_expanded$gene_name <- 
  factor(candidate_sequences_df_final_expanded$gene_name,
         levels = order_genes_plot)

candidate_sequences_df_final_expanded$corrected_gene_type <- 
  factor(candidate_sequences_df_final_expanded$corrected_gene_type,
         levels = c("Complete", "Incomplete", "Pseudogene_Homozygous", "Pseudogene_Heterozygous", "Not_found")) 


species_tree_len_woDr <- drop.tip(species_tree_len, "Danio_rerio") 

p1 <- 
  ggtree(species_tree_len_woDr, size=2) %<+% 
  ace_values_labels_all +
  aes(color=as.numeric(Cave)) +
  scale_color_viridis() + 
  ggnewscale::new_scale_colour() + 
  geom_tiplab(size = 3, offset = 0.01) + 
  theme(legend.position='none') +
  geom_facet(panel = "Gene", 
             data = candidate_sequences_df_final_expanded, 
             geom = geom_tile, 
             aes(x = as.integer(gene_name), 
                 fill = corrected_gene_type), 
             width = 1, 
             color="black") +
  scale_fill_manual(values =
                      c("Complete"="#009E73",
                        "Incomplete"="#56B4E9",
                        "Pseudogene_Homozygous"="#D55E00",
                        "Pseudogene_Heterozygous" = "#E69F00",
                        "Not_found" = "gray"))  +
  theme(legend.position = "none",
        strip.text = element_blank()) +
  xlim_tree(450) 


#Check which genes have been pseudogenized in internal branches of the tree (common LOF)


lof_mult_species_df <- 
  as.data.frame(
    vision_lof_df %>%
      filter(Genotype != "no_LoF_Homozygous") %>%
      group_by(gene_name, lof_type, cds_pos) %>%
      summarise(count = n()) %>%
      filter(count >= 2)
  )

lof_mult_species_df_filt <- as.data.frame(NULL)
for(curr_line in 1:nrow(lof_mult_species_df)){
  curr_gene <- lof_mult_species_df[curr_line,]$gene_name
  curr_type <- lof_mult_species_df[curr_line,]$lof_type
  curr_pos <- lof_mult_species_df[curr_line,]$cds_pos
  
  curr_species <- vision_lof_df %>% 
    filter(gene_name == curr_gene & lof_type == curr_type & cds_pos == curr_pos) %>% pull(species)
  
  curr_species_str <- paste(curr_species, collapse = ',')
  
  curr_df <- lof_mult_species_df[curr_line,] %>% mutate(species = curr_species_str)
  lof_mult_species_df_filt <- rbind(lof_mult_species_df_filt, curr_df)
  
}

lof_mult_species_df_filt %>%
  dplyr::select(gene_name, species) %>%
  distinct() %>% arrange(species)

#### Synthesis - ALL cavefishes  ---------------------------------

all_cavefish_df <- 
  read.table("All_cavefishes_table.vf.tsv",
             header=TRUE,
             sep="\t")

all_cavefish_df <- 
  all_cavefish_df %>%
  mutate(Summary_all_sp = if_else(if_any(-Gene, ~ . == "Pseudogene"), "Pseudogene", "Functional"))

all_cavefish_df_long <- 
  all_cavefish_df %>%
  pivot_longer(!Gene, names_to = "species", values_to = "gene_type")

all_cavefish_df_long$gene_type %>% unique()
all_cavefish_df_long$species %>% unique()
all_cavefish_df_long$Gene %>% unique()


order_species_plot <- 
  rev(c("Astyanax_mexicanus", "Lucifuga_dentata", "Lucifuga_gibarensis",
    "Sinocychlocheilus_anshuiensis", "Sinocychlocheilus_rhinocerous", 
    "Troglichthys_rosae", "Typhlichthys_subterraneus", "Typhlichthys_eigenmanni", "Speoplatyrhinus_poulsoni",
    "Amblyopsis_hoosieri", "Amblyopsis_spelaea",
    "Prietella_phreatophila", "Trichomycterus_rosablanca", "CHM6", "CUL4", "CSV83", "CUL9",
    "Others", "Summary_all_sp"))


order_genes_plot <- 
  c("rh1.1", "rh1.2", "rh2.1", "rh2.2", "lws", "sws1", "sws2", "exorh", "opn3", "opn4m1", 
    "opn4m2", "opn4m3", "opn4x1", "opn4x2", "opn5", "opn6a", "opn6b", "opn7a", "opn7b",
    "opn7c", "opn7d", "opn8a", "opn8b", "opn8c", "opn9", "parapinopsin-1", "parapinopsin-2",
    "parietopsin", "rgr1", "rgr2", "rrh", "tmt1a", "tmt1b", "tmt2a", "tmt2b", "tmt3a", "tmt3b",
    "va1", "va2", "cryaa", "cryba1b", "cryba1l1", "cryba2a", "cryba2b", "cryba4", 
    "crybb1", "crybb1l1", "crybb1l2", "crybb1l3", "crybgx", "crygm5", "crygn2", 
    "rpe65a", "rpe65b", "arr3a", "arr3b", "saga", "sagb", "gcap1", "gcap2", "gcap3", "gcap4",
    "gcap5", "gcap7", "grk1a", "grk1b", "grk7a", "grk7b", "rcv1a", "rcv1b", "rcv2", "gc2", "gc3", 
    "gucy2f", "pde6a", "pde6b", "pde6c", "pde6ga", "pde6gb", "pde6ha", "pde6hb",
    "gnat1", "gnat2", "gnb1a", "gnb1b", "gnb3a", "gnb3b", "gngt1", "gngt2", "gja8a", "gja8b")

all_cavefish_df_long$species <- 
  factor(all_cavefish_df_long$species,
         levels = order_species_plot)

all_cavefish_df_long$Gene <- 
  factor(all_cavefish_df_long$Gene,
         levels = order_genes_plot)


all_cavefish_df_long %>%
  ggplot(., aes(x = Gene, y = species, fill = gene_type)) +
  geom_tile(color = "black") +
  scale_fill_manual(values =
                      c("Functional"="#009E73",
                        "Pseudogene"="#D55E00",
                        "NO_DATA" = "white",
                        "Unknown" = "white",
                        "Missing" = "gray")) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic")) +
  labs(x = "Gene", y = "Species", fill = "Gene Type")


all_cavefish_df_long %>%
  filter(species == "Summary_all_sp") %>%
  group_by(gene_type) %>%
  summarise(n())

##### Genome-Wide dN/dS   ---------------------------------

GW_dNdS_df <- 
  as.data.frame(
    fread("dNdS_per_orthogroup.csv",
        header=FALSE,
        sep=",")
  )

colnames(GW_dNdS_df) <- c("HOG","species_gene", "LB", "MLE", "UB", "dN", "dS")

#Remove saturated genes or genes with not enough mutations
GW_dNdS_df <- 
  GW_dNdS_df %>%  
  filter(dN < 1) %>% 
  filter(dS < 1) %>% 
  filter(dS > 0.01) %>% 
  filter(MLE < 5) 

#Add the species name
GW_dNdS_df_filt <- as.data.frame(NULL)
for(curr_species in species_list){
  
  curr_species_df <- 
    GW_dNdS_df %>%
    filter(grepl(curr_species, species_gene)) %>%
    mutate(species = curr_species)
  
  GW_dNdS_df_filt <- rbind(GW_dNdS_df_filt, curr_species_df)
}

#Compute a mean dN/dS per species

meandNdS_per_sp_df <- 
  GW_dNdS_df_filt %>%
  group_by(species) %>%
  summarise(mean_dNdS = mean(MLE),
            median_dNdS = median(MLE))
meandNdS_per_sp_df <- as.data.frame(meandNdS_per_sp_df)


#Perform a pGLS between the habitat and the mean dNdS
species_dNdS_df <- 
  left_join(species_df, meandNdS_per_sp_df, by="species")


caper_dNdS <- 
  comparative.data(phy = species_tree, 
                   data = species_dNdS_df,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)

dNdS_vs_Habitat_mean <-
  pgls(mean_dNdS ~ Habitat, 
       data = caper_dNdS, 
       lambda = "ML")
summary(dNdS_vs_Habitat_mean)

dNdS_vs_Habitat_median <-
  pgls(median_dNdS ~ Habitat, 
       data = caper_dNdS, 
       lambda = "ML")
summary(dNdS_vs_Habitat_median)




species_dNdS_df %>%
  ggplot(., aes(x=Habitat, y=mean_dNdS)) +
  geom_boxplot()
  

GW_dNdS_df_filt %>%
  ggplot(., aes(x=species, y=MLE)) +
  geom_boxplot() 

env_df <- ace_values_labels_all %>% filter(! grepl("Node", label)) %>% dplyr::select(Cave, label)
colnames(env_df) <- c("Cave", "species")
GW_dNdS_df_filt <- left_join(GW_dNdS_df_filt, env_df, by="species")


GW_dNdS_df_filt %>%
  ggplot(., aes(x=species, y=MLE)) +
  geom_violin() +
  geom_jitter(position=position_jitter(0.2)) +

  ylab("dN/dS") +
  theme_classic() +
  theme(axis.text=element_text(size=5),
        axis.title=element_text(size=5),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")



GW_dNdS_df_filt_plot <- GW_dNdS_df_filt %>% dplyr::select(species, MLE, Cave)
colnames(GW_dNdS_df_filt_plot) <- c("species", "MLE", "Habitat")
p <- 
  ggtree(species_tree, size=2) %<+% 
  ace_values_labels_all +
  aes(color=as.numeric(Cave)) +
  scale_color_viridis() + 
  ggnewscale::new_scale_colour() + 
  geom_tiplab(size = 3, offset = 0.01) +
  theme(legend.position="none") +
  new_scale_fill() +
  new_scale_color()

p +
  geom_fruit(
    data=GW_dNdS_df_filt_plot,
    geom=geom_boxplot,
    mapping = aes(
      x=MLE,
      y=species,
      fill=as.factor(Habitat),
      color=as.factor(Habitat)
    ),
    size=1,
    pwidth=3,
    offset = 0.8,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3,
      limits = c(0, 2)    
      ),
    grid.params=list()
  ) +
  scale_fill_manual(values = c("0" = "gray", "1" = "#FFC107")) +
  scale_color_manual(values = c("0" = "gray", "1" = "#FFC107")) +
  xlim_tree(0.5)



#Now make a boxplot and a make an ANOVA


GW_dNdS_df_filt %>%
  ggplot(., aes(x=reorder(species, MLE, mean), y=MLE, fill=as.factor(Cave))) +
  geom_boxplot() +
  scale_fill_manual(values = c("0" = "gray", "1" = "#FFC107")) +
  ylab("dN/dS") +
  xlab("Species") +
  theme_classic() +
  theme(axis.text=element_text(size=3),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")


res.aov.dNdS_GW <- GW_dNdS_df_filt %>% anova_test(MLE ~ as.factor(Cave))
pwc.dNdS_GW <- GW_dNdS_df_filt %>% tukey_hsd(MLE ~ as.factor(Cave))


##### Compare GW and Vision dN/dS   ---------------------------------

colnames(species_vision_dNdS_df) <- c("species", "Habitat", "mean_dNdS_vision", "median_dNdS_vision")
colnames(species_dNdS_df) <- c("species", "Habitat", "mean_GW_vision", "median_GW_vision")

dNdS_df <- left_join(species_dNdS_df, species_vision_dNdS_df, by=c("species", "Habitat"))

dNdS_df <- dNdS_df %>% mutate(delta_dNdS_intra = mean_dNdS_vision - mean_GW_vision)

dNdS_df %>% arrange(delta_dNdS_intra)

##### Grantham distances results   ---------------------------------

grantham_df <- 
  as.data.frame(
    fread("Table_all_substitutions.grantham.txt",
          header=FALSE,
          sep=",")
  )

colnames(grantham_df) <- c("OGG", "gene_name", "substitution", "ancestral_aa", "new_aa",
                           "ancestral_aa_full", "new_aa_full", "grantham_score")

#Add species corresponding to branches
temp_df <- as.data.frame(NULL)
for(curr_species in species_list){
  
  curr_species_df <- 
    grantham_df %>%
    filter(grepl(curr_species, gene_name)) %>%
    mutate(species = curr_species)
  
  temp_df <- rbind(temp_df, curr_species_df)
}
grantham_df <- temp_df



meanGrantham_per_sp_df <- 
  grantham_df %>%
  group_by(species) %>%
  summarise(mean_grantham = mean(grantham_score),
            median_grantham = median(grantham_score))
meanGrantham_per_sp_df <- as.data.frame(meanGrantham_per_sp_df)


#Perform a pGLS between the habitat and the mean dNdS
species_grantham_df <- 
  left_join(species_df, meanGrantham_per_sp_df, by="species")


caper_grantham <- 
  comparative.data(phy = species_tree, 
                   data = species_grantham_df,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)

Grantham_vs_Habitat_mean <-
  pgls(mean_grantham ~ Habitat, 
       data = caper_grantham, 
       lambda = "ML")
summary(Grantham_vs_Habitat_mean)

Grantham_vs_Habitat_median <-
  pgls(median_grantham ~ Habitat, 
       data = caper_grantham, 
       lambda = "ML")
summary(Grantham_vs_Habitat_median)





env_df <- ace_values_labels_all %>% filter(! grepl("Node", label)) %>% dplyr::select(Cave, label)
colnames(env_df) <- c("Cave", "species")
grantham_df_filt <- left_join(grantham_df %>% dplyr::select(species, grantham_score), env_df, by="species")


colnames(grantham_df_filt) <- c("species", "grantham", "Habitat")
p <- 
  ggtree(species_tree, size=2) %<+% 
  ace_values_labels_all +
  aes(color=as.numeric(Cave)) +
  scale_color_viridis() + 
  ggnewscale::new_scale_colour() + 
  geom_tiplab(size = 3, offset = 0.01) +
  theme(legend.position="none") +
  new_scale_fill() +
  new_scale_color()

p +
  geom_fruit(
    data=grantham_df_filt,
    geom=geom_boxplot,
    mapping = aes(
      x=grantham,
      y=species,
      fill=as.factor(Habitat),
      color=as.factor(Habitat)
    ),
    size=1,
    pwidth=3,
    offset = 0.8,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3#,
      #limits = c(0, 85)    
    ),
    grid.params=list()
  ) +
  scale_fill_manual(values = c("0" = "gray", "1" = "#FFC107")) +
  scale_color_manual(values = c("0" = "gray", "1" = "#FFC107")) +
  xlim_tree(0.5)

#Now make a boxplot and a make an ANOVA


grantham_df_filt %>%
  ggplot(., aes(x=reorder(species, grantham, mean), y=grantham, fill=as.factor(Habitat))) +
  geom_boxplot() +
  scale_fill_manual(values = c("0" = "gray", "1" = "#FFC107")) +
  ylab("Grantham values") +
  xlab("Species") +
  theme_classic() +
  theme(axis.text=element_text(size=3),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")




##### Mutpred2 results   ---------------------------------

mutpred_df <- 
  as.data.frame(
    fread("All_Mutpred2_scores.csv",
          header=FALSE,
          sep=",")
  )

colnames(mutpred_df) <- c("OGG", "gene_name", "substitution", "mp_score")

#Add species corresponding to branches
temp_df <- as.data.frame(NULL)
for(curr_species in species_list){
  
  curr_species_df <- 
    mutpred_df %>%
    filter(grepl(curr_species, gene_name)) %>%
    mutate(species = curr_species)
  
  temp_df <- rbind(temp_df, curr_species_df)
}
mutpred_df <- temp_df



meanMutpred_per_sp_df <- 
  mutpred_df %>%
  group_by(species) %>%
  summarise(mean_mp = mean(mp_score),
            median_mp = median(mp_score))
meanMutpred_per_sp_df <- as.data.frame(meanMutpred_per_sp_df)


#Perform a pGLS between the habitat and the mean dNdS
species_mutpred_df <- 
  left_join(species_df, meanMutpred_per_sp_df, by="species")


caper_mutpred <- 
  comparative.data(phy = species_tree, 
                   data = species_mutpred_df,
                   names.col = species, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)

Mutpred_vs_Habitat_mean <-
  pgls(mean_mp ~ Habitat, 
       data = caper_mutpred, 
       lambda = "ML")
summary(Mutpred_vs_Habitat_mean)

Mutpred_vs_Habitat_median <-
  pgls(median_mp ~ Habitat, 
       data = caper_mutpred, 
       lambda = "ML")
summary(Mutpred_vs_Habitat_median)





env_df <- ace_values_labels_all %>% filter(! grepl("Node", label)) %>% dplyr::select(Cave, label)
colnames(env_df) <- c("Cave", "species")
mutpred_df_filt <- left_join(mutpred_df %>% dplyr::select(species, mp_score), env_df, by="species")


colnames(mutpred_df_filt) <- c("species", "mp_score", "Habitat")
p <- 
  ggtree(species_tree, size=2) %<+% 
  ace_values_labels_all +
  aes(color=as.numeric(Cave)) +
  scale_color_viridis() + 
  ggnewscale::new_scale_colour() + 
  geom_tiplab(size = 3, offset = 0.01) +
  theme(legend.position="none") +
  new_scale_fill() +
  new_scale_color()

p +
  geom_fruit(
    data=mutpred_df_filt,
    geom=geom_boxplot,
    mapping = aes(
      x=mp_score,
      y=species,
      fill=as.factor(Habitat),
      color=as.factor(Habitat)
    ),
    size=1,
    pwidth=3,
    offset = 1.4,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3#,
      #limits = c(0, 85)    
    ),
    grid.params=list()
  ) +
  scale_fill_manual(values = c("0" = "gray", "1" = "#FFC107")) +
  scale_color_manual(values = c("0" = "gray", "1" = "#FFC107")) +
  xlim_tree(0.5)

#Look at the densities

species_mp2 <- mutpred_df$species %>% unique()

#Now make a boxplot 


mutpred_df_filt %>%
  ggplot(., aes(x=reorder(species, mp_score, mean), y=mp_score, fill=as.factor(Habitat))) +
  geom_boxplot() +
  scale_fill_manual(values = c("0" = "gray", "1" = "#FFC107")) +
  ylab("Mutpred2 values") +
  xlab("Species") +
  theme_classic() +
  theme(axis.text=element_text(size=3),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")


## Go term enrichment (wilcoxon-rank test) of mutpred2 values

GO_annots_df <- read.table("OGG_GOterms.tsv",sep="\t",header=FALSE)
colnames(GO_annots_df) <- c("OGG", "GO")
GO_classes_df <- read.table("goterm_class.csv",sep="\t",header=FALSE)
colnames(GO_classes_df) <- c("GO", "class")
GO_annots_df <- left_join(GO_annots_df, GO_classes_df, by="GO")
BP_GO_annots_df <- GO_annots_df %>% filter(class == "biological_process") %>% dplyr::select(-class)
BP_GO_annots_list <- split(BP_GO_annots_df$GO, BP_GO_annots_df$OGG)

mutpred_df <- left_join(mutpred_df, species_df, by="species")

list_goterms <- BP_GO_annots_df %>% pull(GO) %>% unique()
GO_wilcox_mp_df <- as.data.frame(NULL)
#for(curr_GO in list_goterms){
#  
#  curr_OGGs <- BP_GO_annots_df %>% filter(GO == curr_GO) %>% pull(OGG)
#  
#  Surface_scores <-
#    mutpred_df %>% filter(OGG %in% curr_OGGs) %>% filter(Habitat == "Surface") %>% pull(mp_score)
#  
#  Cave_scores <-
#    mutpred_df %>% filter(OGG %in% curr_OGGs) %>% filter(Habitat == "Cave") %>% pull(mp_score)
#  
#  if(length(Surface_scores) >= 10 & length(Cave_scores) >= 10){
#    curr_wilc <- wilcox.test(Surface_scores,Cave_scores)
#    curr_pvalue <- curr_wilc$p.value
#    mean_Surface <- mean(Surface_scores)
#    mean_Cave <- mean(Cave_scores)
#    
#    curr_df <- as.data.frame(cbind(curr_GO,mean_Surface,mean_Cave,curr_pvalue))
#    colnames(curr_df) <- c("GO", "surface_mean", "cave_mean", "pvalue")
#    
#    GO_wilcox_mp_df <- rbind(GO_wilcox_mp_df, curr_df)
#  }
#  
#}

#write.table(GO_wilcox_mp_df, "GO_wilcox_mp_df",
#            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#if one don't want to re-run the loop, load the results directly
GO_wilcox_mp_df <- read.table("GO_wilcox_mp_df", header=TRUE, sep="\t")

#Transform columns as numeric
GO_wilcox_mp_df$pvalue <- as.numeric(GO_wilcox_mp_df$pvalue)
GO_wilcox_mp_df$cave_mean <- as.numeric(GO_wilcox_mp_df$cave_mean)
GO_wilcox_mp_df$surface_mean <- as.numeric(GO_wilcox_mp_df$surface_mean)

#Add GO-term descriptions
GO_term_desc <- read.table("goterm_desc.csv",sep="\t",header=FALSE)
colnames(GO_term_desc) <- c("GO", "term")
GO_wilcox_mp_df <- left_join(GO_wilcox_mp_df, GO_term_desc, by='GO')

BP_GO_annots_df_desc <- left_join(BP_GO_annots_df, GO_term_desc, by="GO")
#Apply a BH correction
raw_pvalues <- GO_wilcox_mp_df$pvalue
corr_pvalues <- p.adjust(raw_pvalues, method = "BH")
GO_wilcox_mp_df$corr_pvalue <- corr_pvalues
  
#Add the delta (cave mean - surface mean)

GO_wilcox_mp_df <- GO_wilcox_mp_df %>% mutate(delta_mean = cave_mean - surface_mean) 

#Plot the delta distribution

mean_delta <- GO_wilcox_mp_df %>% pull(delta_mean) %>% mean()
quantile1 <- quantile(GO_wilcox_mp_df %>% pull(delta_mean), probs = 0.01)
quantile99 <- quantile(GO_wilcox_mp_df %>% pull(delta_mean), probs = 0.99)


GO_wilcox_mp_df %>%
  ggplot(., aes(x=delta_mean)) +
  geom_histogram(bins= 60, color="black", fill="gray") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")  +
  geom_vline(xintercept = mean_delta, color="black") +
  geom_vline(xintercept = quantile1, color="#D55E00", linetype="dashed") +
  geom_vline(xintercept = quantile99, color="#D55E00", linetype="dashed") 
  
#Plot the p-value distribution

mean_p <- GO_wilcox_mp_df %>% mutate(log_p = -log10(corr_pvalue)) %>% pull(log_p) %>% mean()
                                       
quantile99_p <- 
  quantile(GO_wilcox_mp_df %>% mutate(log_p = -log10(corr_pvalue)) %>% pull(log_p), probs = 0.99)


GO_wilcox_mp_df %>%
  ggplot(., aes(x=-log10(corr_pvalue))) +
  geom_histogram(bins= 60, color="black", fill="gray") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")  +
  xlim(0, 20) +
  #geom_vline(xintercept = mean_p, color="black") +
  geom_vline(xintercept = quantile99_p, color="#D55E00", linetype="dashed") +
  xlab("log10(p.adj)") +
  ylab("# of GO term")


nrow(GO_wilcox_mp_df %>% filter(corr_pvalue < 0.05))
nrow(GO_wilcox_mp_df %>% filter(corr_pvalue >= 0.05))

nrow(GO_wilcox_mp_df %>% filter(corr_pvalue < 0.05) %>% filter(delta_mean > 0))
nrow(GO_wilcox_mp_df %>% filter(corr_pvalue < 0.05) %>% filter(delta_mean < 0))

#Extract the distributions tails

GO_wilcox_mp_df %>% 
  mutate(log_p = -log10(corr_pvalue)) %>%
  filter(log_p >= quantile99_p) %>%
  filter(delta_mean < 0) %>%
  arrange(desc(log_p))
  
GO_wilcox_mp_df %>% 
  mutate(log_p = -log10(corr_pvalue)) %>%
  filter(log_p >= quantile99_p) %>%
  filter(delta_mean > 0) %>%
  arrange(desc(log_p))


gt(
  GO_wilcox_mp_df %>% 
    mutate(log_p = -log10(corr_pvalue)) %>%
    filter(log_p >= quantile99_p) %>%
    filter(delta_mean > 0) %>%
    arrange(desc(log_p)) %>%
    dplyr::select(GO, term, corr_pvalue)
) 



curr_GO <- "GO:0007601"
curr_GO <- "GO:0002088"
curr_OGGs <- BP_GO_annots_df %>% filter(GO == curr_GO) %>% pull(OGG)
mutpred_df %>% 
  filter(OGG %in% curr_OGGs) %>%
  ggplot(., aes(x=Habitat, y=mp_score)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("Mutpred2 score") 

mutpred_df %>% 
  filter(OGG %in% curr_OGGs) %>%
  group_by(Habitat) %>%
  summarise(n())


gt(
  GO_wilcox_mp_df %>% 
    mutate(log_p = -log10(corr_pvalue)) %>%
    filter(log_p >= quantile99_p) %>%
    filter(delta_mean < 0) %>%
    arrange(desc(log_p))
)


curr_GO <- "GO:0071939"
curr_OGGs <- BP_GO_annots_df %>% filter(GO == curr_GO) %>% pull(OGG)
mutpred_df %>% 
  filter(OGG %in% curr_OGGs) %>%
  group_by(Habitat) %>%
  summarise(n())

#Make a nice figure for the GO-terms please

GO_wilcox_mp_df_fig <- GO_wilcox_mp_df
  
GO_wilcox_mp_df_fig[(GO_wilcox_mp_df_fig$term == "vitamin A import into cell"),"corr_pvalue"] <- 1e-100
GO_wilcox_mp_df_fig <- 
  GO_wilcox_mp_df_fig %>% mutate(delta_col = if_else(
  delta_mean < 0,
  "sf_higher",
  "cf_higher"
))


GO_wilcox_mp_df_fig %>% 
  mutate(log_p = -log10(corr_pvalue)) %>%
  filter(log_p >= quantile99_p) %>%
  ggplot(., aes(x=abs(delta_mean), y=-log10(corr_pvalue), color=delta_col, label = term)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("sf_higher" = "#6689E4", "cf_higher" = "#CCAC67")) +
  theme_classic() +
  geom_text_repel(size = 2, box.padding = 0.3, point.padding = 0.5, max.overlaps = Inf) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  xlab("Mean Cave-Surface MP2 Score Difference") +
  ylab("-log10(adj.pval)") +
  xlim(0, 0.3)
  

curr_GO <- "GO:0007601"
curr_OGGs <- BP_GO_annots_df %>% filter(GO == curr_GO) %>% pull(OGG)
mutpred_df %>% 
  filter(OGG %in% curr_OGGs) %>%
  ggplot(., aes(y=Habitat, x=mp_score)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  xlab("MP2 score") 


##### Mutpred2 enrichment test with only P. p    ---------------------------------

list_goterms <- BP_GO_annots_df %>% pull(GO) %>% unique()
#GO_wilcox_mp_df_onlyPP <- as.data.frame(NULL)
#for(curr_GO in list_goterms){
#  
#  curr_OGGs <- BP_GO_annots_df %>% filter(GO == curr_GO) %>% pull(OGG)
#  
#  Surface_scores <-
#    mutpred_df %>% filter(OGG %in% curr_OGGs) %>% filter(Habitat == "Surface") %>% pull(mp_score)
#  
#  Cave_scores <-
#    mutpred_df %>% filter(OGG %in% curr_OGGs) %>% filter(species == "Prietella_phreatophila") %>% pull(mp_score)
#  
#  if(length(Surface_scores) >= 10 & length(Cave_scores) >= 10){
#    curr_wilc <- wilcox.test(Surface_scores,Cave_scores)
#    curr_pvalue <- curr_wilc$p.value
#    mean_Surface <- mean(Surface_scores)
#    mean_Cave <- mean(Cave_scores)
#    
#    curr_df <- as.data.frame(cbind(curr_GO,mean_Surface,mean_Cave,curr_pvalue))
#    colnames(curr_df) <- c("GO", "surface_mean", "cave_mean", "pvalue")
#    
#    GO_wilcox_mp_df_onlyPP <- rbind(GO_wilcox_mp_df_onlyPP, curr_df)
#  }
#  
#}
#
#write.table(GO_wilcox_mp_df_onlyPP, "GO_wilcox_mp_df_onlyPP",
#            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#if one don't want to re-run the loop, load the results directly
GO_wilcox_mp_df_onlyPP <- read.table("GO_wilcox_mp_df_onlyPP", header=TRUE, sep="\t")

#Transform columns as numeric
GO_wilcox_mp_df_onlyPP$pvalue <- as.numeric(GO_wilcox_mp_df_onlyPP$pvalue)
GO_wilcox_mp_df_onlyPP$cave_mean <- as.numeric(GO_wilcox_mp_df_onlyPP$cave_mean)
GO_wilcox_mp_df_onlyPP$surface_mean <- as.numeric(GO_wilcox_mp_df_onlyPP$surface_mean)

#Add GO-term descriptions
GO_term_desc <- read.table("goterm_desc.csv",sep="\t",header=FALSE)
colnames(GO_term_desc) <- c("GO", "term")
GO_wilcox_mp_df_onlyPP <- left_join(GO_wilcox_mp_df_onlyPP, GO_term_desc, by='GO')

BP_GO_annots_df_desc <- left_join(BP_GO_annots_df, GO_term_desc, by="GO")
#Apply a BH correction
raw_pvalues <- GO_wilcox_mp_df_onlyPP$pvalue
corr_pvalues <- p.adjust(raw_pvalues, method = "BH")
GO_wilcox_mp_df_onlyPP$corr_pvalue <- corr_pvalues

#Add the delta (cave mean - surface mean)

GO_wilcox_mp_df_onlyPP <- GO_wilcox_mp_df_onlyPP %>% mutate(delta_mean = cave_mean - surface_mean) 


#Plot the delta distribution

mean_delta <- GO_wilcox_mp_df_onlyPP %>% pull(delta_mean) %>% mean()
quantile1 <- quantile(GO_wilcox_mp_df_onlyPP %>% pull(delta_mean), probs = 0.01)
quantile99 <- quantile(GO_wilcox_mp_df_onlyPP %>% pull(delta_mean), probs = 0.99)


GO_wilcox_mp_df_onlyPP %>%
  ggplot(., aes(x=delta_mean)) +
  geom_histogram(bins= 60, color="black", fill="gray") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")  +
  geom_vline(xintercept = mean_delta, color="black") +
  geom_vline(xintercept = quantile1, color="#D55E00", linetype="dashed") +
  geom_vline(xintercept = quantile99, color="#D55E00", linetype="dashed") 

#Plot the p-value distribution

mean_p <- GO_wilcox_mp_df_onlyPP %>% mutate(log_p = -log10(corr_pvalue)) %>% pull(log_p) %>% mean()

quantile99_p <- 
  quantile(GO_wilcox_mp_df_onlyPP %>% mutate(log_p = -log10(corr_pvalue)) %>% pull(log_p), probs = 0.99)

GO_wilcox_mp_df_onlyPP %>%
  ggplot(., aes(x=-log10(corr_pvalue))) +
  geom_histogram(bins= 60, color="black", fill="gray") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")  +
  xlim(0, 20) +
  #geom_vline(xintercept = mean_p, color="black") +
  geom_vline(xintercept = quantile99_p, color="#D55E00", linetype="dashed") +
  xlab("log10(p.adj)") +
  ylab("# of GO term")



nrow(GO_wilcox_mp_df_onlyPP %>% filter(corr_pvalue < 0.05))
nrow(GO_wilcox_mp_df_onlyPP %>% filter(corr_pvalue >= 0.05))

nrow(GO_wilcox_mp_df_onlyPP %>% filter(corr_pvalue < 0.05) %>% filter(delta_mean > 0))
nrow(GO_wilcox_mp_df_onlyPP %>% filter(corr_pvalue < 0.05) %>% filter(delta_mean < 0))

#Extract the distributions tails

GO_wilcox_mp_df_onlyPP %>% 
  mutate(log_p = -log10(corr_pvalue)) %>%
  filter(log_p >= quantile99_p) %>%
  filter(delta_mean < 0) %>%
  arrange(desc(log_p))

GO_wilcox_mp_df_onlyPP %>% 
  mutate(log_p = -log10(corr_pvalue)) %>%
  filter(log_p >= quantile99_p) %>%
  filter(delta_mean > 0) %>%
  arrange(desc(log_p))


gt(
  GO_wilcox_mp_df_onlyPP %>% 
    mutate(log_p = -log10(corr_pvalue)) %>%
    filter(log_p >= quantile99_p) %>%
    filter(delta_mean > 0) %>%
    arrange(desc(log_p)) %>%
    dplyr::select(GO, term, corr_pvalue)
) 



gt(
  GO_wilcox_mp_df_onlyPP %>% 
    mutate(log_p = -log10(corr_pvalue)) %>%
    filter(log_p >= quantile99_p) %>%
    filter(delta_mean < 0) %>%
    arrange(desc(log_p))
)



#Make a nice figure for the GO-terms 

GO_wilcox_mp_df_fig_onlyPP <- GO_wilcox_mp_df_onlyPP

GO_wilcox_mp_df_fig_onlyPP[(GO_wilcox_mp_df_fig_onlyPP$term == "vitamin A import into cell"),"corr_pvalue"] <- 1e-100
GO_wilcox_mp_df_fig_onlyPP <- 
  GO_wilcox_mp_df_fig_onlyPP %>% mutate(delta_col = if_else(
    delta_mean < 0,
    "sf_higher",
    "cf_higher"
  ))



GO_wilcox_mp_df_fig_onlyPP %>% 
  mutate(log_p = -log10(corr_pvalue)) %>%
  filter(log_p >= quantile99_p) %>%
  ggplot(., aes(x=abs(delta_mean), y=-log10(corr_pvalue), color=delta_col, label = term)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("sf_higher" = "#6689E4", "cf_higher" = "#CCAC67")) +
  theme_classic() +
  geom_text_repel(size = 2, box.padding = 0.3, point.padding = 0.5, max.overlaps = Inf) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  xlab("Mean Cave-Surface MP2 Score Difference") +
  ylab("-log10(adj.pval)") +
  xlim(0, 0.3)


##### Mutpred2 enrichment test without P. p    ---------------------------------

#list_goterms <- BP_GO_annots_df %>% pull(GO) %>% unique()
#GO_wilcox_mp_df_noPP <- as.data.frame(NULL)
#for(curr_GO in list_goterms){
#  
#  curr_OGGs <- BP_GO_annots_df %>% filter(GO == curr_GO) %>% pull(OGG)
#  
#  Surface_scores <-
#    mutpred_df %>% filter(OGG %in% curr_OGGs) %>% filter(Habitat == "Surface") %>% pull(mp_score)
#  
#  Cave_scores <-
#    mutpred_df %>% filter(OGG %in% curr_OGGs) %>% filter(Habitat == "Cave") %>% filter(species != "Prietella_phreatophila") %>% pull(mp_score)
#  
#  if(length(Surface_scores) >= 10 & length(Cave_scores) >= 10){
#    curr_wilc <- wilcox.test(Surface_scores,Cave_scores)
#    curr_pvalue <- curr_wilc$p.value
#    mean_Surface <- mean(Surface_scores)
#    mean_Cave <- mean(Cave_scores)
#    
#    curr_df <- as.data.frame(cbind(curr_GO,mean_Surface,mean_Cave,curr_pvalue))
#    colnames(curr_df) <- c("GO", "surface_mean", "cave_mean", "pvalue")
#    
#    GO_wilcox_mp_df_noPP <- rbind(GO_wilcox_mp_df_noPP, curr_df)
#  }
#  
#}

#write.table(GO_wilcox_mp_df_noPP, "GO_wilcox_mp_df_noPP",
#            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#if one don't want to re-run the loop, load the results directly
GO_wilcox_mp_df_noPP <- read.table("GO_wilcox_mp_df_noPP", header=TRUE, sep="\t")

#Transform columns as numeric
GO_wilcox_mp_df_noPP$pvalue <- as.numeric(GO_wilcox_mp_df_noPP$pvalue)
GO_wilcox_mp_df_noPP$cave_mean <- as.numeric(GO_wilcox_mp_df_noPP$cave_mean)
GO_wilcox_mp_df_noPP$surface_mean <- as.numeric(GO_wilcox_mp_df_noPP$surface_mean)

#Add GO-term descriptions
colnames(GO_term_desc) <- c("GO", "term")
GO_wilcox_mp_df_noPP <- left_join(GO_wilcox_mp_df_noPP, GO_term_desc, by='GO')

BP_GO_annots_df_desc <- left_join(BP_GO_annots_df, GO_term_desc, by="GO")
#Apply a BH correction
raw_pvalues <- GO_wilcox_mp_df_noPP$pvalue
corr_pvalues <- p.adjust(raw_pvalues, method = "BH")
GO_wilcox_mp_df_noPP$corr_pvalue <- corr_pvalues

#Add the delta (cave mean - surface mean)

GO_wilcox_mp_df_noPP <- GO_wilcox_mp_df_noPP %>% mutate(delta_mean = cave_mean - surface_mean) 




#Plot the delta distribution

mean_delta <- GO_wilcox_mp_df_noPP %>% pull(delta_mean) %>% mean()
quantile1 <- quantile(GO_wilcox_mp_df_noPP %>% pull(delta_mean), probs = 0.01)
quantile99 <- quantile(GO_wilcox_mp_df_noPP %>% pull(delta_mean), probs = 0.99)


GO_wilcox_mp_df_noPP %>%
  ggplot(., aes(x=delta_mean)) +
  geom_histogram(bins= 60, color="black", fill="gray") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none")  +
  geom_vline(xintercept = mean_delta, color="black") +
  geom_vline(xintercept = quantile1, color="#D55E00", linetype="dashed") +
  geom_vline(xintercept = quantile99, color="#D55E00", linetype="dashed") 

#Plot the p-value distribution

mean_p <- GO_wilcox_mp_df_noPP %>% mutate(log_p = -log10(corr_pvalue)) %>% pull(log_p) %>% mean()

quantile99_p <- 
  quantile(GO_wilcox_mp_df_noPP %>% mutate(log_p = -log10(corr_pvalue)) %>% pull(log_p), probs = 0.99)



GO_wilcox_mp_df_noPP %>%
  ggplot(., aes(x=-log10(corr_pvalue))) +
  geom_histogram(bins= 60, color="black", fill="gray") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none")  +
  xlim(0, 20) +
  #geom_vline(xintercept = mean_p, color="black") +
  geom_vline(xintercept = quantile99_p, color="#D55E00", linetype="dashed") +
  xlab("log10(p.adj)") +
  ylab("# of GO term")


nrow(GO_wilcox_mp_df_noPP %>% filter(corr_pvalue < 0.05))
nrow(GO_wilcox_mp_df_noPP %>% filter(corr_pvalue >= 0.05))

nrow(GO_wilcox_mp_df_noPP %>% filter(corr_pvalue < 0.05) %>% filter(delta_mean > 0))
nrow(GO_wilcox_mp_df_noPP %>% filter(corr_pvalue < 0.05) %>% filter(delta_mean < 0))

#Extract the distributions tails

GO_wilcox_mp_df_noPP %>% 
  mutate(log_p = -log10(corr_pvalue)) %>%
  filter(log_p >= quantile99_p) %>%
  filter(delta_mean < 0) %>%
  arrange(desc(log_p))

GO_wilcox_mp_df_noPP %>% 
  mutate(log_p = -log10(corr_pvalue)) %>%
  filter(log_p >= quantile99_p) %>%
  filter(delta_mean > 0) %>%
  arrange(desc(log_p))


gt(
  GO_wilcox_mp_df_noPP %>% 
    mutate(log_p = -log10(corr_pvalue)) %>%
    filter(log_p >= quantile99_p) %>%
    filter(delta_mean > 0) %>%
    arrange(desc(log_p)) %>%
    dplyr::select(GO, term, corr_pvalue)
) 


curr_GO <- "GO:0007601"
curr_GO <- "GO:0002088"
curr_OGGs <- BP_GO_annots_df %>% filter(GO == curr_GO) %>% pull(OGG)
mutpred_df %>% 
  filter(OGG %in% curr_OGGs) %>%
  ggplot(., aes(x=Habitat, y=mp_score)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("Mutpred2 score") 

mutpred_df %>% 
  filter(OGG %in% curr_OGGs) %>%
  group_by(Habitat) %>%
  summarise(n())


gt(
  GO_wilcox_mp_df_noPP %>% 
    mutate(log_p = -log10(corr_pvalue)) %>%
    filter(log_p >= quantile99_p) %>%
    filter(delta_mean < 0) %>%
    arrange(desc(log_p))
)




curr_GO <- "GO:0071939"
curr_OGGs <- BP_GO_annots_df %>% filter(GO == curr_GO) %>% pull(OGG)
mutpred_df %>% 
  filter(OGG %in% curr_OGGs) %>%
  group_by(Habitat) %>%
  summarise(n())


#Make a nice figure for the GO-terms 

GO_wilcox_mp_df_fig_noPP <- GO_wilcox_mp_df_noPP

GO_wilcox_mp_df_fig_noPP <- 
  GO_wilcox_mp_df_fig_noPP %>% mutate(delta_col = if_else(
    delta_mean < 0,
    "sf_higher",
    "cf_higher"
  ))


GO_wilcox_mp_df_fig_noPP %>% 
  mutate(log_p = -log10(corr_pvalue)) %>%
  filter(log_p >= quantile99_p) %>%
  ggplot(., aes(x=abs(delta_mean), y=-log10(corr_pvalue), color=delta_col, label = term)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("sf_higher" = "#6689E4", "cf_higher" = "#CCAC67")) +
  theme_classic() +
  geom_text_repel(size = 2, box.padding = 0.3, point.padding = 0.5, max.overlaps = Inf) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  xlab("Mean Cave-Surface MP2 Score Difference") +
  ylab("-log10(adj.pval)") +
  xlim(0, 0.3)


##### Figure - Boxplot GW dN/dS and GW Mutpred2   ---------------------------------


GW_dNdS_df_filt %>%
  mutate(sp_color = if_else(
    Cave == 1,
    species,
    "Surface"
  )) %>%
  ggplot(., aes(x=reorder(species, MLE, mean), y=MLE, fill=as.factor(sp_color))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9",
    "Surface" = "lightgray"
  )) +
  ylab("dN/dS") +
  xlab("Species") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  ylim(0, 1)


sp_order_box <- 
  GW_dNdS_df_filt %>%
  group_by(species) %>%
  summarise(mean_MLE = mean(MLE)) %>%
  arrange(mean_MLE) %>%
  pull(species)

mutpred_df_filt$species <- 
  factor(mutpred_df_filt$species, levels = sp_order_box)


mutpred_df_filt %>%
  mutate(sp_color = if_else(
    Habitat == 1,
    species,
    "Surface"
  )) %>%
  ggplot(., aes(x=species, y=mp_score, fill=as.factor(sp_color))) +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9",
    "Surface" = "lightgray"
  )) +
  ylab("Mutpred2 values") +
  xlab("Species") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  ylim(0, 1)



#Make the link between dN/dS and mp2 scores

species_mp_dNdS <- left_join(species_mutpred_df, meandNdS_per_sp_df, by=c("species"))


cor_test_result <- 
  cor.test(
    species_mp_dNdS %>% pull(mean_dNdS),
    species_mp_dNdS %>% pull(mean_mp),
    method="pearson"
  )
estimate_p <- round(cor_test_result$estimate, 3)
pvalue_p <- cor_test_result$p.value
pvalue_p <- format.pval(pvalue_p, digits = 2)

species_mp_dNdS %>%
  mutate(sp_color = if_else(
    Habitat == "Cave",
    species,
    "Surface"
  )) %>%
  ggplot(., aes(x=mean_dNdS, y=mean_mp, color=as.factor(sp_color))) +
  geom_smooth(method = "lm", color = "black", alpha=0.1) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9",
    "Surface" = "lightgray"
  )) +
  xlab("Mean dN/dS") +
  ylab("Mean MutPred2 score") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") 


##### Compute cave only (wCO) dN/dS - Genome wide    ---------------------------------

comp_table <- 
  as.data.frame(
  rbind(
  cbind("CUL9", "RHP1", (8.188 * 1000000)),
  cbind("CSV83", "RUI2", (14.6389 * 1000000)),
  cbind("CUL4", "RUI2", (14.6389 * 1000000)),
  cbind("CHM6", "RSS1", (10.5187 * 1000000)),
  cbind("Trichomycterus_rosablanca", "RSS1", (15.1434 * 1000000)),
  cbind("Prietella_phreatophila", "Ameiurus_melas", (12.9362 * 1000000))
  )
)
  
comp_table$V3 <- as.numeric(comp_table$V3)

cave_only_GW_metrics_df <- as.data.frame(NULL)
for(curr_species in cave_species){
  
  close_surface_sp <- comp_table %>% filter(V1 == curr_species) %>% pull(V2)
    
  Tp <- Summary_dNdS_datation_noclean %>%
    filter(species == curr_species) %>%
    filter(bootstrap_nb == "original_aln") %>%
    pull(Tp)
  
  dNdScf <- species_dNdS_df %>% filter(species == curr_species) %>% pull(mean_dNdS)
  dNdSsf <- species_dNdS_df %>% filter(species == close_surface_sp) %>% pull(mean_dNdS)
  
  Tcs <- comp_table %>% filter(V1 == curr_species) %>% pull(V3)
  
  dNdS_co_100 <- (Tcs * (dNdScf - dNdSsf) + (Tp * dNdSsf))/Tp
  dNdS_co_80 <- (Tcs * (dNdScf - dNdSsf) + (0.8*Tp * dNdSsf))/(0.8 * Tp)
  

  curr_df <- as.data.frame(cbind(curr_species, Tp,Tcs, dNdScf, dNdSsf, dNdS_co_100, dNdS_co_80))
  colnames(curr_df) <- c("species", "Tp","Tcs", "wCF", "wSF","wCO_100", "wCO_80")
  
  cave_only_GW_metrics_df <- rbind(cave_only_GW_metrics_df, curr_df)
}


##### Cave age vs Heterozygosity    ---------------------------------

cave_ages_df <- 
  Summary_datation %>% filter(datatation_type == "normal") %>% dplyr::select(species, Best_GenNb)

H_df_age <- left_join(cave_ages_df, genomescope_df, by="species")

H_df_age %>%
  ggplot(., aes(x=Best_GenNb, y=hetero_max, color=species)) +
  geom_point() +
  scale_color_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9"
  )) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  xlab("Generation number") +
  ylab("Heterozygosity")




##### Cave age vs Genome-wide dN/dS   ---------------------------------

cave_ages_df <- 
  Summary_datation %>% filter(datatation_type == "normal") %>% dplyr::select(species, Best_GenNb)

dNdS_df_age <- left_join(cave_ages_df, meandNdS_per_sp_df, by="species")


dNdS_df_age %>%
  ggplot(., aes(x=Best_GenNb, y=mean_dNdS, color=species)) +
  geom_point() +
  scale_color_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9"
  )) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  xlab("Generation number") +
  ylab("G-W dN/dS")


pearson_dNdS_age_GW <- 
  cor.test(delta_dNdS_df_age$Best_GenNb, delta_dNdS_df_age$delta_dNdS, method = 'pearson')


##### Cave age vs Vision dN/dS   ---------------------------------


cave_ages_df <- 
  Summary_datation %>% filter(datatation_type == "normal") %>% dplyr::select(species, Best_GenNb)

dNdS_df_age <- left_join(cave_ages_df, species_vision_dNdS_df, by="species")




dNdS_df_age %>%
  ggplot(., aes(x=Best_GenNb, y=mean_dNdS, color=species)) +
  geom_point() +
  scale_color_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9"
  )) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  xlab("Generation number") +
  ylab("vision dN/dS")


pearson_dNdS_age_vision <- 
  cor.test(delta_dNdS_df_age$Best_GenNb, delta_dNdS_df_age$delta_dNdS, method = 'pearson')


##### Load RERconverge results   ---------------------------------

#Import trees and phenotypes
library(RERconverge)
cavefish_species <- c("Prietella_phreatophila","CHM6","CSV83","CUL4","CUL9","Trichomycterus_rosablanca")
RERw.all <- readRDS("RERw.all.RDS")
GW_trees.all <- readRDS("AllTrees.Phangorn.RDS")

all.anc <- foreground2Tree(treesObj = GW_trees.all,foreground = cavefish_species,clade = "ancestral")
all.all <- foreground2Tree(treesObj = GW_trees.all,foreground = cavefish_species,clade = "all")
all.term <- foreground2Tree(treesObj = GW_trees.all,foreground = cavefish_species,clade = "terminal")

all.anc.phenv <- tree2Paths(all.anc, GW_trees.all)
all.all.phenv <- tree2Paths(all.all, GW_trees.all)
all.term.phenv <- tree2Paths(all.term, GW_trees.all)
dev.off()

#Load RER associations results
RER_all_df <- read.table("all.all.res.perm", sep=",", header=TRUE)
RER_anc_df <- read.table("all.anc.res.perm", sep=",", header=TRUE)
RER_term_df <- read.table("all.term.res.perm", sep=",", header=TRUE)

RER_all_df <- RER_all_df %>% filter(! is.na(Rho))
RER_anc_df <- RER_anc_df %>% filter(! is.na(Rho))
RER_term_df <- RER_term_df %>% filter(! is.na(Rho))

#Combine tables with gene names

colnames(OGG_names_df) <- c("HOG", "gene_ID", "gene_name", "evalue_name")
RER_all_df <- left_join(RER_all_df, OGG_names_df, by="HOG")
RER_anc_df <- left_join(RER_anc_df, OGG_names_df, by="HOG")
RER_term_df <- left_join(RER_term_df, OGG_names_df, by="HOG")


#Load enrichment results
RER_all_enrich <- readRDS("enrichment.all.perm.RDS")
RER_anc_enrich <- readRDS("enrichment.anc.perm.RDS")
RER_term_enrich <- readRDS("enrichment.term.perm.RDS")


#Convert enrichment list to dataframes
nb_enrich_group <- length(RER_all_enrich)
RER_all_enrich_df <- as.data.frame(NULL)
for(count in 1:nb_enrich_group){
  curr_database <- names(RER_all_enrich[count])
  curr_df <- as.data.frame(RER_all_enrich[[count]]) %>% mutate(database = curr_database)
  RER_all_enrich_df <- rbind(RER_all_enrich_df, curr_df) 
}




nb_enrich_group <- length(RER_anc_enrich)
RER_anc_enrich_df <- as.data.frame(NULL)
for(count in 1:nb_enrich_group){
  curr_database <- names(RER_anc_enrich[count])
  curr_df <- as.data.frame(RER_anc_enrich[[count]])  %>% mutate(database = curr_database)
  RER_anc_enrich_df <- rbind(RER_anc_enrich_df, curr_df) 
}


#Add goterm names to enrichment tables

rownames_curr_df <- row.names(RER_all_enrich_df)
rownames_curr_df <- gsub("_.*", "", rownames_curr_df)
RER_all_enrich_df <- RER_all_enrich_df %>% mutate(GO_term_name = rownames_curr_df)

rownames_curr_df <- row.names(RER_anc_enrich_df)
rownames_curr_df <- gsub("_.*", "", rownames_curr_df)
RER_anc_enrich_df <- RER_anc_enrich_df %>% mutate(GO_term_name = rownames_curr_df)



##### Analyse RERconverge associations results   ---------------------------------
 
#Add the permulation results 
RER_all_df <- 
  RER_all_df %>%
  mutate(perm_category = case_when(
      permpval < 0.05 & P < 0.05 & Rho > 0 ~ "Signif_accelerated",
      permpval < 0.05 & P < 0.05 & Rho < 0 ~ "Signif_deccelerated",
      permpval > 0.05 ~ "Non_signif",
      P > 0.05 ~ "Non_signif"
    ))

RER_all_df$perm_category[is.na(RER_all_df$perm_category)] <- "Non_signif"

#Mark genes involved in vision 

GO_annots_df <- read.table("OGG_GOterms.tsv",sep="\t",header=FALSE)
colnames(GO_annots_df) <- c("OGG", "GO")
GO_classes_df <- read.table("goterm_class.csv",sep="\t",header=FALSE)
colnames(GO_classes_df) <- c("GO", "class")
GO_annots_df <- left_join(GO_annots_df, GO_classes_df, by="GO")
BP_GO_annots_df <- GO_annots_df %>% filter(class == "biological_process") %>% dplyr::select(-class)
BP_GO_annots_list <- split(BP_GO_annots_df$GO, BP_GO_annots_df$OGG)
GO_term_desc <- read.table("goterm_desc.csv",sep="\t",header=FALSE)
colnames(GO_term_desc) <- c("GO", "term")
BP_GO_annots_df_desc <- left_join(BP_GO_annots_df, GO_term_desc, by="GO")


BP_GO_annots_df_desc %>% 
  filter(
    grepl(" visual", term) |
      grepl("visual", term) | 
      grepl("eye", term) | 
      grepl("photorecep", term)    
  ) %>% pull(GO) %>% unique()



visual_OGGs <- 
  BP_GO_annots_df_desc %>% 
  filter(
    grepl(" visual", term) |
    grepl("visual", term) | 
    grepl("eye", term) | 
    grepl("photorecep", term)    
    ) %>%
  pull(OGG) %>%
  unique()

RER_all_df <- 
  RER_all_df %>%
  mutate(vision_related = if_else(
    HOG %in% visual_OGGs,
    "Yes", 
    "No"
  ))



#First plot the results


RER_all_df %>%
  ggplot(., aes(x=Rho, y=-log10(P), fill=perm_category, shape=vision_related, color=vision_related)) + 
  geom_point() +
  scale_fill_manual(values = c("Non_signif" = "black", "Signif_accelerated" = "#E69F00", 
                               "Signif_deccelerated" = "#56B4E9")) +
  scale_shape_manual(values = c("Yes" = 24, "No" = 21)) +
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) +
  geom_hline(yintercept = -log10(0.05), color="black", linetype="dashed") +
  geom_vline(xintercept = 0, color="black", linetype="dashed") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")



RER_all_df %>%
  ggplot(., aes(x=Rho, y=-log10(P))) + 
  geom_point(data = subset(RER_all_df, vision_related == "No"), 
             aes(fill = perm_category, shape = vision_related), 
             size = 1) + 
  geom_point(data = subset(RER_all_df, vision_related == "Yes"), 
             aes(fill = perm_category, shape = vision_related), 
             size = 3) + 
  scale_fill_manual(values = c("Non_signif" = "black", "Signif_accelerated" = "#E69F00", 
                                "Signif_deccelerated" = "#56B4E9")) +
  scale_shape_manual(values = c("Yes" = 24, "No" = 21)) +
  #scale_color_manual(values = c("Yes" = "red", "No" = "black")) +
  geom_hline(yintercept = -log10(0.05), color="black", linetype="dashed") +
  geom_vline(xintercept = 0, color="black", linetype="dashed") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

RER_all_df %>%
  filter(vision_related == "Yes") %>%
  ggplot(., aes(x=Rho, y=-log10(P), color=perm_category, shape=vision_related)) + 
  geom_point() +
  scale_color_manual(values = c("Non_signif" = "black", "Signif_accelerated" = "#E69F00", 
                                "Signif_deccelerated" = "#56B4E9")) +
  geom_hline(yintercept = -log10(0.05), color="black", linetype="dashed") +
  geom_vline(xintercept = 0, color="black", linetype="dashed") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") 


#Check visual related genes significant

RER_all_df %>%
  filter(vision_related == "Yes") %>%
  filter(Rho > 0.5)


RER_all_df %>%
  filter(vision_related == "Yes") %>%
  filter(P < 0.05 & permpval < 0.05) %>%
  arrange(Rho)

#Count significant genes

nrow(RER_all_df %>%
  filter(! is.na(P)))

nrow(RER_all_df %>%
  filter(P < 0.05 & permpval < 0.05) %>%
  filter(Rho <= 0))

nrow(RER_all_df %>%
       filter(P < 0.05 & permpval < 0.05) %>%
       filter(Rho > 0))

RER_accelerated_genes <- 
  RER_all_df %>%
  filter(P < 0.05 & permpval < 0.05) %>%
  filter(Rho > 0) %>% pull(HOG)

#Check which genes are the most significant

RER_all_df %>%
  filter(P < 0.05 & permpval < 0.05) %>%
  filter(Rho > 0) %>%
  arrange(desc(Rho)) %>%
  #arrange((P)) %>%
  dplyr::select(HOG, Rho, P, permpval, gene_ID, gene_name) %>%
  head(5)


plotRers(RERw.all, "N1.HOG0023808", phenv=all.all.phenv)


RER_all_df %>%
  filter(P < 0.05 & permpval < 0.05) %>%
  filter(Rho < 0) %>%
  arrange(Rho) %>%
  #arrange(P) %>%
  dplyr::select(HOG, Rho, P, permpval, gene_ID, gene_name) %>%
  head(5)

plotRers(RERw.all, "N1.HOG0022635", phenv=all.all.phenv)

##### Analyse RERconverge enrichments results   ---------------------------------


RER_all_enrich_df %>% pull(database) %>% unique()


#Biological processes 
RER_all_enrich_df %>% 
  filter(pval < 0.05 & permpval < 0.05) %>% 
  filter(stat < 0) %>%
  arrange(desc(stat)) %>%
  filter(database == "biological_process") %>%
  dplyr::select(stat, pval, permpval, GO_term_name) %>%
  pull(GO_term_name)


gt(
  RER_all_enrich_df %>% 
  filter(pval < 0.05 & permpval < 0.05) %>% 
  filter(stat > 0) %>%
  arrange(desc(stat)) %>%
  filter(database == "biological_process") %>%
  dplyr::select(stat, pval, permpval, GO_term_name)
)
  

gt(
  RER_all_enrich_df %>% 
    filter(pval < 0.05 & permpval < 0.05) %>% 
    filter(stat < 0) %>%
    arrange(stat) %>%
    filter(database == "biological_process") %>%
    dplyr::select(stat, pval, permpval, GO_term_name)
)




#PANTHER

gt(
  RER_all_enrich_df %>% 
    filter(pval < 0.05 & permpval < 0.05) %>% 
    filter(stat > 0) %>%
    arrange(desc(stat)) %>%
    filter(database == "PANTHER") %>%
    dplyr::select(stat, pval, permpval, GO_term_name)
)

gt(
  RER_all_enrich_df %>% 
    filter(pval < 0.05 & permpval < 0.05) %>% 
    filter(stat < 0) %>%
    arrange(stat) %>%
    filter(database == "PANTHER") %>%
    dplyr::select(stat, pval, permpval, GO_term_name)
)



##### Import and analyse RELAX results   ---------------------------------


#Import RELAX results and compute p-values
RELAX_df <- 
  read.table("RELAX_results.csv",
             header=FALSE,
             sep=",")
colnames(RELAX_df) <- c("OGG", "LRT", "K")
RELAX_df <- 
  as.data.frame(RELAX_df %>%
                  rowwise() %>%
                  mutate(pvalue = pchisq(LRT, df=1, lower.tail = FALSE)))

RELAX_df <- RELAX_df %>% filter(! is.na(pvalue))

list_pvalues <- RELAX_df$pvalue
corrected_pvalues <- p.adjust(list_pvalues, method = "BH")
RELAX_df <- 
  RELAX_df %>%
  mutate(BH_pvalue = corrected_pvalues)



OGG_names_df <- 
  read.table("OGG_Names_Product.tsv",
             header=FALSE,
             sep="\t")
colnames(OGG_names_df) <- c("OGG", "gene_ID", "gene_name", "evalue_name")

RELAX_df <- left_join(RELAX_df, OGG_names_df, by="OGG") %>% dplyr::select(-evalue_name)

#Add a column with RELAX results

RELAX_df <-
  RELAX_df %>%
  mutate(RELAX_rslt = case_when(
    pvalue <= 0.05 & K >= 1 ~ "Itensification",
    pvalue <= 0.05 & K < 1 ~ "Relaxation",
    pvalue > 0.05 ~ "Non_significant"
  ))
  
#Select genes found under accelerated evolution by 

RELAX_RERacc_df <- 
  RELAX_df %>%
  filter(OGG %in% RER_accelerated_genes)
  

#Look at if these genes are under relaxed or positive selection


RELAX_RERacc_df %>% 
  group_by(RELAX_rslt) %>%
  summarise(count = n())

RELAX_RERacc_df <- 
  RELAX_RERacc_df %>%
  mutate(K_normal = if_else(
    K > 5,
    5,
    K
  ))


# Add genes belonign to visual go-terms

RELAX_RERacc_df <- 
  RELAX_RERacc_df %>%
  mutate(vision_related = if_else(
    OGG %in% visual_OGGs,
    "Yes", 
    "No"
  ))



#Make a plot of the results  


RELAX_RERacc_df %>%
  ggplot(., aes(x=K_normal, y=-log10(pvalue), fill=RELAX_rslt, shape=vision_related)) + 
  geom_jitter(width = 0.1, height = 0.1) + 

  scale_fill_manual(values = c("Non_significant" = "black", "Itensification" = "#E69F00", 
                               "Relaxation" = "#56B4E9")) +
  scale_shape_manual(values = c("Yes" = 24, "No" = 21)) +
  geom_hline(yintercept = -log10(0.05), color="black", linetype="dashed") +
  geom_vline(xintercept = 1, color="black", linetype="dashed") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("K") 


RELAX_RERacc_df %>% filter(vision_related == "Yes")

#Look at genes evolving under relaxed selection

RELAX_RERacc_df %>%
  filter(RELAX_rslt == "Relaxation")

#Look at genes evolving under intensified selection

RELAX_RERacc_df %>%
  filter(RELAX_rslt == "Itensification")

##### Import cSUBST results -- Gene Tree  ---------------------------------

#Import branch informations
csubst_b.genetree <- 
  as.data.frame(
    fread("all_OGG_csubst_b.GeneTree.tsv",
          header=TRUE,
          sep="\t")
  )





#Import K = 2 results

csubst_K2.genetree <- 
  as.data.frame(
    fread("all_OGG_csubst_cb_2.GeneTree.tsv",
          header=TRUE,
          sep="\t")
  )





#Import K = 3 results


csubst_K3.genetree <- 
  as.data.frame(
    fread("all_OGG_csubst_cb_3.GeneTree.tsv",
          header=TRUE,
          sep="\t")
  )




#Add branch names to K2 and K3 tables

#Gene Tree
branch_id_name_df <- csubst_b.genetree  %>% dplyr::select(OGG, branch_name, branch_id)
colnames(branch_id_name_df) <- c("OGG", "branch_name_1", "branch_id_1")
csubst_K2.genetree <- left_join(csubst_K2.genetree, branch_id_name_df, by=c("OGG", "branch_id_1"))
csubst_K3.genetree <- left_join(csubst_K3.genetree, branch_id_name_df, by=c("OGG", "branch_id_1"))
colnames(branch_id_name_df) <- c("OGG", "branch_name_2", "branch_id_2")
csubst_K2.genetree <- left_join(csubst_K2.genetree, branch_id_name_df, by=c("OGG", "branch_id_2"))
csubst_K3.genetree <- left_join(csubst_K3.genetree, branch_id_name_df, by=c("OGG", "branch_id_2"))
colnames(branch_id_name_df) <- c("OGG", "branch_name_3", "branch_id_3")
csubst_K3.genetree <- left_join(csubst_K3.genetree, branch_id_name_df, by=c("OGG", "branch_id_3"))




#For CSUBST on gene trees, remove node values

csubst_K2.genetree <- 
  csubst_K2.genetree %>%
  filter(! grepl("Node", branch_name_1)) %>%
  filter(! grepl("Node", branch_name_2))

csubst_K3.genetree <- 
  csubst_K3.genetree %>%
  filter(! grepl("Node", branch_name_1)) %>%
  filter(! grepl("Node", branch_name_2)) %>%
  filter(! grepl("Node", branch_name_3))




#Add the species names

csubst_K2.genetree <- 
  csubst_K2.genetree %>%
  mutate(species1 = case_when(
    grepl("Ameiurus_melas", branch_name_1) ~ "Ameiurus_melas",
    grepl("Prietella_phreatophila", branch_name_1) ~ "Prietella_phreatophila",
    grepl("Ictalurus_punctatus", branch_name_1) ~ "Ictalurus_punctatus",
    grepl("Cranoglanis_bouderius", branch_name_1) ~ "Cranoglanis_bouderius",
    grepl("Pangasianodon_hypophthalmus", branch_name_1) ~ "Pangasianodon_hypophthalmus",
    grepl("Pangasius_djambal", branch_name_1) ~ "Pangasius_djambal",
    grepl("Gogo_arcuatus", branch_name_1) ~ "Gogo_arcuatus",
    grepl("Neoarius_graeffei", branch_name_1) ~ "Neoarius_graeffei",
    grepl("Notoglanidium_macrostoma", branch_name_1) ~ "Notoglanidium_macrostoma",
    grepl("Synodontis_membranacea", branch_name_1) ~ "Synodontis_membranacea",
    grepl("Doras_micropoeus", branch_name_1) ~ "Doras_micropoeus",
    grepl("Scorpiodoras_heckelii", branch_name_1) ~ "Scorpiodoras_heckelii",
    grepl("Bagarius_yarrelli", branch_name_1) ~ "Bagarius_yarrelli",
    grepl("Hemibagrus_wyckioides", branch_name_1) ~ "Hemibagrus_wyckioides",
    grepl("Tachysurus_fulvidraco", branch_name_1) ~ "Tachysurus_fulvidraco",
    grepl("Duopalatinus_emarginatus", branch_name_1) ~ "Duopalatinus_emarginatus",
    grepl("Silurus_meridionalis", branch_name_1) ~ "Silurus_meridionalis",
    grepl("Plotosus_lineatus", branch_name_1) ~ "Plotosus_lineatus",
    grepl("Clarias_gariepinus", branch_name_1) ~ "Clarias_gariepinus",
    grepl("Ancistrus_triradiatus", branch_name_1) ~ "Ancistrus_triradiatus",
    grepl("CUL9", branch_name_1) ~ "CUL9",
    grepl("RHP1", branch_name_1) ~ "RHP1",
    grepl("Corydoras_maculifer", branch_name_1) ~ "Corydoras_maculifer",
    grepl("CHM6", branch_name_1) ~ "CHM6",
    grepl("RHH2", branch_name_1) ~ "RHH2",
    grepl("RSS1", branch_name_1) ~ "RSS1",
    grepl("Trichomycterus_rosablanca", branch_name_1) ~ "Trichomycterus_rosablanca",
    grepl("CSV83", branch_name_1) ~ "CSV83",
    grepl("CUL4", branch_name_1) ~ "CUL4",
    grepl("RUI2", branch_name_1) ~ "RUI2",
    grepl("Microcambeva_barbata", branch_name_1) ~ "Microcambeva_barbata",
    grepl("Danio_rerio", branch_name_1) ~ "Danio_rerio"
  )) %>% 
  mutate(species2 = case_when(
    grepl("Ameiurus_melas", branch_name_2) ~ "Ameiurus_melas",
    grepl("Prietella_phreatophila", branch_name_2) ~ "Prietella_phreatophila",
    grepl("Ictalurus_punctatus", branch_name_2) ~ "Ictalurus_punctatus",
    grepl("Cranoglanis_bouderius", branch_name_2) ~ "Cranoglanis_bouderius",
    grepl("Pangasianodon_hypophthalmus", branch_name_2) ~ "Pangasianodon_hypophthalmus",
    grepl("Pangasius_djambal", branch_name_2) ~ "Pangasius_djambal",
    grepl("Gogo_arcuatus", branch_name_2) ~ "Gogo_arcuatus",
    grepl("Neoarius_graeffei", branch_name_2) ~ "Neoarius_graeffei",
    grepl("Notoglanidium_macrostoma", branch_name_2) ~ "Notoglanidium_macrostoma",
    grepl("Synodontis_membranacea", branch_name_2) ~ "Synodontis_membranacea",
    grepl("Doras_micropoeus", branch_name_2) ~ "Doras_micropoeus",
    grepl("Scorpiodoras_heckelii", branch_name_2) ~ "Scorpiodoras_heckelii",
    grepl("Bagarius_yarrelli", branch_name_2) ~ "Bagarius_yarrelli",
    grepl("Hemibagrus_wyckioides", branch_name_2) ~ "Hemibagrus_wyckioides",
    grepl("Tachysurus_fulvidraco", branch_name_2) ~ "Tachysurus_fulvidraco",
    grepl("Duopalatinus_emarginatus", branch_name_2) ~ "Duopalatinus_emarginatus",
    grepl("Silurus_meridionalis", branch_name_2) ~ "Silurus_meridionalis",
    grepl("Plotosus_lineatus", branch_name_2) ~ "Plotosus_lineatus",
    grepl("Clarias_gariepinus", branch_name_2) ~ "Clarias_gariepinus",
    grepl("Ancistrus_triradiatus", branch_name_2) ~ "Ancistrus_triradiatus",
    grepl("CUL9", branch_name_2) ~ "CUL9",
    grepl("RHP1", branch_name_2) ~ "RHP1",
    grepl("Corydoras_maculifer", branch_name_2) ~ "Corydoras_maculifer",
    grepl("CHM6", branch_name_2) ~ "CHM6",
    grepl("RHH2", branch_name_2) ~ "RHH2",
    grepl("RSS1", branch_name_2) ~ "RSS1",
    grepl("Trichomycterus_rosablanca", branch_name_2) ~ "Trichomycterus_rosablanca",
    grepl("CSV83", branch_name_2) ~ "CSV83",
    grepl("CUL4", branch_name_2) ~ "CUL4",
    grepl("RUI2", branch_name_2) ~ "RUI2",
    grepl("Microcambeva_barbata", branch_name_2) ~ "Microcambeva_barbata",
    grepl("Danio_rerio", branch_name_2) ~ "Danio_rerio"
  )) 



csubst_K3.genetree <- 
  csubst_K3.genetree %>%
  mutate(species1 = case_when(
    grepl("Ameiurus_melas", branch_name_1) ~ "Ameiurus_melas",
    grepl("Prietella_phreatophila", branch_name_1) ~ "Prietella_phreatophila",
    grepl("Ictalurus_punctatus", branch_name_1) ~ "Ictalurus_punctatus",
    grepl("Cranoglanis_bouderius", branch_name_1) ~ "Cranoglanis_bouderius",
    grepl("Pangasianodon_hypophthalmus", branch_name_1) ~ "Pangasianodon_hypophthalmus",
    grepl("Pangasius_djambal", branch_name_1) ~ "Pangasius_djambal",
    grepl("Gogo_arcuatus", branch_name_1) ~ "Gogo_arcuatus",
    grepl("Neoarius_graeffei", branch_name_1) ~ "Neoarius_graeffei",
    grepl("Notoglanidium_macrostoma", branch_name_1) ~ "Notoglanidium_macrostoma",
    grepl("Synodontis_membranacea", branch_name_1) ~ "Synodontis_membranacea",
    grepl("Doras_micropoeus", branch_name_1) ~ "Doras_micropoeus",
    grepl("Scorpiodoras_heckelii", branch_name_1) ~ "Scorpiodoras_heckelii",
    grepl("Bagarius_yarrelli", branch_name_1) ~ "Bagarius_yarrelli",
    grepl("Hemibagrus_wyckioides", branch_name_1) ~ "Hemibagrus_wyckioides",
    grepl("Tachysurus_fulvidraco", branch_name_1) ~ "Tachysurus_fulvidraco",
    grepl("Duopalatinus_emarginatus", branch_name_1) ~ "Duopalatinus_emarginatus",
    grepl("Silurus_meridionalis", branch_name_1) ~ "Silurus_meridionalis",
    grepl("Plotosus_lineatus", branch_name_1) ~ "Plotosus_lineatus",
    grepl("Clarias_gariepinus", branch_name_1) ~ "Clarias_gariepinus",
    grepl("Ancistrus_triradiatus", branch_name_1) ~ "Ancistrus_triradiatus",
    grepl("CUL9", branch_name_1) ~ "CUL9",
    grepl("RHP1", branch_name_1) ~ "RHP1",
    grepl("Corydoras_maculifer", branch_name_1) ~ "Corydoras_maculifer",
    grepl("CHM6", branch_name_1) ~ "CHM6",
    grepl("RHH2", branch_name_1) ~ "RHH2",
    grepl("RSS1", branch_name_1) ~ "RSS1",
    grepl("Trichomycterus_rosablanca", branch_name_1) ~ "Trichomycterus_rosablanca",
    grepl("CSV83", branch_name_1) ~ "CSV83",
    grepl("CUL4", branch_name_1) ~ "CUL4",
    grepl("RUI2", branch_name_1) ~ "RUI2",
    grepl("Microcambeva_barbata", branch_name_1) ~ "Microcambeva_barbata",
    grepl("Danio_rerio", branch_name_1) ~ "Danio_rerio"
  )) %>% 
  mutate(species2 = case_when(
    grepl("Ameiurus_melas", branch_name_2) ~ "Ameiurus_melas",
    grepl("Prietella_phreatophila", branch_name_2) ~ "Prietella_phreatophila",
    grepl("Ictalurus_punctatus", branch_name_2) ~ "Ictalurus_punctatus",
    grepl("Cranoglanis_bouderius", branch_name_2) ~ "Cranoglanis_bouderius",
    grepl("Pangasianodon_hypophthalmus", branch_name_2) ~ "Pangasianodon_hypophthalmus",
    grepl("Pangasius_djambal", branch_name_2) ~ "Pangasius_djambal",
    grepl("Gogo_arcuatus", branch_name_2) ~ "Gogo_arcuatus",
    grepl("Neoarius_graeffei", branch_name_2) ~ "Neoarius_graeffei",
    grepl("Notoglanidium_macrostoma", branch_name_2) ~ "Notoglanidium_macrostoma",
    grepl("Synodontis_membranacea", branch_name_2) ~ "Synodontis_membranacea",
    grepl("Doras_micropoeus", branch_name_2) ~ "Doras_micropoeus",
    grepl("Scorpiodoras_heckelii", branch_name_2) ~ "Scorpiodoras_heckelii",
    grepl("Bagarius_yarrelli", branch_name_2) ~ "Bagarius_yarrelli",
    grepl("Hemibagrus_wyckioides", branch_name_2) ~ "Hemibagrus_wyckioides",
    grepl("Tachysurus_fulvidraco", branch_name_2) ~ "Tachysurus_fulvidraco",
    grepl("Duopalatinus_emarginatus", branch_name_2) ~ "Duopalatinus_emarginatus",
    grepl("Silurus_meridionalis", branch_name_2) ~ "Silurus_meridionalis",
    grepl("Plotosus_lineatus", branch_name_2) ~ "Plotosus_lineatus",
    grepl("Clarias_gariepinus", branch_name_2) ~ "Clarias_gariepinus",
    grepl("Ancistrus_triradiatus", branch_name_2) ~ "Ancistrus_triradiatus",
    grepl("CUL9", branch_name_2) ~ "CUL9",
    grepl("RHP1", branch_name_2) ~ "RHP1",
    grepl("Corydoras_maculifer", branch_name_2) ~ "Corydoras_maculifer",
    grepl("CHM6", branch_name_2) ~ "CHM6",
    grepl("RHH2", branch_name_2) ~ "RHH2",
    grepl("RSS1", branch_name_2) ~ "RSS1",
    grepl("Trichomycterus_rosablanca", branch_name_2) ~ "Trichomycterus_rosablanca",
    grepl("CSV83", branch_name_2) ~ "CSV83",
    grepl("CUL4", branch_name_2) ~ "CUL4",
    grepl("RUI2", branch_name_2) ~ "RUI2",
    grepl("Microcambeva_barbata", branch_name_2) ~ "Microcambeva_barbata",
    grepl("Danio_rerio", branch_name_2) ~ "Danio_rerio"
  )) %>%
  mutate(species3 = case_when(
    grepl("Ameiurus_melas", branch_name_3) ~ "Ameiurus_melas",
    grepl("Prietella_phreatophila", branch_name_3) ~ "Prietella_phreatophila",
    grepl("Ictalurus_punctatus", branch_name_3) ~ "Ictalurus_punctatus",
    grepl("Cranoglanis_bouderius", branch_name_3) ~ "Cranoglanis_bouderius",
    grepl("Pangasianodon_hypophthalmus", branch_name_3) ~ "Pangasianodon_hypophthalmus",
    grepl("Pangasius_djambal", branch_name_3) ~ "Pangasius_djambal",
    grepl("Gogo_arcuatus", branch_name_3) ~ "Gogo_arcuatus",
    grepl("Neoarius_graeffei", branch_name_3) ~ "Neoarius_graeffei",
    grepl("Notoglanidium_macrostoma", branch_name_3) ~ "Notoglanidium_macrostoma",
    grepl("Synodontis_membranacea", branch_name_3) ~ "Synodontis_membranacea",
    grepl("Doras_micropoeus", branch_name_3) ~ "Doras_micropoeus",
    grepl("Scorpiodoras_heckelii", branch_name_3) ~ "Scorpiodoras_heckelii",
    grepl("Bagarius_yarrelli", branch_name_3) ~ "Bagarius_yarrelli",
    grepl("Hemibagrus_wyckioides", branch_name_3) ~ "Hemibagrus_wyckioides",
    grepl("Tachysurus_fulvidraco", branch_name_3) ~ "Tachysurus_fulvidraco",
    grepl("Duopalatinus_emarginatus", branch_name_3) ~ "Duopalatinus_emarginatus",
    grepl("Silurus_meridionalis", branch_name_3) ~ "Silurus_meridionalis",
    grepl("Plotosus_lineatus", branch_name_3) ~ "Plotosus_lineatus",
    grepl("Clarias_gariepinus", branch_name_3) ~ "Clarias_gariepinus",
    grepl("Ancistrus_triradiatus", branch_name_3) ~ "Ancistrus_triradiatus",
    grepl("CUL9", branch_name_3) ~ "CUL9",
    grepl("RHP1", branch_name_3) ~ "RHP1",
    grepl("Corydoras_maculifer", branch_name_3) ~ "Corydoras_maculifer",
    grepl("CHM6", branch_name_3) ~ "CHM6",
    grepl("RHH2", branch_name_3) ~ "RHH2",
    grepl("RSS1", branch_name_3) ~ "RSS1",
    grepl("Trichomycterus_rosablanca", branch_name_3) ~ "Trichomycterus_rosablanca",
    grepl("CSV83", branch_name_3) ~ "CSV83",
    grepl("CUL4", branch_name_3) ~ "CUL4",
    grepl("RUI2", branch_name_3) ~ "RUI2",
    grepl("Microcambeva_barbata", branch_name_3) ~ "Microcambeva_barbata",
    grepl("Danio_rerio", branch_name_3) ~ "Danio_rerio"
  )) 




#Remove if Inf values = computation not possible
#Filter only significant results

csubst_K2.genetree <- 
  csubst_K2.genetree %>%
  #filter(omegaCany2spe > 5.0) %>%
  #filter(OCNany2spe > 2.0) %>%
  filter(omegaCany2spe != Inf)  %>%
  filter(OCNany2spe != Inf) 

csubst_K3.genetree <- 
  csubst_K3.genetree %>%
  #filter(omegaCany2spe > 5.0) %>%
  #filter(OCNany2spe > 2.0) %>%
  filter(omegaCany2spe != Inf)  %>%
  filter(OCNany2spe != Inf) 


#Add a column to say if the result is significant 

csubst_K2.genetree <- 
  csubst_K2.genetree %>%
  mutate(significant = if_else(
    omegaCany2spe > 5.0 & OCNany2spe > 2.0,
    "yes",
    "no"
  ))


csubst_K3.genetree <- 
  csubst_K3.genetree %>%
  mutate(significant = if_else(
    omegaCany2spe > 5.0 & OCNany2spe > 2.0,
    "yes",
    "no"
  ))


##### Analyze cSUBST results -- Gene Tree   ---------------------------------

#Do cave/cave have more convergent proteins than surface/surface and surface/cave ?

species_list_silu <- species_tree_len_woDr$tip.label

uniq_comb <- combn(species_list_silu,2)
uniq_comb <- data.frame(species1 = uniq_comb[1, ], species2 = uniq_comb[2, ])
uniq_comb <- left_join(uniq_comb, species_df, by=c("species1" = "species"))
uniq_comb <- left_join(uniq_comb, species_df, by=c("species2" = "species"))
colnames(uniq_comb) <- c("species1", "species2", "Habitat_1", "Habitat_2")


summary_pairs_K2_df <- as.data.frame(NULL)
for(curr_comb in 1:nrow(uniq_comb)){
  curr_comb_l <- uniq_comb[curr_comb,]
  curr_sp1 <- uniq_comb[curr_comb,]$species1
  curr_sp2 <- uniq_comb[curr_comb,]$species2
  
  curr_K2 <- 
    csubst_K2.genetree %>% 
    filter((species1 == curr_sp1 & species2 == curr_sp2) | (species1 == curr_sp2 & species2 == curr_sp1))
  
  non_signif_nb <- nrow(curr_K2 %>% filter(significant == "no"))
  signif_nb <- nrow(curr_K2 %>% filter(significant == "yes"))
  
  curr_df <- cbind(curr_comb_l, signif_nb, non_signif_nb)
  
  summary_pairs_K2_df <- rbind(summary_pairs_K2_df, curr_df)
}

is.numeric(summary_pairs_K2_df$signif_nb)
is.numeric(summary_pairs_K2_df$non_signif_nb)

summary_pairs_K2_df <- 
  summary_pairs_K2_df %>%
  mutate(prop_signif = signif_nb/(signif_nb + non_signif_nb))


summary_pairs_K2_df <- 
  summary_pairs_K2_df %>%
  mutate(comp_habitat = case_when(
    Habitat_1 == "Surface" & Habitat_2 == "Surface" ~ "S_S",
    Habitat_1 == "Surface" & Habitat_2 == "Cave" ~ "S_C",
    Habitat_1 == "Cave" & Habitat_2 == "Surface" ~ "S_C",
    Habitat_1 == "Cave" & Habitat_2 == "Cave" ~ "C_C"
  ))


summary_pairs_K2_df %>%
  filter(comp_habitat == "C_C") %>%
  arrange(prop_signif) 

#look at which genes are convergent between cave species




#Compare csubst between environments at the genome-wide level

summary_pairs_K2_df %>%
  ggplot(., aes(x=comp_habitat, y=prop_signif)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  ylab("Proportion of convergent genes") +
  xlab("Habitat comparison")


#Compare csubst between environments per go-terms


GO_annots_df <- read.table("OGG_GOterms.tsv",sep="\t",header=FALSE)
colnames(GO_annots_df) <- c("OGG", "GO")
GO_classes_df <- read.table("goterm_class.csv",sep="\t",header=FALSE)
colnames(GO_classes_df) <- c("GO", "class")
GO_annots_df <- left_join(GO_annots_df, GO_classes_df, by="GO")
BP_GO_annots_df <- GO_annots_df %>% filter(class == "biological_process") %>% dplyr::select(-class)
BP_GO_annots_list <- split(BP_GO_annots_df$GO, BP_GO_annots_df$OGG)


list_goterms <- BP_GO_annots_df %>% pull(GO) %>% unique()
K2_wilcox_GO_df <- as.data.frame(NULL)
for(curr_GO in list_goterms){
  
  curr_OGGs <- BP_GO_annots_df %>% filter(GO == curr_GO) %>% pull(OGG)
  
  summary_pairs_K2_df_go <- as.data.frame(NULL)
  for(curr_comb in 1:nrow(uniq_comb)){
    curr_comb_l <- uniq_comb[curr_comb,]
    curr_sp1 <- uniq_comb[curr_comb,]$species1
    curr_sp2 <- uniq_comb[curr_comb,]$species2
    
    curr_K2 <- 
      csubst_K2.genetree %>% 
      filter(OGG %in% curr_OGGs) %>%
      filter((species1 == curr_sp1 & species2 == curr_sp2) | (species1 == curr_sp2 & species2 == curr_sp1))
    
    non_signif_nb <- nrow(curr_K2 %>% filter(significant == "no"))
    signif_nb <- nrow(curr_K2 %>% filter(significant == "yes"))
    total_nb <- non_signif_nb + signif_nb
    #keep if >100 genes are computed
    if(total_nb >= 100){
      curr_df <- cbind(curr_comb_l, signif_nb, non_signif_nb)
      summary_pairs_K2_df_go <- rbind(summary_pairs_K2_df_go, curr_df)
    }
  }
  
  if(nrow(summary_pairs_K2_df_go) >= 6){
    
    summary_pairs_K2_df_go <- summary_pairs_K2_df_go %>% mutate(prop_signif = signif_nb/(signif_nb + non_signif_nb))
   
    summary_pairs_K2_df_go <- 
      summary_pairs_K2_df_go %>%
      mutate(comp_habitat = case_when(
        Habitat_1 == "Surface" & Habitat_2 == "Surface" ~ "S_S",
        Habitat_1 == "Surface" & Habitat_2 == "Cave" ~ "S_C",
        Habitat_1 == "Cave" & Habitat_2 == "Surface" ~ "S_C",
        Habitat_1 == "Cave" & Habitat_2 == "Cave" ~ "C_C"
      ))
    
    nb_CC <- summary_pairs_K2_df_go %>% filter(comp_habitat == "C_C")
    nb_SC <- summary_pairs_K2_df_go %>% filter(comp_habitat == "S_C")
    nb_SS <- summary_pairs_K2_df_go %>% filter(comp_habitat == "S_S")
    
    if(nb_SC >= 3 & nb_SS >= 3 & nb_CC >= 3){
      
      mean_S_S <- summary_pairs_K2_df_go %>% filter(comp_habitat == "S_S") %>% pull(prop_signif) %>% mean()
      mean_S_C <- summary_pairs_K2_df_go %>% filter(comp_habitat == "S_C") %>% pull(prop_signif) %>% mean()
      mean_C_C <- summary_pairs_K2_df_go %>% filter(comp_habitat == "C_C") %>% pull(prop_signif) %>% mean()
      
      C_C_vs_S_C <-
        wilcox.test(summary_pairs_K2_df_go %>% filter(comp_habitat == "C_C") %>% pull(prop_signif),
                    summary_pairs_K2_df_go %>% filter(comp_habitat == "S_C") %>% pull(prop_signif))
      P_C_C_vs_S_C <- C_C_vs_S_C$p.value
      
      
      C_C_vs_S_S <-
        wilcox.test(summary_pairs_K2_df_go %>% filter(comp_habitat == "C_C") %>% pull(prop_signif),
                    summary_pairs_K2_df_go %>% filter(comp_habitat == "S_S") %>% pull(prop_signif))
      P_C_C_vs_S_S <- C_C_vs_S_S$p.value
      
      curr_df_wilcox <- as.data.frame(cbind(curr_GO, mean_S_S, mean_S_C, mean_C_C, P_C_C_vs_S_C, P_C_C_vs_S_S))
      
      K2_wilcox_GO_df <- rbind(K2_wilcox_GO_df, curr_df_wilcox)
      
      
    }
    
  }

  
}

K2_wilcox_GO_df %>%
  filter(mean_C_C > mean_S_S)




# Look at which genes are convergent between cave species

species_df1 <- species_df
species_df2 <- species_df

colnames(species_df1) <- c("species1", "Habitat1")
colnames(species_df2) <- c("species2", "Habitat2")
csubst_K2.genetree <- left_join(csubst_K2.genetree, species_df1, by="species1")
csubst_K2.genetree <- left_join(csubst_K2.genetree, species_df2, by="species2")

cave_csubst_K2.genetree <-
  csubst_K2.genetree %>%
  filter(Habitat1 == "Cave" & Habitat2 == "Cave") %>%
  filter(significant == "yes") %>%
  filter(species1 != species2)


colnames(OGG_names_df) <- c("OGG", "gene_ID", "gene_name", "evalue_name")
cave_csubst_K2.genetree <- left_join(cave_csubst_K2.genetree, OGG_names_df, by="OGG")



cave_csubst_K2.genetree_dist <- 
  cave_csubst_K2.genetree %>% 
  dplyr::select(OGG, species1, species2, gene_ID, gene_name) %>%
  distinct()

signif_OGG <- cave_csubst_K2.genetree %>% pull(OGG) %>% unique()

signif_names <- cave_csubst_K2.genetree %>% pull(gene_name) %>% unique()



## is there an enrichment ?

#Define the gene universe and significant genes
background_gene_list <- csubst_K2.genetree %>% pull(OGG) %>% unique()
my_geneList <- factor(as.integer(background_gene_list %in% signif_OGG))
names(my_geneList) <- background_gene_list


GOdata_BP <- new("topGOdata", ontology = "BP", 
                 allGenes = my_geneList, 
                 annot = annFUN.gene2GO,
                 gene2GO = BP_GO_annots_list)


result_enrich_fisher_BP <- runTest(GOdata_BP, algorithm = "classic", statistic = "fisher")

adjusted = p.adjust(score(result_enrich_fisher_BP), method = "BH")
adjustedP = data.frame(GO.ID = names(adjusted), correctedFisher = adjusted)
colnames(adjustedP) <- c("GO", "adj.pvalue")

my_results_BP <- 
  GenTable(GOdata_BP, classicFisher = result_enrich_fisher_BP, orderBy = "classicFisher",
                          ranksOf = "classicFisher", topNodes = 60)


colnames(my_results_BP) <- c("GO", "term", "annotated", "significant", "expected", "pvalue")


my_results_BP_adj <- left_join(my_results_BP, adjustedP, by="GO")


##### Import cSUBST results -- Species Tree   ---------------------------------

#Import branch informations

csubst_b.speciestree <- 
  as.data.frame(
    fread("all_OGG_csubst_b.SpeciesTree.tsv",
          header=TRUE,
          sep="\t")
  )



#Import K = 2 results


csubst_K2.speciestree <- 
  as.data.frame(
    fread("all_OGG_csubst_cb_2.SpeciesTree.tsv",
          header=TRUE,
          sep="\t")
  )



#Import K = 3 results

csubst_K3.speciestree <- 
  as.data.frame(
    fread("all_OGG_csubst_cb_3.SpeciesTree.tsv",
          header=TRUE,
          sep="\t")
  ) #empty dataframe - No conergent at k=3




#Add branch names to K2 and K3 tables


#Species Tree
branch_id_name_df <- csubst_b.speciestree  %>% dplyr::select(OGG, branch_name, branch_id)
colnames(branch_id_name_df) <- c("OGG", "branch_name_1", "branch_id_1")
csubst_K2.speciestree <- left_join(csubst_K2.speciestree, branch_id_name_df, by=c("OGG", "branch_id_1"))
colnames(branch_id_name_df) <- c("OGG", "branch_name_2", "branch_id_2")
csubst_K2.speciestree <- left_join(csubst_K2.speciestree, branch_id_name_df, by=c("OGG", "branch_id_2"))


#For CSUBST on gene trees, remove node values

csubst_K2.speciestree <- 
  csubst_K2.speciestree %>%
  filter(! grepl("Node", branch_name_1)) %>%
  filter(! grepl("Node", branch_name_2))




#Add the species names

csubst_K2.speciestree <- 
  csubst_K2.speciestree %>%
  mutate(species1 = branch_name_1) %>% 
  mutate(species2 = branch_name_2) 




#Remove if Inf values = computation not possible
#Filter only significant results

csubst_K2.speciestree <- 
  csubst_K2.speciestree %>%
  #filter(omegaCany2spe > 5.0) %>%
  #filter(OCNany2spe > 2.0) %>%
  filter(omegaCany2spe != Inf)  %>%
  filter(OCNany2spe != Inf) 


#Add a column to say if the result is significant 

csubst_K2.speciestree <- 
  csubst_K2.speciestree %>%
  mutate(significant = if_else(
    omegaCany2spe > 5.0 & OCNany2spe > 2.0,
    "yes",
    "no"
  ))





##### Analyze cSUBST results -- Species Tree   ---------------------------------

#Do cave/cave have more convergent proteins than surface/surface and surface/cave ?


uniq_comb <- combn(species_list_silu,2)
uniq_comb <- data.frame(species1 = uniq_comb[1, ], species2 = uniq_comb[2, ])
uniq_comb <- left_join(uniq_comb, species_df, by=c("species1" = "species"))
uniq_comb <- left_join(uniq_comb, species_df, by=c("species2" = "species"))
colnames(uniq_comb) <- c("species1", "species2", "Habitat_1", "Habitat_2")


summary_pairs_K2_df_sp <- as.data.frame(NULL)
for(curr_comb in 1:nrow(uniq_comb)){
  curr_comb_l <- uniq_comb[curr_comb,]
  curr_sp1 <- uniq_comb[curr_comb,]$species1
  curr_sp2 <- uniq_comb[curr_comb,]$species2
  
  curr_K2 <- 
    csubst_K2.speciestree %>% 
    filter((species1 == curr_sp1 & species2 == curr_sp2) | (species1 == curr_sp2 & species2 == curr_sp1))
  
  non_signif_nb <- nrow(curr_K2 %>% filter(significant == "no"))
  signif_nb <- nrow(curr_K2 %>% filter(significant == "yes"))
  
  curr_df <- cbind(curr_comb_l, signif_nb, non_signif_nb)
  
  summary_pairs_K2_df_sp <- rbind(summary_pairs_K2_df_sp, curr_df)
}

is.numeric(summary_pairs_K2_df_sp$signif_nb)
is.numeric(summary_pairs_K2_df_sp$non_signif_nb)

summary_pairs_K2_df_sp <- 
  summary_pairs_K2_df_sp %>%
  mutate(prop_signif = signif_nb/(signif_nb + non_signif_nb))


summary_pairs_K2_df_sp <- 
  summary_pairs_K2_df_sp %>%
  mutate(comp_habitat = case_when(
    Habitat_1 == "Surface" & Habitat_2 == "Surface" ~ "S_S",
    Habitat_1 == "Surface" & Habitat_2 == "Cave" ~ "S_C",
    Habitat_1 == "Cave" & Habitat_2 == "Surface" ~ "S_C",
    Habitat_1 == "Cave" & Habitat_2 == "Cave" ~ "C_C"
  ))



#Compare csubst between environments at the genome-wide level

summary_pairs_K2_df_sp %>%
  ggplot(., aes(x=comp_habitat, y=prop_signif)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        axis.title.x=element_blank()) +
  theme(legend.position="none") +
  ylab("Proportion of convergent genes") +
  xlab("Habitat comparison")

