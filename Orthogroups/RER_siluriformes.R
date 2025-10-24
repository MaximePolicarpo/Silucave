##### Libraries  ---------------------------------

rm(list=ls())

set.seed(2712)


library("caper")
library("ape")
library(dplyr)
library(ggplot2)
library(tidyverse)
library(data.table)
library(phytools)
library(purrr)
library("RERconverge")
library(ggtree)
library(ggtreeExtra)
split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}




#### My functions + Colors palettes  ---------------------------------

args = commandArgs(trailingOnly=TRUE)


PGLS_pvalue <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  pvalue =  formatC(sum_cor$coefficients[8], digits = 3)
  if (pvalue == "   0"){ pvalue = 2.2e-16}
  return(pvalue)
}

PGLS_R2 <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  r2 = formatC(sum_cor$r.squared, digits = 2)
  return(r2)
}

PGLS_lambda <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  PGLS_lambda = sum_cor$param[2]
  return(PGLS_lambda)
}


GLS_function <- function(gls_rslt) {
  GLS_cc <- coef(gls_rslt)
  GLS_f <- function(x) GLS_cc[1] + GLS_cc[2]*x
  return(GLS_f)
}


##### Data load -  GO terms   -------------------

GO_annots.BP <- read.gmt("goterms_annot.gmt.BP")
GO_annots.CC <- read.gmt("goterms_annot.gmt.CP")
GO_annots.MF <- read.gmt("goterms_annot.gmt.MF")
Panther_annots <-  read.gmt("panther_annot.gmt")
annots.KEGG <- read.gmt("KEGG.all.HOG.gmt")

GO_annots_list <- list(GO_annots.BP, GO_annots.CC, GO_annots.MF, Panther_annots, annots.KEGG)
names(GO_annots_list) <- 
  c("biological_process", "cellular_component", "molecular_function", "PANTHER","KEGG")




##### Data load - Species trees   ---------------------------------

species_tree <- 
  read.tree("AMAS_concatenated_alignment_BUSCO.fa.timetree.nwk")

species_list <- species_tree$tip.label

##### Read trees and form a MasterTree   ---------------------------------

GW_trees.all <- 
  RERconverge::readTrees(
    "AllTrees.Phangorn.nwk",
    reestimateBranches = FALSE,
    minSpecs=6)


saveRDS(GW_trees.all, file ="AllTrees.Phangorn.RDS")


##### RERconverge -- Binary  ---------------------------------

RERw.all <-
 getAllResiduals(GW_trees.all,
                 useSpecies=GW_trees.all$masterTree$tip.label,
                 transform = "sqrt",
                 weighted = T,
                 scale = T,
                 min.sp = 6)



saveRDS(RERw.all, file="RERw.all.RDS")

## Generate a binary trait tree

cavefish_species <- 
  c("Prietella_phreatophila",
    "CHM6",
    "CSV83",
    "CUL4",
    "CUL9",
    "Trichomycterus_rosablanca")


all.anc <- 
 foreground2Tree(
   treesObj = GW_trees.all,
   foreground = cavefish_species,
   clade = "ancestral")

all.all <- 
 foreground2Tree(
   treesObj = GW_trees.all,
   foreground = cavefish_species,
   clade = "all")

all.term <- 
 foreground2Tree(
   treesObj = GW_trees.all,
   foreground = cavefish_species,
   clade = "terminal")

dev.off()

## Generate paths

all.anc.phenv <- tree2Paths(all.anc, GW_trees.all)
all.all.phenv <- tree2Paths(all.all, GW_trees.all)
all.term.phenv <- tree2Paths(all.term, GW_trees.all)

## Correlate gene evolution with our binary trait 

all.anc.res <- 
  correlateWithBinaryPhenotype(RERw.all, all.anc.phenv, min.sp=6, min.pos=2, weighted="auto")
all.anc.res$HOG <- row.names(all.anc.res)
all.anc.res$stat = -log10(all.anc.res$P) * sign(all.anc.res$Rho)
write.table(all.anc.res, "all.anc.res", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)


all.all.res <- 
  correlateWithBinaryPhenotype(RERw.all, all.all.phenv, min.sp=6, min.pos=2, weighted="auto")
all.all.res$HOG <- row.names(all.all.res)
all.all.res$stat = -log10(all.all.res$P) * sign(all.all.res$Rho)
write.table(all.all.res, "all.all.res", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

all.term.res <- 
  correlateWithBinaryPhenotype(RERw.all, all.term.phenv, min.sp=6, min.pos=2, weighted="auto")
all.term.res$HOG <- row.names(all.term.res)
all.term.res$stat = -log10(all.term.res$P) * sign(all.term.res$Rho)
write.table(all.all.res, "all.term.res", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)


print("RERconverge in binary mode ran successfully")


print("Launching binary permulation analysis")



#Find which species are the outgroup in the master tree.

curr_pruned_tree <- GW_trees.all$masterTree

root_edge <- curr_pruned_tree$edge[curr_pruned_tree$edge[, 1] == Ntip(curr_pruned_tree) + 1,]
outgroup_indice_1 <- root_edge[, 2][1]
outgroup_indice_2 <- root_edge[, 2][2]
ntip_ind1 <- Ntip(extract.clade(curr_pruned_tree, outgroup_indice_1))
ntip_ind2 <- Ntip(extract.clade(curr_pruned_tree, outgroup_indice_2))
if(ntip_ind1 <= ntip_ind2){
 good_outgroup_ind <- outgroup_indice_1
} else {
 good_outgroup_ind <- outgroup_indice_2
}
outgroup_species <- extract.clade(curr_pruned_tree, good_outgroup_ind)$tip.label
root_MRCA <- getMRCA(curr_pruned_tree, outgroup_species)
rooted_master_tree <- root(curr_pruned_tree, node=root_MRCA, resolve.root= TRUE)


#Generate a list of species which have the phenotype of interest
tip_labels_df <- as.data.frame(rooted_master_tree$tip.label)
tip_labels_df <- as.data.frame(tip_labels_df)
colnames(tip_labels_df) <- "species"

#Generate the list of monophyletic clades with >= 2 species sharing the phenotype of interest
monophyletic_clades <- list()
for (node in 1:rooted_master_tree$Nnode + length(rooted_master_tree$tip.label)) {
 curr_clade <- extract.clade(rooted_master_tree, node = node)
 curr_species <- curr_clade$tip.label
 
 
 nb_non_pheno <- 
   nrow(
     tip_labels_df %>%
       filter(species %in% curr_species) %>%
       filter(! species %in% cavefish_species)
   )
 
 
 if (nb_non_pheno == 0) {
   monophyletic_clades[[length(monophyletic_clades) + 1]] <- curr_clade$tip.label
 }
}

non_redundant_clades <- list()
for (i in seq_along(monophyletic_clades)) {
 current_clade <- monophyletic_clades[[i]]
 is_redundant <- FALSE
 
 for (j in seq_along(monophyletic_clades)) {
   if (i != j) {
     other_clade <- monophyletic_clades[[j]]
     if (all(current_clade %in% other_clade)) {
       is_redundant <- TRUE
       break
     }
   }
 }
 
 if (!is_redundant) {
   non_redundant_clades[[length(non_redundant_clades) + 1]] <- current_clade
 }
}

#Perform enrichment analysis on the raw results

my_stat.anc <- sign(all.anc.res$Rho) * (-log10(all.anc.res$P))
names(my_stat.anc) <- rownames(all.anc.res)
genenames <- names(my_stat.anc)
my_stat.anc <- my_stat.anc[!is.na(my_stat.anc)]
my_enrichment.anc <- fastwilcoxGMTall(my_stat.anc, GO_annots_list, outputGeneVals=T, num.g=5)

saveRDS(my_enrichment.anc, file="enrichment_all_raw.RDS")

my_stat.term <- sign(all.term.res$Rho) * (-log10(all.term.res$P))
names(my_stat.term) <- rownames(all.term.res)
genenames <- names(my_stat.term)
my_stat.term <- my_stat.term[!is.na(my_stat.term)]
my_enrichment.term <- fastwilcoxGMTall(my_stat.term, GO_annots_list, outputGeneVals=T, num.g=5)

saveRDS(my_enrichment.anc, file="enrichment_all_raw.RDS")



my_stat.all <- sign(all.all.res$Rho) * (-log10(all.all.res$P))
names(my_stat.all) <- rownames(all.all.res)
genenames <- names(my_stat.all)
my_stat.all <- my_stat.all[!is.na(my_stat.all)]
my_enrichment.all <- fastwilcoxGMTall(my_stat.all, GO_annots_list, outputGeneVals=T, num.g=5)


saveRDS(my_enrichment.all, file="enrichment_all_raw.RDS")


#Lets launch 1,000 permulations with the CC mode ! 

my_perms_binary <- 
  getPermsBinary(numperms = 1000, 
                 fg_vec = cavefish_species, 
                 root_sp= outgroup_species,
                 sisters_list = non_redundant_clades, 
                 RERmat = RERw.all, 
                 trees = GW_trees.all, 
                 mastertree = rooted_master_tree,
                 permmode="cc",
                 calculateenrich = F,
                 annotlist=GO_annots_list)

#Compute permulations enrichment statistics and p-values
permswithenrich.term <- getEnrichPerms(my_perms_binary, my_enrichment.term, GO_annots_list)
permswithenrich.anc <- getEnrichPerms(my_perms_binary, my_enrichment.anc, GO_annots_list)
permswithenrich.all <- getEnrichPerms(my_perms_binary, my_enrichment.all, GO_annots_list)

#Compute permulations  p-values
permpvalCC.term <- permpvalcor(all.term.res, my_perms_binary)
permpvalCC.anc <- permpvalcor(all.anc.res, my_perms_binary)
permpvalCC.all <- permpvalcor(all.all.res, my_perms_binary)

#Compute permulations enrichment p-values
enrichpermpvals.term <- permpvalenrich(my_enrichment.term, permswithenrich.term)
enrichpermpvals.anc <- permpvalenrich(my_enrichment.anc, permswithenrich.anc)
enrichpermpvals.all <- permpvalenrich(my_enrichment.all, permswithenrich.all)

#Add the enrichment permulations p-values to the real enrichment list

count=1
while(count<=length(my_enrichment.all)){
  my_enrichment.all[[count]]$permpval <- enrichpermpvals.all[[count]][match(rownames(my_enrichment.all[[count]]),names(enrichpermpvals.all[[count]]))]
  my_enrichment.all[[count]]$permpvaladj <- p.adjust(my_enrichment.all[[count]]$permpval, method="BH")
 count=count+1 
}

count=1
while(count<=length(my_enrichment.term)){
  my_enrichment.term[[count]]$permpval <- enrichpermpvals.term[[count]][match(rownames(my_enrichment.term[[count]]),names(enrichpermpvals.term[[count]]))]
  my_enrichment.term[[count]]$permpvaladj <- p.adjust(my_enrichment.term[[count]]$permpval, method="BH")
  count=count+1 
}


count=1
while(count<=length(my_enrichment.anc)){
  my_enrichment.anc[[count]]$permpval <- enrichpermpvals.anc[[count]][match(rownames(my_enrichment.anc[[count]]),names(enrichpermpvals.anc[[count]]))]
  my_enrichment.anc[[count]]$permpvaladj <- p.adjust(my_enrichment.anc[[count]]$permpval, method="BH")
  count=count+1 
}

#Add the permulations p-values and statistics (just -log10(p-values) * sign of Rho) to the raw results 

all.all.res.perm <- 
  left_join(
    rownames_to_column(all.all.res), rownames_to_column(permpvalCC.all), by="rowname") %>%
  column_to_rownames(var = "rowname")

all.term.res.perm <- 
 left_join(
   rownames_to_column(all.term.res), rownames_to_column(permpvalCC.term), by="rowname") %>%
 column_to_rownames(var = "rowname")

all.anc.res.perm <- 
  left_join(
    rownames_to_column(all.anc.res), rownames_to_column(permpvalCC.anc), by="rowname") %>%
  column_to_rownames(var = "rowname")



#Compute the adjusted p-value :) 
all.all.res.perm$permpvaladj <- p.adjust(all.all.res.perm$permpval, method="BH")
all.anc.res.perm$permpvaladj <- p.adjust(all.anc.res.perm$permpval, method="BH")
all.term.res.perm$permpvaladj <- p.adjust(all.term.res.perm$permpval, method="BH")

#write the data.frame to file , as well as the enrichment list


write.table(all.all.res.perm, "all.all.res.perm", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(all.anc.res.perm, "all.anc.res.perm", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(all.term.res.perm, "all.term.res.perm", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)


saveRDS(my_enrichment.all, file="enrichment.all.perm.RDS")
saveRDS(my_enrichment.term, file="enrichment.term.perm.RDS")
saveRDS(my_enrichment.anc, file="enrichment.anc.perm.RDS")

print("Binary permulation statistics successfully computed !")


