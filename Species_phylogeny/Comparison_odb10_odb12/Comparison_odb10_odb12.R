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

##### Import dated trees  ---------------------------------

odb10_dated_tree <- read.nexus("ODB10.AMAS_concatenated_alignment_BUSCO.fa.timetree.nex")
odb12_dated_tree <- read.nexus("ODB12.AMAS_concatenated_alignment_BUSCO.fa.timetree.nex")

odb10_dated_tree <- makeNodeLabel(odb10_dated_tree, method = "number", prefix = "Node")
odb12_dated_tree <- makeNodeLabel(odb12_dated_tree, method = "number", prefix = "Node")


odb10_data <- ggtree(odb10_dated_tree)$data 
odb10_data <- as.data.frame(odb10_data %>% dplyr::select(label, branch.length))

odb12_data <- ggtree(odb12_dated_tree)$data 
odb12_data <- as.data.frame(odb12_data %>% dplyr::select(label, branch.length))

combined_data <- left_join(odb10_data, odb12_data, by="label") %>% filter(label != "Danio_rerio")

cor.test(combined_data$branch.length.x, combined_data$branch.length.y, method="pearson")
lm_date <- lm(branch.length.y ~ branch.length.x, data=combined_data)
summary(lm_date)


combined_data %>%
  ggplot(., aes(x=branch.length.x, y=branch.length.y)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_classic() +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  xlab("odb10 - Branch lengths") +
  ylab("odb12 - Branch lengths") 


#### Compare dN/dS dating results  ---------------------------------

neutral.dNdS <- 1

## Prietella_phreatophila

Prietella_dNdS_bs <- 
  read.table("~/SiluCave/Meredith_method/Prietella_phreatophila/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(Prietella_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

Prietella_dNdS_bs.ODB10 <- Prietella_dNdS_bs
Prietella_dNdS_bs.ODB12 <- Prietella_dNdS_bs



DivergenceTime_P_A.ODB10 <- 12.9362 * 1000000
DivergenceTime_P_A.ODB12 <- 11.6734 * 1000000


Prietella_dNdS_bs.ODB10 <- 
  Prietella_dNdS_bs.ODB10 %>%
  mutate(Tp = 
           DivergenceTime_P_A.ODB10 - ((DivergenceTime_P_A.ODB10 * (Cave_dNdS - neutral.dNdS)) / 
                                         (Surface_dNdS - neutral.dNdS))) 

Prietella_dNdS_bs.ODB10 <- Prietella_dNdS_bs.ODB10 %>% mutate(Gp = Tp/3)


Prietella_dNdS_bs.ODB12 <- 
  Prietella_dNdS_bs.ODB12 %>%
  mutate(Tp = 
           DivergenceTime_P_A.ODB12 - ((DivergenceTime_P_A.ODB12 * (Cave_dNdS - neutral.dNdS)) / 
                                         (Surface_dNdS - neutral.dNdS)))

Prietella_dNdS_bs.ODB12 <- Prietella_dNdS_bs.ODB12 %>% mutate(Gp = Tp/3)


## CUL9

CUL9_dNdS_bs <- 
  read.table("~/SiluCave/Meredith_method/CUL9/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(CUL9_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

CUL9_dNdS_bs.ODB10 <- CUL9_dNdS_bs
CUL9_dNdS_bs.ODB12 <- CUL9_dNdS_bs


DivergenceTime_CUL9_RHP1.ODB10 <- 8.188 * 1000000
DivergenceTime_CUL9_RHP1.ODB12 <- 7.21 * 1000000

CUL9_dNdS_bs.ODB10 <- 
  CUL9_dNdS_bs.ODB10 %>%
  mutate(Tp = 
           DivergenceTime_CUL9_RHP1.ODB10 - ((DivergenceTime_CUL9_RHP1.ODB10 * (Cave_dNdS - neutral.dNdS)) / 
                                               (Surface_dNdS - neutral.dNdS))) 
CUL9_dNdS_bs.ODB10 <-  CUL9_dNdS_bs.ODB10 %>% mutate(Gp = Tp/3)


CUL9_dNdS_bs.ODB12 <- 
  CUL9_dNdS_bs.ODB12 %>%
  mutate(Tp = 
           DivergenceTime_CUL9_RHP1.ODB12 - ((DivergenceTime_CUL9_RHP1.ODB12 * (Cave_dNdS - neutral.dNdS)) / 
                                               (Surface_dNdS - neutral.dNdS))) 
CUL9_dNdS_bs.ODB12 <-  CUL9_dNdS_bs.ODB12 %>% mutate(Gp = Tp/3)

## CSV83

CSV83_dNdS_bs <- 
  read.table("~/SiluCave/Meredith_method/CSV83/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(CSV83_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

CSV83_dNdS_bs.ODB10 <- CSV83_dNdS_bs
CSV83_dNdS_bs.ODB12 <- CSV83_dNdS_bs


DivergenceTime_CSVandCUL4_RUI2.ODB10 <- 14.6389 * 1000000
DivergenceTime_CSVandCUL4_RUI2.ODB12 <- 12.06 * 1000000

CSV83_dNdS_bs.ODB10 <- 
  CSV83_dNdS_bs.ODB10 %>%
  mutate(Tp = 
           DivergenceTime_CSVandCUL4_RUI2.ODB10 - ((DivergenceTime_CSVandCUL4_RUI2.ODB10 * (Cave_dNdS - neutral.dNdS)) / 
                                                     (Surface_dNdS - neutral.dNdS))) 
CSV83_dNdS_bs.ODB10 <- CSV83_dNdS_bs.ODB10 %>% mutate(Gp = Tp/3)



CSV83_dNdS_bs.ODB12 <- 
  CSV83_dNdS_bs.ODB12 %>%
  mutate(Tp = 
           DivergenceTime_CSVandCUL4_RUI2.ODB12 - ((DivergenceTime_CSVandCUL4_RUI2.ODB12 * (Cave_dNdS - neutral.dNdS)) / 
                                                     (Surface_dNdS - neutral.dNdS))) 
CSV83_dNdS_bs.ODB12 <- CSV83_dNdS_bs.ODB12 %>% mutate(Gp = Tp/3)


## CUL4

CUL4_dNdS_bs <- 
  read.table("~/SiluCave/Meredith_method/CUL4/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(CUL4_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

CUL4_dNdS_bs.ODB10 <- CUL4_dNdS_bs
CUL4_dNdS_bs.ODB12 <- CUL4_dNdS_bs


DivergenceTime_CSVandCUL4_RUI2.ODB10 <- 14.6389 * 1000000
DivergenceTime_CSVandCUL4_RUI2.ODB12 <- 12.06 * 1000000

CUL4_dNdS_bs.ODB10 <- 
  CUL4_dNdS_bs.ODB10 %>%
  mutate(Tp = 
           DivergenceTime_CSVandCUL4_RUI2.ODB10 - ((DivergenceTime_CSVandCUL4_RUI2.ODB10 * (Cave_dNdS - neutral.dNdS)) / 
                                                     (Surface_dNdS - neutral.dNdS))) 
CUL4_dNdS_bs.ODB10 <- CUL4_dNdS_bs.ODB10 %>% mutate(Gp = Tp/3)


CUL4_dNdS_bs.ODB12 <- 
  CUL4_dNdS_bs.ODB12 %>%
  mutate(Tp = 
           DivergenceTime_CSVandCUL4_RUI2.ODB12 - ((DivergenceTime_CSVandCUL4_RUI2.ODB12 * (Cave_dNdS - neutral.dNdS)) / 
                                                     (Surface_dNdS - neutral.dNdS))) 
CUL4_dNdS_bs.ODB12 <- CUL4_dNdS_bs.ODB12 %>% mutate(Gp = Tp/3)

## CHM6

CHM6_dNdS_bs <- 
  read.table("~/SiluCave/Meredith_method/CHM6/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(CHM6_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

CHM6_dNdS_bs.ODB10 <- CHM6_dNdS_bs
CHM6_dNdS_bs.ODB12 <- CHM6_dNdS_bs

DivergenceTime_CHM6_RSS1.ODB10 <- 10.5187 * 1000000
DivergenceTime_CHM6_RSS1.ODB12 <- 8.482 * 1000000

CHM6_dNdS_bs.ODB10 <- 
  CHM6_dNdS_bs.ODB10 %>%
  mutate(Tp = 
           DivergenceTime_CHM6_RSS1.ODB10 - ((DivergenceTime_CHM6_RSS1.ODB10 * (Cave_dNdS - neutral.dNdS)) / 
                                               (Surface_dNdS - neutral.dNdS))) 
CHM6_dNdS_bs.ODB10 <- CHM6_dNdS_bs.ODB10 %>% mutate(Gp = Tp/3)


CHM6_dNdS_bs.ODB12 <- 
  CHM6_dNdS_bs.ODB12 %>%
  mutate(Tp = 
           DivergenceTime_CHM6_RSS1.ODB12 - ((DivergenceTime_CHM6_RSS1.ODB12 * (Cave_dNdS - neutral.dNdS)) / 
                                               (Surface_dNdS - neutral.dNdS))) 
CHM6_dNdS_bs.ODB12 <- CHM6_dNdS_bs.ODB12 %>% mutate(Gp = Tp/3)


## Trichomycterus_rosablanca

Trichomycterus_rosablanca_dNdS_bs <- 
  read.table("~/SiluCave/Meredith_method/Trichomycterus_rosablanca/dNdS_values_bootstraps.csv",
             sep=",",
             header=FALSE)
colnames(Trichomycterus_rosablanca_dNdS_bs) <- c("bootstrap_nb", "Surface_dNdS", "Cave_dNdS")

Trichomycterus_rosablanca_dNdS_bs.ODB10 <- Trichomycterus_rosablanca_dNdS_bs
Trichomycterus_rosablanca_dNdS_bs.ODB12 <- Trichomycterus_rosablanca_dNdS_bs

DivergenceTime_T_RSS1.ODB10 <- 15.1434 * 1000000
DivergenceTime_T_RSS1.ODB12 <- 12.985 * 1000000

Trichomycterus_rosablanca_dNdS_bs.ODB10 <- 
  Trichomycterus_rosablanca_dNdS_bs.ODB10 %>%
  mutate(Tp = 
           DivergenceTime_T_RSS1.ODB10 - ((DivergenceTime_T_RSS1.ODB10 * (Cave_dNdS - neutral.dNdS)) / 
                                            (Surface_dNdS - neutral.dNdS))) 
Trichomycterus_rosablanca_dNdS_bs.ODB10 <- Trichomycterus_rosablanca_dNdS_bs.ODB10 %>% mutate(Gp = Tp/3)



Trichomycterus_rosablanca_dNdS_bs.ODB12 <- 
  Trichomycterus_rosablanca_dNdS_bs.ODB12 %>%
  mutate(Tp = 
           DivergenceTime_T_RSS1.ODB12 - ((DivergenceTime_T_RSS1.ODB12 * (Cave_dNdS - neutral.dNdS)) / 
                                            (Surface_dNdS - neutral.dNdS))) 
Trichomycterus_rosablanca_dNdS_bs.ODB12 <- Trichomycterus_rosablanca_dNdS_bs.ODB12 %>% mutate(Gp = Tp/3)



## Summary 

Summary_dNdS_datation.ODB10 <- 
  rbind(
    Prietella_dNdS_bs.ODB10 %>% mutate(species = "Prietella_phreatophila"),
    CUL9_dNdS_bs.ODB10 %>% mutate(species = "CUL9"),
    CSV83_dNdS_bs.ODB10 %>% mutate(species = "CSV83"),
    CUL4_dNdS_bs.ODB10 %>% mutate(species = "CUL4"),
    CHM6_dNdS_bs.ODB10 %>% mutate(species = "CHM6"),
    Trichomycterus_rosablanca_dNdS_bs.ODB10 %>% mutate(species = "Trichomycterus_rosablanca")
  )

Summary_dNdS_datation.ODB12 <- 
  rbind(
    Prietella_dNdS_bs.ODB12 %>% mutate(species = "Prietella_phreatophila"),
    CUL9_dNdS_bs.ODB12 %>% mutate(species = "CUL9"),
    CSV83_dNdS_bs.ODB12 %>% mutate(species = "CSV83"),
    CUL4_dNdS_bs.ODB12 %>% mutate(species = "CUL4"),
    CHM6_dNdS_bs.ODB12 %>% mutate(species = "CHM6"),
    Trichomycterus_rosablanca_dNdS_bs.ODB12 %>% mutate(species = "Trichomycterus_rosablanca")
  )

Summary_dNdS_datation.ODB10 %>% filter(bootstrap_nb == "original_aln")
Summary_dNdS_datation.ODB12 %>% filter(bootstrap_nb == "original_aln")




Summary_dNdS_datation.ODB10 %>%
  ggplot(., aes(x=Tp, color=species)) +
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
        plot.subtitle=element_text(size=16)) +
  theme(legend.position="none") +
  xlim(0, 5000000) +
  geom_density(data = Summary_dNdS_datation.ODB12, linetype = "dashed") +
  scale_color_manual(values = c(
    "Prietella_phreatophila" = "#000000",
    "Trichomycterus_rosablanca" = "#E69F00",
    "CHM6" = "#D55E00",
    "CSV83" = "#CC79A7",
    "CUL4" = "#009E73",
    "CUL9" = "#56B4E9"
  )) +
  xlab("Years ago")




