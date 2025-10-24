#!/bin/bash


#SBATCH --job-name=ASTRAL   # Job name


java -jar Astral/astral.5.7.8.jar -i All_genes_tree.concatnwk -o Astral_nolength.nwk