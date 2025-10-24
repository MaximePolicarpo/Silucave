#!/bin/bash


module purge ; module load R/4.4.1-foss-2023b

current_tree=$1 
current_tree_rooted=$2


Rscript root_at_midpoint.R $current_tree $current_tree_rooted