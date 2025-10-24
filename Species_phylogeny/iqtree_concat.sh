#!/bin/bash


#SBATCH --job-name=IQTREE   # Job name


iqtree -s AMAS_concatenated_alignment_BUSCO.fa --seqtype AA -m LG+F+G4 -nt 10 -bb 1000
