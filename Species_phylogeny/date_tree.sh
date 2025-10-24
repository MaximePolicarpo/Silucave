#!/bin/bash


#SBATCH --job-name=ML_model_iqtree   # Job name


iqtree2 -s AMAS_concatenated_alignment_BUSCO.fa --date Actino_calibration.txt -te Astral_unrooted_nolength.nwk --date-options "-u 0.1" --date-tip 0 --date-ci 100 -o "Danio_rerio" -nt 30 -m LG+F+G4
