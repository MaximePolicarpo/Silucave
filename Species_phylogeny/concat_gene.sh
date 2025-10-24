#!/bin/bash

#SBATCH --job-name=concat_gene  # Job name


i=$1

for species in `cat species_list.txt` ; do
	if grep -q "^$i$" Common_busco_genes/$species.singlecopy.list ; then
		sed "s/>.*/>$species/g" BUSCO_$species/run_actinopterygii_odb10/busco_sequences/single_copy_busco_sequences/$i >> Alignements/$i
	fi
done
