library(ape)
library(phytools)

args = commandArgs(trailingOnly=TRUE)


mytree <- read.tree(args[1])
species_tree <- read.tree("AMAS_concatenated_alignment_BUSCO.fa.timetree.nwk")
pruned_species_tree <- keep.tip(species_tree, mytree$tip.label)
tree_rooted_labeled <- makeNodeLabel(pruned_species_tree, method="number", prefix="Node")


write.tree(tree_rooted_labeled, file = args[2])