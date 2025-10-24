library(ape)
library(phytools)

args = commandArgs(trailingOnly=TRUE)


mytree <- read.tree(args[1])
tree_rooted <- midpoint.root(mytree)
tree_rooted$node.label <- NULL
tree_rooted_labeled <- makeNodeLabel(tree_rooted, method="number", prefix="Node")


write.tree(tree_rooted_labeled, file = args[2])