library(ape)

args = commandArgs(trailingOnly=TRUE)


mytree <- read.tree(args[1])
Badnodes <- which(as.numeric(mytree$node.label) <= 10) + length(mytree$tip.label)
Badnodes_indexes <- c()
for(node in Badnodes){ Badnodes_indexes <- c(Badnodes_indexes, which(mytree$edge[,2] == node)) }

mytree$edge.length[Badnodes_indexes] <- 0 
tree_multi <- di2multi(mytree) 
write.tree(tree_multi, file = args[2])