library("ape")
library("phytools")

args = commandArgs(trailingOnly=TRUE)

mytree <- read.tree(args[1])
mytree$node.label <- NULL
write.tree(mytree, args[2])