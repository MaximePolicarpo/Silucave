library(ape)


args = commandArgs(trailingOnly=TRUE)

mytree <- read.tree(args[1])
mytree$node.label <- rep(" {reference}",  mytree$Nnode)
write.tree(mytree, args[2])