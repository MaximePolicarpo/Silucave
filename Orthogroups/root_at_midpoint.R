library("ape")
library("phytools")

args <- commandArgs(trailingOnly = TRUE)

current_tree <- read.tree(args[1])

#Root tree at midpoint 
current_tree_rooted <- midpoint.root(current_tree)

write.tree(current_tree_rooted, args[2])
