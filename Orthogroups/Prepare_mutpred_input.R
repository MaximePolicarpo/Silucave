library(treeio)
library(Biostrings)
library(dplyr)
library(ggtree)
library(ape)
library(phytools)

args = commandArgs(trailingOnly=TRUE)


my_rst <- read.paml_rst(args[1])

all_tips_v <- my_rst@phylo$tip.label

### Make a corrrespondance table between nodes numbers and labels

species_test <- "misc"
misc_table <- as.data.frame(species_test) %>% mutate(value = 1)
colnames(misc_table) <- c("node","value")
misc_table$node <- as.numeric(misc_table$node)
node_label_corresp <- 
  left_join(my_rst@phylo, misc_table, by = 'node') %>%
  dplyr::select(node, label)
node_nb <- 
  node_label_corresp %>%
  filter(., grepl("Node", label)) %>%
  pull(node)
species_names <- 
  node_label_corresp %>%
  filter(., !grepl("Node", label)) %>%
  pull(label)
species_and_nodes <- c(species_names, node_nb)
node_label_corresp <- node_label_corresp %>% mutate(tree_names = species_and_nodes)
node_label_corresp <- as.data.frame(node_label_corresp)

## Now for each tips, extract the ancestral sequence and amino acid substitutions

AAML_summary_df <- as.data.frame(NULL)
for(curr_tip in all_tips_v){
  curr_node <- node_label_corresp %>% filter(label == curr_tip) %>% pull(node)
  curr_parent <- getParent(my_rst@phylo, curr_node)
  curr_AA_subst <- my_rst@data %>% filter(node == curr_node) %>% pull(subs)
   
  
  
  labels <- names(my_rst@anc_seq)
  parent_sequence_index <- which(labels == curr_parent)

  
  curr_df <- as.data.frame(cbind(curr_tip, curr_AA_subst, curr_parent))
  colnames(curr_df) <- c("tip", "AA_subst", "parent")
  AAML_summary_df <- rbind(AAML_summary_df, curr_df)
}

write.table(AAML_summary_df,
            args[2],
            sep=",",
            quote=FALSE,
            col.names=FALSE,
            row.names=FALSE)