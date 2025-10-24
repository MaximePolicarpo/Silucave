library("ape")
library("phytools")
library("castor")
library("adephylo")
library("dplyr")


args <- commandArgs(trailingOnly = TRUE)


species_list <- scan("species_list.txt", what="character")

current_tree <- read.tree(args[1])

distances_vector <- distRoot(current_tree, tips = "all", method = "patristic")

distances_df <- data.frame(tiplabel = names(distances_vector), distance = as.vector(distances_vector))

tips_to_retain_vector <- c()
for(species in species_list){

	curr_sp_df <- distances_df %>% filter(grepl(species, tiplabel))
	nb_rows <- nrow(curr_sp_df)

	if(nb_rows == 1){
		tips_to_retain_vector <- c(tips_to_retain_vector, curr_sp_df %>% pull(tiplabel))
	} else if (nb_rows > 1) {
		tips_to_retain_vector <- c(tips_to_retain_vector, curr_sp_df %>% arrange(distance) %>% slice_head(n = 1) %>% pull(tiplabel))
	} 

}


write(tips_to_retain_vector, args[2])