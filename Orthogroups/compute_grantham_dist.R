library(dplyr)
library(grantham)
library(seqinr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)


aa_subst_df <- 
  as.data.frame(
    fread(args[1],
          header=FALSE,
          sep=" ")
  )

colnames(aa_subst_df) <- c("HOG","seq_name", "aa_subst")

aa_subst_df <- 
  aa_subst_df %>%
  rowwise() %>%
  mutate(aa_ancestral = substr(aa_subst, 1, 1)) %>%
  mutate(aa_new = substr(aa_subst, nchar(aa_subst)-1+1, nchar(aa_subst))) %>%
  mutate(aa_ancestral_t = aaa(aa_ancestral)) %>%
  mutate(aa_new_t = aaa(aa_new)) %>%
  mutate(grantham_dist = grantham_distance(aa_ancestral_t, aa_new_t) %>% pull(d))

write.table(aa_subst_df, args[2], sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)
