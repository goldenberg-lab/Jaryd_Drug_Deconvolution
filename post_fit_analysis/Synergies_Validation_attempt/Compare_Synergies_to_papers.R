# check correlation of synergy scores found in paper: 
# "Construction of synergy networks from gene expression data related to disease"

library(tidyverse)
library(parallel)

Ind_Syn <- read_csv('./Papers/Construction of synergy networks from gene expression data related to disease found synergies.csv')

My_Syn <- read_csv(paste0('./Coefs/Penalty/CellNum/CellNumMinDrgs:4_alpha:0_',
                          'pairs:TRUE_Inter:FALSE_avgreps:TRUE_concs:10,20',
                          'lambda1se_rmvdups2019-03-25.csv')) %>%
  filter(group_id != group_id2) %>%
  select(-genekey) %>%
  rename(Syn_coef = x)

grp1020 <- read_csv('./grouping/latestgroups10,20CellNum.csv') %>%
  mutate(group_sym = if_else(is.na(group_sym), ortho_sym, group_sym))

#####
# Q: How many individual genes are in our genes and this Independent study.
Ind_genes <- unique(c(Ind_Syn$`Gene 1`, Ind_Syn$`Gene 2`)) %>%
  str_split(., '///') %>% Reduce(c, .) %>%
  trimws(.) %>% unique(.)
My_genes <- unique(c(My_Syn$group_sym, My_Syn$group_sym2)) %>%
  enframe(name = NULL) %>%
  left_join(., grp1020, by = c('value' = 'group_sym')) 

My_genes <- c(My_genes$value, My_genes$ortho_sym) %>%
  unique(.)

intersect(My_genes, Ind_genes) %>% 
  # a little silly to assign it this way since it is very unorthodox, but I 
  # thought it would be interesting to do it anyway.
  assign("g_in_both", ., envir = globalenv()) %>% 
  length(.)
# [1] 118
# A: When also considering all group members for genes that were included in at 
#   least a synergy pair, there are 118 genes overlapping between our two sets.
#####
# Join on the pairs that are in common

# Helper function to split the synonyms of a "gene" listed in a string apart.
Split_Synonyms <- function(dat, column, sep, id_col = 'Synonym_id', name_col = 'Synonym', id_suffix = 'id_'){
  col = enquo(column)
  max_synonyms <- str_count(select(dat, !!col)[[1]], sep) %>% max(.)
  dat %>% separate(!!col, paste0(id_suffix, 1:(max_synonyms + 1)), sep, 
                   fill = 'right') %>%
    gather(!!id_col, !!name_col, starts_with(!!id_suffix)) %>%
    mutate_at(vars(!!name_col), trimws) %>%
    filter(!is.na(!!sym(name_col)))
}

Split_Ind_Syn <- Ind_Syn %>%
  mutate(Synergy_id = 1:n()) %>%
  Split_Synonyms(., `Gene 1`, '///', id_col = 'LHSyn_id', name_col = 'LHSyn') %>%
  Split_Synonyms(., `Gene 2`, '///', id_col = 'RHSyn_id', name_col = 'RHSyn') %>%
  mutate(Gene_1 = if_else(LHSyn < RHSyn, LHSyn, RHSyn),
         Gene_2 = if_else(LHSyn < RHSyn, RHSyn, LHSyn)) %>%
  select(-LHSyn_id, -LHSyn, -RHSyn_id, -RHSyn, -Synergy_id)

# This table will be useful for later analysis also, since it has all synergies
# split into separate rows for each group member and sorted each row where the
# first gene will be a lower value for the character string than the second gene.
Split_My_Syn <- My_Syn %>%
  left_join(., grp1020, by = c('group_sym', 'group_id')) %>%
  mutate(gene_sym = if_else(is.na(ortho_sym), group_sym, ortho_sym)) %>%
  select(-gene_id, -group_num, -group_id, -group_sym, -num_drugs, -ortho_id, -ortho_sym) %>%
  rename(LHS = gene_sym) %>%
  left_join(., grp1020, by = c('group_sym2' = 'group_sym', 'group_id2' = 'group_id')) %>%
  mutate(gene_sym = if_else(is.na(ortho_sym), group_sym2, ortho_sym)) %>%
  select(-gene_id, -group_num, -group_id2, -group_sym2,  -num_drugs, -ortho_id, -ortho_sym) %>%
  rename(RHS = gene_sym) %>%
  mutate(Gene_1 = if_else(LHS < RHS, LHS, RHS),
         Gene_2 = if_else(LHS < RHS, RHS, LHS)) %>%
  select(-LHS, -RHS)

test <- left_join(Split_Ind_Syn, Split_My_Syn) %>% filter(!is.na(Syn_coef))
ggplot(test, aes(x = Syn_coef, y = `Synergy score`)) + geom_point()

# Looks like our method did not find similar things to this paper.

#####
# Look at enrichment of known Genetic interactions from GeneMania with our 
# ordered synergies.

# First remove objects from comparing our Synergies to the paper which are not 
#needed for this next part.
rm(list = setdiff(ls(), 'Split_My_Syn'))

Id2Sym <- read_tsv('../../CellSummerProjectData/NCBI/Homo_sapiens.gene_info') %>%
  select(GeneID, Symbol)

Ensb2Gid <- read_tsv('../../CellSummerProjectData/NCBI/gene2ensembl.tsv') %>%
  filter(`#tax_id` == 9606) %>%
  left_join(., Id2Sym) %>%
  select(GeneID, Ensembl_gene_identifier, Symbol) %>%
  rename(gene_id = GeneID, gene_sym = Symbol, Ensembl = Ensembl_gene_identifier) %>%
  distinct()

GM <- read_tsv('../../CellSummerProjectData/GeneMania/Genetic_interactions.tsv')

# Q: how many ensembl ids have more than one gene_sym in the NCBI data
GM %>% left_join(., Ensb2Gid, by = c('Gene_A' = 'Ensembl')) %>%
  select(Gene_A, gene_sym) %>%
  distinct() %>%
  group_by(Gene_A) %>%
  add_count(Gene_A) %>% filter(n != 1) -> duplicate_Syms
# A: about 81 of them, see tibble for what they are.

Synergies_Tbl <- GM %>% 
  left_join(., Ensb2Gid, by = c('Gene_A' = 'Ensembl')) %>%
  left_join(., Ensb2Gid, by = c('Gene_B' = 'Ensembl'), suffix = c('_A', '_B')) %>%
  select(gene_sym_A, gene_sym_B) %>%
  # Again, order the genes so that the lowest name is the first one. I don't 
  # like the naming scheme of Gene_1 and Gene_2, will go back and change one day.
  mutate(Gene_1 = if_else(gene_sym_A < gene_sym_B, gene_sym_A, gene_sym_B),
         Gene_2 = if_else(gene_sym_A < gene_sym_B, gene_sym_B, gene_sym_A),
         known = T) %>%
  select(-gene_sym_A, -gene_sym_B) %>%
  full_join(., Split_My_Syn, by = c('Gene_1', 'Gene_2')) %>%
  filter(!is.na(Syn_coef)) %>%
  mutate(known = if_else(is.na(known), F, known))

# Prior notes on how I sped up naive implementation of Thresh_Fish
#####
# # There might be potential for speedups since the fisher test is only a couple 
# # milliseconds to run, but I don't think there is any way faster than tabulate 
# # to get all the counts out of the data frame.
# # A sample of tests done to create Faster_TF
# library(microbenchmark)
# tmpdat <- tibble(rnorm(100000))
# tmpdat <- tmpdat %>% mutate(in_set_col = as.logical(runif(nrow(.))))
# microbenchmark::microbenchmark(tmpdat %>% mutate(pos = c(rep(2, 7), rep(4, dim(tmpdat)[[1]] - 7)),
#                                                  pos = pos - in_set_col) %>%
#                                  .$pos %>%
#                                  tabulate(., nbins = 4), 
#                                tmpdat %>% mutate(rnk = seq.int(nrow(.)),
#                                                  pos = if_else(rnk <= 100, 
#                                                                if_else(in_set_col == T, 1, 2),
#                                                                if_else(in_set_col == T, 3, 4)),
#                                                  pos = factor(pos, levels = c(1,2,3,4))) %>%
#                                  .$pos %>%
#                                  table(.))
# 
# # This didn't make any difference to the speed.
# xint <- 1L:100000L
# yint <- rep(1L, 100000)
# ylog <- rep(T, 100000)
# microbenchmark(xint - yint, xint - ylog)
# 
# # check that for a set number of thresholds they produce the same pvalues
# thresholds <- 1:100
# fishers1 <- lapply(thresholds, function(x){print(x);Faster_TF(Synergies_Tbl, x, Syn_coef, known)})
# fishers2 <- lapply(thresholds, function(x){print(x);Thresholded_fishers(Synergies_Tbl, x, Syn_coef, known)})
# 
# pvals1 <- sapply(fishers1, function(x){x$p.value})
# pvals2 <- sapply(fishers2, function(x){x$p.value})
# # the two vectors are identical looks like faster is as correct as the slow 
# # implementation
#####

# Fast implementation of Thresholded Fishers, assumes in_set_vec is presorted.
# 1 if <T & known, 2 if <T & !known, 3 if > T & known, 4 if > T & !known
Thresh_Fish <- function(threshold, in_set_vec){
  {c(rep(2L, threshold), rep(4L, length(in_set_vec) - threshold)) - in_set_vec} %>%
    tabulate(., nbins = 4) %>%
    matrix(., ncol = 2) %>%
    fisher.test(.)
}

# Calculates the pvalue for each permutation given.
Perm_Fish <- function(threshold, in_set_vecs){
  pvals <- numeric(length = length(in_set_vecs))
  for (i in 1:length(in_set_vecs)){
    pvals[[i]] <- Thresh_Fish(threshold, in_set_vecs[[i]])$p.value
  }
  pvals
}


thresholds <- 1:dim(Synergies_Tbl)[[1]]
nperm <- 10
permutations <- vector("list", length = nperm)
permutations[[1]] <- arrange(Synergies_Tbl, desc(abs(Syn_coef))) %>% .$known
for (i in 2:nperm){
  permutations[[i]] <- sample(permutations[[1]], length(permutations[[1]]))
}
names(permutations) <- paste0('Starting Permutation: ', 1:nperm)
fisher_pvals <- mcmapply(permutations, names(permutations), FUN = 
                          function(p, id){print(id); 
                            sapply(thresholds, 
                                   function(t){if(t%%1000 == 0){print(t)}; 
                                     Thresh_Fish(t, p)$p.value})},
                    mc.cores = 7)

# get the pvalue given the list of fisher.test results, across thresholds, 
# across permutations.
best_pvals <- sapply(fisher_pvals, min)

pval <- (1 + sum(best_pvals < best_pvals[[1]]))/(nperm)

save.image('./post fit analysis/Synergies_Validation_attempt/Fishers_results_test.RData')








