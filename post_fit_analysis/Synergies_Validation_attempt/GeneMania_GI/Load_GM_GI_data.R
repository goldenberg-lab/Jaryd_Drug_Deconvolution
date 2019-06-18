library(tidyverse)

# Load the coefficients and GeneMania Genetic Interactions information, where
# all pairs are sorted and groups are split into separate rows.

My_Syn <- read_csv(paste0('./Coefs/Penalty/CellNum/CellNumMinDrgs:4_alpha:0_',
                          'pairs:TRUE_Inter:FALSE_avgreps:TRUE_concs:10,20',
                          'lambda1se_rmvdups2019-03-25_Interp.csv')) %>%
  filter(group_id != group_id2 & isinterp1020 == T) %>%
  select(-genekey, -isinterp1020) %>%
  rename(Syn_coef = x)

fill_NA_ids <- function(x){
  cur_num <- max(x, na.rm = T)
  sapply(x, function(num){if(is.na(num)){cur_num <<- cur_num + 1; cur_num}else{num}})
}

grp1020 <- read_csv('./grouping/latestgroups10,20CellNum.csv') %>%
  select(-gene_id, -gene_sym) %>%
  distinct() %>%
  mutate(group_sym = if_else(is.na(group_sym), ortho_sym, group_sym),
         group_num = fill_NA_ids(group_num),
         group_id = if_else(ortho_sym == group_sym, ortho_id, group_id))

Split_My_Syn <- My_Syn %>%
  left_join(., grp1020, 
            by = c('group_sym', 'group_id')) %>%
  select(-group_id, -group_sym, -num_drugs, -ortho_id, -group_num) %>%
  rename(LHS = ortho_sym) %>%
  left_join(., grp1020, 
            by = c('group_sym2' = 'group_sym', 'group_id2' = 'group_id')) %>%
  select(-group_id2, -group_sym2,  -num_drugs, -ortho_id, -group_num) %>%
  rename(RHS = ortho_sym) %>%
  mutate(Gene_1 = if_else(LHS < RHS, LHS, RHS),
         Gene_2 = if_else(LHS < RHS, RHS, LHS)) %>%
  select(-LHS, -RHS)

Id2Sym <- read_tsv('../../CellSummerProjectData/NCBI/Homo_sapiens.gene_info') %>%
  select(GeneID, Symbol)

Ensb2Gid <- read_tsv('../../CellSummerProjectData/NCBI/gene2ensembl.tsv') %>%
  filter(`#tax_id` == 9606) %>%
  left_join(., Id2Sym) %>%
  select(GeneID, Ensembl_gene_identifier, Symbol) %>%
  rename(gene_id = GeneID, gene_sym = Symbol, Ensembl = Ensembl_gene_identifier) %>%
  distinct()

GM <- read_tsv('../../CellSummerProjectData/GeneMania/Genetic_interactions.tsv') %>%
  mutate(Weight = as.numeric(Weight))
  
# Q: how many ensembl ids have more than one gene_sym in the NCBI data
GM %>% left_join(., Ensb2Gid, by = c('Gene_A' = 'Ensembl')) %>%
  select(Gene_A, gene_sym) %>%
  distinct() %>%
  group_by(Gene_A) %>%
  add_count(Gene_A) %>% filter(n != 1) -> duplicate_Syms
# A: 81 of them, see tibble for what they are.

Synergies_Tbl <- GM %>% 
  left_join(., Ensb2Gid, by = c('Gene_A' = 'Ensembl')) %>%
  left_join(., Ensb2Gid, by = c('Gene_B' = 'Ensembl'), suffix = c('_A', '_B')) %>%
  select(gene_sym_A, gene_sym_B, Weight) %>%
  # Again, order the genes so that the lowest name is the first one. I don't 
  # like the naming scheme of Gene_1 and Gene_2, will go back and change one day.
  mutate(Gene_1 = if_else(gene_sym_A < gene_sym_B, gene_sym_A, gene_sym_B),
         Gene_2 = if_else(gene_sym_A < gene_sym_B, gene_sym_B, gene_sym_A),
         known = T) %>%
  select(-gene_sym_A, -gene_sym_B) %>%
  full_join(., Split_My_Syn, by = c('Gene_1', 'Gene_2')) %>%
  filter(!is.na(Syn_coef)) %>%
  mutate(known = if_else(is.na(known), F, known))
