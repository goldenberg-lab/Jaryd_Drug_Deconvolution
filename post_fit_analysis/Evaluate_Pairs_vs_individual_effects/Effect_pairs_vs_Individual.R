library(tidyverse)

# read in coefs. Reading in the best CellNum Models coefficients, which is the 
# model we validated.
coefs <- read_csv(paste0('./Coefs/Penalty/CellNum/CellNumMinDrgs:4_',
                         'alpha:0_pairs:TRUE_Inter:FALSE_avgreps:TRUE_',
                         'concs:10,20lambda1se_rmvdups2018-12-19.csv'))

source('./import_data.R')
lst <- import_grpd(rmv_cell_count_2 = F, two_conc = c(10,20), 
                   average_replicates = T, add_pair_features = T)
JoinedTable <- lst[[1]]
# rmv duplicate drugs, here until I implement this part in import_grpd.
JoinedTable <- JoinedTable %>% group_by(inchi_key) %>% 
  filter(Cpds.conc.uM == max(as.integer(as.character(Cpds.conc.uM)))) %>% 
  ungroup()

# get total sum for coefficients of individual effects and the total sum for 
# coefficients of pair effects.
Totaleffects <- coefs %>% 
  mutate(pair_coef = group_id != group_id2) %>% 
  group_by(pair_coef) %>%
  summarise(totaleffect = sum(x))

# Many drugs have no pair effect since there were no pair coefficients fit to 
# those drugs since they probably weren't targeted by enough unique drugs to be 
# included.
Perdrug <- JoinedTable %>% 
  select(inchi_key, group_id, group_sym, group_id2, group_sym2) %>% 
  distinct() %>% 
  left_join(., coefs, by = c('group_id', 'group_sym', 'group_id2','group_sym2')) %>%
  mutate(pair = if_else(group_id != group_id2, T, F)) %>%
  group_by(inchi_key, pair) %>%
  summarise(effectSep = sum(x)) %>%
  group_by(inchi_key) %>%
  mutate(effectTotal = sum(effectSep)) %>%
  left_join(., JoinedTable %>% select(inchi_key, CellNum) %>% distinct())

# get list of drugs where there were no pairs fitter.
NoPairs <- Perdrug %>% 
  spread(pair, effectSep) %>%
  filter(is.na(`TRUE`)) %>% select(inchi_key) %>% distinct()

# count number of unique targets
test <- JoinedTable %>%
  select(inchi_key, group_id, group_id2) %>% distinct() %>% select(inchi_key) %>% 
  table() %>% as.data.frame() %>%
  rename('inchi_key' = '.', 'num_targets' = 'Freq') %>%
  left_join(Perdrug, .)

i <- 1
plots <- list() 

plots[[i]] <- test %>% ggplot(., aes(x = effectSep, y = num_targets, colour = pair)) + geom_point()

pdf('./post fit analysis/Evaluate_Pairs_vs_individual_effects/Pair_v_Ind.pdf')
plots
dev.off()








