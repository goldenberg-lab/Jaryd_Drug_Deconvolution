library(tidyverse)
library(cowplot)

# Compare the mTOR synergies found in Ceryl's experiments with the synergies
# found in the linear model.

# Make the folder to hold figures
fig_dir <- './post_fit_analysis/Synergies_Validation_attempt/CerylmTOR/figures/'

if(!dir.exists(fig_dir)) {
  dir.create(fig_dir)
  if(!dir.exists(fig_dir)){
    stop('Failed to create directory for figures. Are you running from directory
         ./CellSummerProject?')
  }
}

mTORSyns <- read_csv('../../CellSummerProjectData/CerylMtor/AllsiRNAData.csv')

tidied <- read_csv('../../CellSummerProjectData/CerylMtor/AllsiRNAData.csv') %>%
  separate(Rap, into = c('Rap_nM'), sep = '(?=[a-zA-Z])', extra = 'drop') %>%
  mutate(Rap_nM = as.numeric(Rap_nM))

# Outcomes histograms.
facet_hists <- function(dat, save_file_name, ...){
  vars <- enquos(...)
  dat %>% select(!!! vars) %>%
    gather() %>% ggplot(., aes(x = value)) + 
    geom_histogram(bins = 30) + facet_wrap('~key', scales = 'free')
  ggsave(save_file_name)
}

str_grp <- c('Area', 'Num', 'Size', 'Props')
lapply(str_grp, 
       function(x){facet_hists(tidied, 
                               paste0(fig_dir, x, '_measures_Hist.pdf'), 
                               contains(x))})

# Very few genes in the G1S phase.
mTORSyns %>% mutate(num_G1S = Props_G1S*CellNum_raw) %>% 
  ggplot(., aes(x = num_G1S)) + geom_histogram(binwidth = 1)

# check for batch effects, left commented out since it takes awhile to run due
# to labelling.
# library(ggrepel)
# tmp <- tidied %>% mutate(label = if_else(CellNum_clean < 1000 & CellNum_raw > 4000, Genes, ''))
# tmp %>% ggplot(., aes(x = CellNum_raw, y = CellNum_clean, label = label)) +
#   geom_point() + geom_text_repel()
# ggsave(paste0(fig_dir, 'Extreme_Cleaning_Pal.pdf'))

# load in synergy coefficients
My_Syn <- read_csv(paste0('./Coefs/Penalty/CellNum/CellNumMinDrgs:4_alpha:0_',
                          'pairs:TRUE_Inter:FALSE_avgreps:TRUE_concs:10,20',
                          'lambda1se_rmvdups2019-03-25_Interp.csv')) %>%
  filter(group_id != group_id2 & isinterp1020 == T) %>%
  select(-genekey, -isinterp1020) %>%
  rename(Syn_coef = x)

grp1020 <- read_csv('./grouping/latestgroups10,20CellNum.csv') %>%
  mutate(group_sym = if_else(is.na(group_sym), ortho_sym, group_sym))

Split_My_Syn <- My_Syn %>%
  left_join(., grp1020 %>% select(-gene_id, -gene_sym) %>% distinct(), 
            by = c('group_sym', 'group_id')) %>%
  mutate(gene_sym = if_else(is.na(ortho_sym), group_sym, ortho_sym)) %>%
  select(-group_num, -group_id, -group_sym, -num_drugs, -ortho_id, -ortho_sym) %>%
  distinct() %>%
  rename(LHS = gene_sym) %>%
  left_join(., grp1020 %>% select(-gene_id, -gene_sym) %>% distinct(), 
            by = c('group_sym2' = 'group_sym', 'group_id2' = 'group_id')) %>%
  mutate(gene_sym = if_else(is.na(ortho_sym), group_sym2, ortho_sym)) %>%
  select(-group_num, -group_id2, -group_sym2,  -num_drugs, -ortho_id, -ortho_sym) %>%
  distinct() %>%
  rename(RHS = gene_sym) %>%
  mutate(gene_1 = if_else(LHS < RHS, LHS, RHS),
         gene_2 = if_else(LHS < RHS, RHS, LHS)) %>%
  select(-LHS, -RHS)

myMtor <- Split_My_Syn %>% filter(gene_1 == 'MTOR' | gene_2 == 'MTOR') %>%
  mutate(gene = if_else(gene_1 == 'MTOR', gene_2, gene_1)) %>%
  select(Syn_coef, gene)

# first attempt, will calculate correlation of formula in notebook with the coefficients.
tmp <- left_join(tidied, myMtor, by = c('Genes' = 'gene')) %>% 
  filter((!is.na(Syn_coef) | Genes == 'DMSO') & (Rap_nM == 30 | Rap_nM == 0)) %>%
  mutate(Syn_coef = if_else(is.na(Syn_coef), 0, Syn_coef)) %>%
  select(Genes, Rap_nM, CellNum_clean) %>%
  group_by(Genes, Rap_nM) %>%
  summarize_all(mean)
  
V_1 <- tmp %>% filter(Genes == 'DMSO' & Rap_nM == 0) %>% .$CellNum_clean
V_2 <- tmp %>% filter(Genes == 'DMSO' & Rap_nM == 30) %>% .$CellNum_clean
tmp <- tmp %>% filter(!Genes == 'DMSO') %>%
  group_by(Genes) %>%
  spread(Rap_nM, CellNum_clean, sep = '_') %>%
  mutate(act_syn = -(V_2) - (Rap_nM_0) + (Rap_nM_30)) %>%
  left_join(., myMtor, by = c('Genes' = 'gene'))

library(ggrepel)
ggplot(tmp, aes(x = Syn_coef, y = act_syn, label = Genes)) + geom_point() +
  ggtitle("Coefficient vs actual values from Ceryls experiments, rho=0.07, pval=0.8") +
  geom_text_repel()
ggsave(paste0(fig_dir, 'Syn_effect.pdf'))

cor.test(tmp$Syn_coef, tmp$act_syn, method = c("spearman"), alternative = 'two.sided')

#####
# Remove CCND1 which could be conceived as an outlier, see if this reveals a good correlation.
rmv_CCND1 <- tmp %>% filter(Genes != 'CCND1')
cor.test(rmv_CCND1$Syn_coef, rmv_CCND1$act_syn, method = c('spearman'), alternative = 'two.sided')

ggplot(rmv_CCND1, aes(x = Syn_coef, y = act_syn, label = Genes)) + geom_point() +
  ggtitle("Rmvd CCND1, because it is an outlier, rho=0.3, pval=0.3") +
  geom_text_repel()
ggsave(paste0(fig_dir, 'Syn_effect_rmv_CCND1.pdf'))

#####
# Not a significant correlation nor a very good looking plot, continue once I hear
# back from Ceryl with answers to questions about what those drugs mean, perhaps
# they will make the plot better. Answer from Ceryl, the drugs are probably multi
# targeting and so she suggests to ignore them I have no reason as of now to use them.
# A: Heard back from Ceryl, the drugs targets are unknown, I will not consider.

# Check if there is even a correlation with the main gene effect coefs and the 
# observations in Ceryls data.

Ind_coefs <- read_csv(paste0('./Coefs/Penalty/CellNum/CellNumMinDrgs:4_alpha:0_',
                         'pairs:TRUE_Inter:FALSE_avgreps:TRUE_concs:10,20',
                         'lambda1se_rmvdups2019-03-25_Interp.csv')) %>%
  filter(group_id == group_id2 & isinterp1020 == T) %>%
  left_join(., grp1020 %>% select(-gene_id, -gene_sym, -group_id) %>% distinct(), 
            by = c('group_sym')) %>%
  mutate(gene_sym = if_else(is.na(ortho_sym), group_sym, ortho_sym)) %>%
  select(-group_num, -group_sym, -num_drugs, -ortho_id, -ortho_sym, -isinterp1020) %>%
  distinct()
  
Ceryl_Ind_eff <- tidied %>% filter(Rap_nM == 0)

test <- left_join(Ceryl_Ind_eff, Ind_coefs, by = c("Genes" = 'gene_sym')) %>%
  select(Genes, CellNum_raw, CellNum_clean, coef = x, group_sym2)

cor.test(test$CellNum_clean, test$coef, method = 'spearman')
ggplot(test, aes(x = coef, y = CellNum_clean)) + geom_point() +
  ggtitle("individual effect coefficients vs individual siRNA effect, rho=0.2, pval = 0.001")
ggsave(paste0(fig_dir, 'Ind_effect.pdf'))

test2 <- test %>% group_by(coef) %>% summarize(mean(CellNum_clean))
ggplot(test2, aes(x = coef, y = `mean(CellNum_clean)`)) + geom_point() +
  ggtitle('Mean CellNum per gene group, rho = 0.5, pval = 0.005')
cor.test(test2$coef, test2$`mean(CellNum_clean)`, method = 'spearman')
ggsave(paste0(fig_dir, "Ind_effect_group_means.pdf"))

























