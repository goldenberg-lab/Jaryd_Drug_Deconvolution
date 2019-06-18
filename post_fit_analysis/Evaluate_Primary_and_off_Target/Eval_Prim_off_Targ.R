# for a given model evaluate how much of the modelled effects come from primary 
# and how much from off targets.
library(tidyverse)

# read in coefs. Reading in the best CellNum Models coefficients, which is the 
# model we validated.
coefs <- read_csv(paste0('./Coefs/Penalty/CellNum/CellNumMinDrgs:4_alpha:0_',
                         'pairs:TRUE_Inter:FALSE_avgreps:TRUE_concs:10,20',
                         'lambda1se_rmvdups2018-12-19.csv'))

# Read in created Primary Target table
PrimaryTargs <- read_csv('./post_fit_analysis/Evaluate_Primary_and_off_Target/PrimaryTargetTbl_10,20.csv')

# Get number of Primary Targets per Drug, i.e. how many genes are tied for 
# lowest IC50 value. 
NumPrimaryTargs <- PrimaryTargs %>%
  group_by(inchi_key, numGenesTargeted) %>%
  summarise(numPrimary = sum(primaryTarg, na.rm = T))

TiedPrimaryonlyDrugs <- NumPrimaryTargs %>% filter(numPrimary != 1)
NumTied <- TiedPrimaryonlyDrugs %>% .[[1]] %>% length()
  
TiedPrimary <- PrimaryTargs %>% filter(inchi_key %in% TiedPrimaryonlyDrugs$inchi_key) %>% 
  arrange(inchi_key)

# write table of drug with Primary target ties
write_csv(TiedPrimary, './post_fit_analysis/Evaluate_Primary_and_off_Target/Tied_Primary_Targets.csv')

source('./import_data.R')
lst <- import_grpd(rmv_cell_count_2 = F, average_replicate = T, 
                   add_pair_features = T, two_conc = c(10,20))
JoinedTable <- lst[[1]]
# rmv duplicate drugs, here until I implement this part in import_grpd.
JoinedTable <- JoinedTable %>% group_by(inchi_key) %>% 
  filter(Cpds.conc.uM == max(as.integer(as.character(Cpds.conc.uM)))) %>% 
  ungroup()

Phenotype <- JoinedTable %>% select(inchi_key, CellNum) %>% distinct()

# need fitted coefficients from the model to compare
coefs <- read.csv(paste0('./Coefs/Penalty/CellNum/CellNumMinDrgs:4_alpha:0',
                         '_pairs:TRUE_Inter:FALSE_avgreps:TRUE_concs:10,20',
                         'lambda1se_rmvdups2019-03-25.csv'))

coefsInd <- coefs %>% filter(group_id == group_id2) %>% select(x, group_id, group_sym)

# TODO: Deal with Pairs, when thinking about proportion perturbance for a gene, 
# is it the sum with the interaction terms is it only it's individual effect? 
# For, now just plot the main effect coefficients.
PerturbanceTbl <- JoinedTable %>% select(inchi_key, group_id, group_sym, 
                                         group_id2, group_sym2, CellNum) %>%
  distinct() %>% #duplicate rows created due to projection, losing "Homol" cols are the cause.
  left_join(., PrimaryTargs, by = c('inchi_key', 'group_id', 'group_sym')) %>% 
  left_join(., select(PrimaryTargs, -micromolar_activity), 
            by = c('inchi_key', 'group_id2' = 'group_id', 
                   'group_sym2' = 'group_sym', 'numGenesTargeted')) %>%
  left_join(., coefs) %>%
  # If either gene in the gene pair is the primary Target label the coefficient as a primary Target.
  mutate(primaryTarg = if_else(is.na(primaryTarg.x & primaryTarg.y), F, primaryTarg.x | primaryTarg.y)) %>%
  select(-primaryTarg.x, -primaryTarg.y) %>% 
  # get coefficient sum, and Residual
  group_by(inchi_key) %>% mutate(Predicted = sum(x)) %>%
  group_by(inchi_key) %>% mutate(Residual = CellNum - Predicted) %>%
  # only sums per drug so that it is the effect the drug is explained by primary
  # and off targets.
  group_by(inchi_key) %>% mutate(SumAbsCoef = sum(abs(x))) %>%
  mutate(AbsSumCoef = abs(sum(x))) %>%
  mutate(ProportionPrediction = abs(x)/((SumAbsCoef)))

# I also tried defining the PredWell to be drugs where 
# 0.5*CellNum >= Predicted >= 2*CellNum, which might be better if wanting to 
# look at drugs that had an actual value close to zero, but it does eliminate 
# some of the drugs where it predicted say -40 but it was actually -100, so it 
# was obviously picking up on some signal but not enough to accurately describe 
# the total drug effect.
RankPredicted <- PerturbanceTbl %>% 
  filter(Predicted >= 2*CellNum & Predicted <= 0.5*CellNum) %>% 
  arrange(desc(abs(Predicted))) %>%
  ungroup() %>% select(inchi_key) %>% distinct() %>% mutate(PredWell = row_number())

PerturbanceTbl <- PerturbanceTbl %>% left_join(., RankPredicted, by = 'inchi_key')

# write to file for use in other scripts that won't have to source this one.
write_csv(PerturbanceTbl, paste0('./post_fit_analysis/',
                                 'Evaluate_Primary_and_off_Target/',
                                 'FeaturePerturbanceTable_CellNum_',
                                 format(Sys.Date(), '%b-%d-%Y' ), '.csv'))

# Plots suggested by Shixuan, latest changes, see draft storage for old versions
library(ggrepel)
plots <- list() 

# How much any single gene is responsible for all predicted effects for a single
# compound, interesting to see if our model learnt a more complex representation
# of the effect for a single compound on the phenotype, or if it mostly 
# attributed all effect to a single gene.
i<- 1
plots[[i]] <- PerturbanceTbl %>% 
  select(inchi_key, Predicted, CellNum, PredWell) %>% distinct() %>%
  mutate(PredWellLab = if_else(PredWell < 21, as.character(PredWell), '')) %>%
  ggplot(., aes(x = CellNum, y = Predicted, colour = PredWell, label = PredWellLab)) + 
  geom_text_repel() +
  geom_point() + geom_abline(intercept = 0, slope = 1) 

i <- i + 1
plots[[i]] <- PerturbanceTbl %>% filter(group_id == group_id2) %>%
  group_by(inchi_key, primaryTarg) %>%
  mutate(ProportionPrediction == max(ProportionPrediction)) %>%
  ggplot(., aes(x = ProportionPrediction)) + 
  geom_histogram() + 
  facet_wrap(~primaryTarg, labeller = label_both) +
  ggtitle('ProportionPrediction = $|coef|/((sum |coefs|*CmpdTargets))$; Individual Effects only') +
  ylab('count drugs')

i <- i + 1
# Looks at the rankings of the best predicted drugs against the Proportion Perturbance.
plots[[i]] <- PerturbanceTbl %>% filter(group_id == group_id2 & !is.na(PredWell)) %>%
  # For this plot remove drugs that I couldn't identify their primary target, 
  # only plot the off target with the largest proportion perturbance per compound.
  group_by(inchi_key) %>% 
  mutate(hasPrimaryTarg = sum(primaryTarg) > 0) %>% 
  filter(hasPrimaryTarg) %>%
  group_by(inchi_key, primaryTarg) %>%
  filter(ProportionPrediction == max(ProportionPrediction)) %>%
  select(inchi_key, PredWell, ProportionPrediction, primaryTarg) %>%
  distinct() %>%
  ggplot(., aes(x = PredWell, y = ProportionPrediction, colour = primaryTarg)) + geom_point() #+
  #facet_grid(~ primaryTarg, labeller=label_both)
i <- i + 1

# just the largest predictions on the drugs that had the most signal.
plots[[i]] <- PerturbanceTbl %>% filter(group_id == group_id2 & !is.na(PredWell)) %>%
  # For this plot remove drugs that I couldn't identify their primary target, 
  # only plot the off target with the largest proportion perturbance per compound.
  group_by(inchi_key) %>% 
  mutate(hasPrimaryTarg = sum(primaryTarg) > 0) %>% 
  filter(hasPrimaryTarg) %>%
  group_by(inchi_key, primaryTarg) %>%
  filter(ProportionPrediction == max(ProportionPrediction) & PredWell <= 20) %>%
  select(inchi_key, PredWell, ProportionPrediction, primaryTarg) %>%
  distinct() %>%
  ggplot(., aes(x = PredWell, y = ProportionPrediction, colour = primaryTarg, label = substr(inchi_key, 1, 5))) + 
  geom_point() + coord_cartesian(xlim = c(0,20)) + geom_text_repel()
i <- i + 1

plots[[i]] <- PerturbanceTbl %>% ungroup() %>% 
  select(group_id, group_id2, x) %>% distinct() %>% arrange(desc(abs(x))) %>% 
  mutate(rankGeneEffect = row_number()) %>% rowwise() %>% 
  mutate(rankGeneEffect = min(200, rankGeneEffect)) %>% select(-x) %>% 
  left_join(PerturbanceTbl, ., by = c('group_id', 'group_id2')) %>%
  group_by(inchi_key) %>% filter(x == max(x)) %>% 
  ggplot(., aes(x = CellNum, y = Predicted, colour = rankGeneEffect)) + 
  geom_point() + geom_abline(intercept = 0, slope = 1)

pdf(paste0('./post fit analysis/Evaluate_Primary_and_off_Target/PropPert_PredPerf_',
    format(Sys.Date(), '%b-%d-%Y' ), '.pdf'))
plots
dev.off()

# Build table of the form: drug, PredPerf, PrimaryTargCoef, largestOffTargetCoef
LrgTargets <- PerturbanceTbl %>% filter(group_id == group_id2) %>% 
  group_by(inchi_key, primaryTarg) %>% add_tally() %>% 
  mutate(offCount = if_else(primaryTarg == F, n, 0L), 
         primCount = if_else(primaryTarg == T, n, 0L)) %>% 
  filter(abs(x) == max(abs(x))) %>%
  select(inchi_key, relativeDiff, group_id, group_sym, primaryTarg, offCount, primCount) %>%
  ungroup() %>%
  mutate(primaryTarg = if_else(primaryTarg, 'primaryTarg', 'lrgstOffTarg')) %>%
  # creating group column as amalgamation of sym and id, so that I can use it as
  # a variable to spread on later, spread uses vars_pull which only returns one 
  # column so I need the identifiers for the genes to be a single column.
  mutate(group = paste0(group_sym, '_', group_id)) %>% select(-group_id, -group_sym) %>%
  group_by(inchi_key) %>%
  mutate(offCount = max(offCount), primCount = max(primCount)) %>%
  spread(primaryTarg, group) %>%
  # need to separate sym and id again to join the actual coefficients into this 
  # table.
  separate(primaryTarg, c("primTargSym", "primTargId"), sep = '_') %>%
  separate(lrgstOffTarg, c("offTargSym", "offTargId"), sep = '_') %>%
  mutate(primTargId = as.integer(primTargId), offTargId = as.integer(offTargId)) %>%
  left_join(., coefsInd, by = c("offTargSym" = "group_sym", "offTargId" = "group_id")) %>%
  rename(offTargX = x) %>%
  left_join(., coefsInd, by = c("primTargSym" = "group_sym", "primTargId" = "group_id")) %>%
  rename(primTargX = x)

ggplot(LrgTargets %>% filter(!is.na(offTargX) & !is.na(primTargX) & relativeDiff < 0.2), 
       aes(x = primTargX, y = offTargX)) + geom_point()
ggsave('./post fit analysis/Evaluate_Primary_and_off_Target/PTvLrgstOT.pdf')

