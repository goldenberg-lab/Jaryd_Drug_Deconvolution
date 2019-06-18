library(tidyverse)
# Assume all stats above have been calculated before stats below i.e. all 
# objects created to calculate above stats exist for lower calculations.

#####
# Get percentage of gene edges with JI>0.5, that have a GMConnection versus not a GMConnection.
# i.e. count_edges(Gmconnection and JI>0.5) / count_edges(JI>0.5)
Edges <- read_csv('./grouping/GroupingsByEdgesComplete.csv')

sum(Edges$GMConnection)/length(Edges$GMConnection)

#####
# Number of gene pairs connected due to JI==1 and have no known GM connection.
# Edges doesn't have duplicate edges, i.e. a -> b and b <- a it only lists one or the other.
length(Edges %>% filter(JaccardIdx == 1 & GMConnection == F) %>% .[[1]])

#####
# Number compound target pairs, and how many are quantitative vs qualitative
# (i.e. given IC50 value vs NaN)
originalPlate <- read_csv('../../CellSummerProjectData/OriginalData/PlateBookID2.txt')
Comp_targ_original <- originalPlate %>% 
  select(inchi_key, gene_id, micromolar_acitivity) %>% 
  distinct() %>%
  mutate(qualitative = is.na(micromolar_acitivity))

sum(Comp_targ_original$qualitative)/length(Comp_targ_original$qualitative)

#####
# Get counts of how many compound target pairs are "activating" vs "inhibiting"
activities_modality <- read_csv('./moa_activities/ActivitiesDirections.csv') %>%
  unite(moa_activity, starts_with('activity'))
orig_activities <- read_csv('../../CellSummerProjectData/OriginalData/Platebook_with_Activity.csv') %>%
  separate(activity_moa, paste0('activity_moa', seq(1,8))) %>%
  unite(moa_activity, starts_with('activity'))

activities_modality <- left_join(orig_activities, activities_modality) %>% 
  select(inchi_key, gene_id, moa_activity, `value reviewed by Shixuan`) %>% 
  distinct()

# get all compound target pairs where a modality is listed
modality_known <- activities_modality[grep('^unknown_NA', activities_modality$moa_activity, invert = T),]

# get percentage of all known modalities which contain inhibition
Inhib_present <- activities_modality[grep('inhibition', activities_modality$moa_activity),]
length(Inhib_present[[1]])/length(modality_known[[1]])

# get percentage of all known modalities which are only inhibition
Inhib_only <- activities_modality[grep('^inhibition_NA', activities_modality$moa_activity),]
length(Inhib_only[[1]])/length(modality_known[[1]])

# get percentage of all modalities which are unknown
unknown <- activities_modality[grep('^unknown_NA', activities_modality$moa_activity),]
length(unknown[[1]])/length(activities_modality[[1]])

# get percentage of known modalities which are negative i.e. same as inhibition
modality_known_summarised <- modality_known %>% count(`value reviewed by Shixuan`)
num_inhib <- (modality_known_summarised %>% filter(`value reviewed by Shixuan` == '-1'))$n
num_unknown <- sum((modality_known_summarised %>% filter(`value reviewed by Shixuan` %in% c('?', NA, '-1?', '1?')))$n)
total <- length(modality_known[[1]])

# number of those with defined modality which are negative
num_inhib/(total-num_unknown)
# number with a known defined modality
(total-num_unknown)/total

#####
# get number of compounds which had a single identifiable primary target on CellNum
PrimaryTargs <- read_csv('./post_fit_analysis/Evaluate_Primary_and_off_Target/PrimaryTargetTable_CellNum.csv')

# number of compounds with a single primary target
PrimaryTargs %>% group_by(inchi_key) %>% 
  filter(any(sum(primaryTarg == 1))) %>%
  .$inchi_key %>% unique() %>% length()

# number of compounds included in the model
PrimaryTargs %>% .$inchi_key %>% unique() %>% length()

#####
# get total number of predictors for final cell number model
source('import_data.R')

lst <- import_grpd(average_replicates = T, rmv_cell_count_2 = F,
                   create_group_files = F, two_conc = c(10,20),
                   min_num_drugs = 4, add_pair_features = T,
                   file_suffix= paste0(paste0(c(10,20), collapse = ','), 'incase'))

JoinedTable <- lst[[1]]
individualgenesother <- JoinedTable %>% select(group_id, group_sym, group_id2, group_sym2) %>%
  filter(group_id == group_id2) %>% distinct()


#####
# get histogram of genes in how many folds during cross validation.

source('import_data.R')
source('glmnetHelperFuns.R')
# CellNum
lst <- import_grpd(average_replicates = T, rmv_cell_count_2 = F, 
                   create_group_files = F, two_conc = c(10,20),
                   min_num_drugs = 4, add_pair_features = T,
                   file_suffix = paste0(paste0(c(10,20), collapse = ','), 'incase'))
JoinedTable <- lst[[1]]

JoinedTable <- JoinedTable %>% group_by(inchi_key) %>% 
  filter(Cpds.conc.uM == max(as.integer(as.character(Cpds.conc.uM)))) %>% ungroup()

JoinedTable <- filter(JoinedTable, !is.na(CellNum))

numbins <- 5
JoinedTable <- BinData(JoinedTable, inchi_key, group_id, num_bins = numbins)

JoinedTable %>% group_by(group_id) %>% 
  mutate(binsum = n_distinct(bin)) %>% 
  select(group_id, binsum) %>% distinct() %>% 
  ggplot() + geom_histogram(aes(x = binsum), bins = 4)

# CellSize
lst <- import_grpd(average_replicates = T, rmv_cell_count_2 = T, 
                   create_group_files = F, two_conc = c(1,2),
                   min_num_drugs = 4, add_pair_features = F,
                   file_suffix = paste0(paste0(c(1,2), collapse = ','), 'incase'))
JoinedTable <- lst[[1]]

JoinedTable <- JoinedTable %>% group_by(inchi_key) %>% 
  filter(Cpds.conc.uM == max(as.integer(as.character(Cpds.conc.uM)))) %>% ungroup()

JoinedTable <- filter(JoinedTable, !is.na(CellNum))

numbins <- 5
JoinedTable <- BinData(JoinedTable, inchi_key, group_id, num_bins = numbins)

JoinedTable %>% group_by(group_id) %>% 
  mutate(binsum = n_distinct(bin)) %>% 
  select(group_id, binsum) %>% distinct() %>% 
  ggplot() + geom_histogram(aes(x = binsum), bins = 4)














