library(tidyverse)

# compare the distributions of the two input matrices, since the primary target 
# model had a higher deviance explained then the all targets model, but there are
# fewer observations in the primary target model so is it reasonable to compare 
# the two models on their deviance explained since both measure have different 
# denominators. Going to check how different those denoms are, etc.. 

source('import_data_baseline_changes.R')
source('glmnetHelperFuns.R')

lst <- import_grpd(average_replicates = T, rmv_cell_count_2 = F,
                   create_group_files = F, two_conc = c(10, 20),
                   min_num_drugs = 4, add_pair_features = T, only_primary_targets = T,
                   file_suffix= paste0(paste0(c(10, 20), collapse = ','), 'PTBaseline'))
# only primary targets
JoinedTable_PT <- lst[[1]]

source('import_data.R')
lst <- import_grpd(average_replicates = T, rmv_cell_count_2 = F,
                   create_group_files = F, two_conc = c(10, 20),
                   min_num_drugs = 4, add_pair_features = T, only_primary_targets = F,
                   file_suffix= paste0(paste0(c(10, 20), collapse = ',')))
# data with all targets
JoinedTable_AT <- lst[[1]]

#1 - (sum(residuals.min[[i]]**2)/sum(Actuals$CellNum**2)
denomPT <- JoinedTable_PT %>% filter(!is.na(CellNum)) %>% .$CellNum %>% unique(.) %>% sum(.**2)
denomAT <- JoinedTable_AT %>% filter(!is.na(CellNum)) %>% .$CellNum %>% unique(.) %>% sum(.**2)

hist(JoinedTable_AT$CellNum %>% unique(.))
hist(JoinedTable_PT$CellNum %>% unique(.))

