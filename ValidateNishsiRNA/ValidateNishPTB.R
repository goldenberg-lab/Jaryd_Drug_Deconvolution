library(tidyverse)

# calculate correlation with the Nish Validation data for the PTBaseline model.

# Calculate the correlation between the siRNA results and the magnitudes of the
#coefficients found by our linear regression.

# Read in the siRNA dataset from Nish
folder <- '../../CellSummerProjectData/ValidNish/Coulter/'
val_p1 <- read_csv(paste0(folder, '20180801_siRNA-reseed_SL.csv'),
                   col_names = c('exp_num', 'well_group', 'gene_sym',
                                 'num_cells_48h', 'num_cells_72h',
                                 'med_diam_48h', 'med_diam_72h'),
                   col_types = 'ccciidd', skip = 3) %>%
  mutate(exp_num = 1) %>% filter(!is.na(gene_sym))

val_p2 <- read_csv(paste0(folder, '20180809_siRNA-reseed_SL.csv'),
                   col_names = c('exp_num', 'well_group', 'gene_sym',
                                 'num_cells_48h', 'num_cells_72h',
                                 'med_diam_48h', 'med_diam_72h'),
                   col_types = 'ccciidd', skip = 3) %>%
  mutate(exp_num = 2) %>% filter(!is.na(gene_sym))

val <- rbind(val_p1, val_p2) %>% group_by(gene_sym) %>%
  summarise_at(vars(num_cells_48h:med_diam_72h), mean)


grpsPTB <- read_csv('./grouping/latestgroups10,20PTBaseline.csv') %>%
  select(gene_sym, group_sym) %>% dplyr::rename(group_symPT = group_sym)
grps <- read_csv('./grouping/latestgroups10,20CellNum.csv') %>%
  select(gene_sym, group_sym)

val <- left_join(val, grpsPTB, by = 'gene_sym')
val <- left_join(val, grps, by = 'gene_sym')

CoefsPT <- read_csv(paste0('./Coefs/Penalty/CellNum/CellNumMinDrgs:4_alpha:0',
                         '_pairs:TRUE_Inter:FALSE_avgreps:TRUE_concs:10,20',
                         'lambda1se_rmvdups2018-11-12PTBaseline.csv'))

Coefs <- read_csv(paste0('./Coefs/Penalty/CellNum/CellNumMinDrgs:4_alpha:0',
                         '_pairs:TRUE_Inter:FALSE_avgreps:TRUE_concs:10,20',
                         'lambda1se_rmvdups2018-08-28.csv'))
#Get only main effect Coefficients
CoefsPT <- CoefsPT %>% filter(group_id == group_id2)
Coefs <- Coefs %>% filter(group_id == group_id2)

# change CoefsPT columns to have PT suffix
names(CoefsPT) <- paste0(names(CoefsPT), 'PT')

val <- val %>% mutate_at(vars(starts_with('group')), funs(if_else(is.na(.), gene_sym, .)))
all_data <- val

all_data <- all_data %>% left_join(., Coefs, by = c('group_sym')) %>% 
  left_join(., CoefsPT, by = c('group_sym' = 'group_symPT'))

all_data_subset <- all_data %>% filter(!is.na(xPT))


correlationPT48 <- cor.test(all_data_subset$num_cells_48h, all_data_subset$xPT, use = 'complete')
correlationPT72 <- cor.test(all_data_subset$num_cells_72h, all_data_subset$xPT, use = 'complete')

correlation <- cor.test(all_data_subset$num_cells_48h, all_data_subset$x)

sink('./ValidateNishsiRNA/PTBCorrelations.txt')
print('Correlations of Primary Target Models with siRNA experiments:')
correlationPT48
correlationPT72
print('Correlations of Our Model with siRNA experiments:')
correlation
sink()



