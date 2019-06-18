library(tidyverse)

# compare PTBaseline to Original Model
# Create plot of coefficients from Baseline against coefficients from Original 
# Model. Calculate correlation also

PT <- read_csv(paste0('/data/Jaryd/R/CellSummerProject/Coefs/Penalty/CellNum/PTBaseline/',
                'CellNumMinDrgs:4_alpha:0_pairs:TRUE_Inter:FALSE_avgreps:TRUE',
                '_concs:10,20lambda1se_rmvdups2019-01-09PTBaseline_Interp.csv')) %>%
  select(-genekey)

#Updated script on June 14th to read the latest coefficient list, which is 
# identical to the previous one imported here, but now it is consistent with 
# other scripts.
Orig <- read_csv(paste0('/data/Jaryd/R/CellSummerProject/Coefs/Penalty/CellNum/',
                        'CellNumMinDrgs:4_alpha:0_pairs:TRUE_Inter:FALSE_',
                        'avgreps:TRUE_concs:10,20lambda1se_rmvdups2019-03-25_',
                        'Interp.csv')) %>%
  select(-genekey)

Both <- left_join(PT, Orig, by = c('group_id', 'group_id2', 'group_sym', 'group_sym2'), suffix = c('.PT', '.Orig'))

ggplot(Both, aes(x = x.Orig, y = x.PT)) + geom_point()

cor(Both$x.PT, Both$x.Orig, method = 'spearman')















