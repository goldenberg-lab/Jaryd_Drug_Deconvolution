library(tidyverse)

# for CellNum
#####
source('import_data.R')
lst <- import_grpd(rmv_cell_count_2 = F, two_conc = c(10,20), 
                   average_replicates = T, file_suffix = 'Baseline_comparison', 
                   create_group_files = F, create_group_files_long = F,
                   add_pair_features = T)

JoinedTable <- lst[[1]]

BaselineCN <- JoinedTable %>% 
  select(group_id, group_sym, group_id2, group_sym2, CellNum) %>% 
  filter(!is.na(CellNum)) %>%
  distinct() %>% 
  group_by(group_id, group_sym, group_id2, group_sym2) %>%
  summarise_all(mean) %>%
  arrange(desc(CellNum))

write_csv(BaselineCN, './Baseline/Baseline_ordering_Num.csv')

# for CellSize
######
source('import_data.R')
lst <- import_grpd(rmv_cell_count_2 = T, two_conc = c(1,2), 
                   average_replicates = T, file_suffix = 'Baseline_comparison', 
                   create_group_files = F, create_group_files_long = F,
                   add_pair_features = F)

JoinedTable <- lst[[1]]

BaselineCS <- JoinedTable %>% 
  select(group_id, group_sym, CellSize_1) %>% 
  filter(!is.na(CellSize_1)) %>%
  distinct() %>% 
  group_by(group_id, group_sym) %>%
  summarise_all(mean) %>%
  arrange(desc(CellSize_1))

write_csv(BaselineCS, './Baseline/Baseline_ordering_Size.csv')
#####


IndBase <- BaselineCN %>% filter(group_id == group_id2)
PairBase <- BaselineCN %>% filter(group_id != group_id2)
# Only using CellNum coefficients for now, not sure if Shixuan even wants to 
# include cellsize analysis in the final paper at this point.
# idea: calculate Hamming distance between our list of coefficients and this avg list.
coefs <- read_csv(paste0('./Coefs/Penalty/CellNum/CellNumMinDrgs:4_alpha:0_',
                           'pairs:TRUE_Inter:FALSE_avgreps:TRUE_concs:10,20',
                           'lambdaMin_rmvdups2018-08-28.csv'))

Indeffects <- coefs %>% filter(group_id == group_id2)
Paireffects <- coefs %>% filter(group_id != group_id2) %>%
  unite(syms, group_sym, group_sym2) %>% 
  unite(ids, group_id, group_id2)

# probably not the best way to compare, best way is probably looking at the 
# correlation with the validation. But, this returns the percentage of entries 
# that n%/%2 indexs of the same value in y compared to x.
HammingDistlocal <- function(x, y, n = 1){
  if (length(x) != length(y)){
    stop('length(x) must equal length(y)')
  }
  if( n > length(x)){
    stop('n can be at most the length of vector x')
  }
  if ( n > 1){
    y <- c(y[(length(y) - (n%/%2 - 1)):length(y)], y[1:(length(y) - (n%/%2))])  # bugs here
  }
  i = 0
  dist = 0
  while (i < n){
    dist <- dist + length(y) - sum(x != y)
    y <- c(y[2:length(y)], y[1])
    i <- i + 1
  }
  dist/(length(y))
}

HammingDistlocal(c(1,2,3,4,5,6), c(2,4,1,6,5,3), n = 1)

vals <- vector('numeric', length(Indeffects$group_id%/%2))
for (i in 1:length(Indeffects$group_id%/%2)){
  vals[i] <- HammingDistlocal(Indeffects$group_id, IndBase$group_id, i)  
}


