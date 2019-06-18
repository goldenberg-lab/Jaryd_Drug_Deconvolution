library(tidyverse)

# Label coefficients as unintepretable/interpretable groups by the definitions 
# agreed upon with Shixuan.
 
# Read in latest coefficients, change if new coefficients are calculated.

dir <- './Coefs/Penalty/'

filepaths <- list.files(dir, recursive = T)
filepaths <- filepaths[grep('.csv', filepaths)]
filepaths <- filepaths[grep(paste0('(old/)|',
                                   '(alpha:1)|',
                                   '(Inter:T)|',
                                   '(CellSize.*pairs:TRUE)|',
                                   '(Interp)|',
                                   '(NoPairs_for_comparison)'), filepaths, invert = T)]
#filepaths <- filepaths[grep('lambdaMin', filepaths)]
filepaths <- paste0(dir, filepaths)

Coefs <- read_csv(filepaths[1]) %>% mutate(Model = filepaths[1])
for (file in filepaths[-1]){
  if (!is_empty(grep('CellNum', file))){
    tmp <- read_csv(file, col_types = 'idiicc') %>% 
      mutate(Model = file)
  } else{
    tmp <- read_csv(file, col_types = 'idic') %>% 
      mutate(Model = file, group_id2 = NA, group_sym2 = NA)
  }
  
  Coefs <- rbind(Coefs, tmp)
}


# Read in Groups
# Define the groups that are interpretable based on the following:
# If there is at least one member added to the group due to inseparability 
# instead of a known physical interaction or shared protein domain.
grps1020 <- read_csv('./grouping/GroupingsByEdges10,20.csv', col_types = 'cicidlci')
grps36 <- read_csv('./grouping/GroupingsByEdges3,6.csv', col_types = 'cicidlci')
grps12 <- read_csv('./grouping/GroupingsByEdges1,2.csv', col_types = 'cicidlci')
grpsNA <- read_csv('./grouping/GroupingsByEdgesNA.csv', col_types = 'cicidlci')

# there is one exception Shixuan found which is group_id 8313, these two genes with no 
# connection in GeneMania are in the same pathway so we can interpret this group.
grps1020 <- grps1020 %>% group_by(group_sym, group_id) %>%
  mutate(numMemb = n_distinct(c(gene_id1, gene_id2)), 
         NoBioRelation = F %in% GMConnection & group_id != 8313, 
         #NoBioRelation implies that there is at least one edge without a 
         #connection in Gene Mania and it is not the one group exception 8313 
         #which Shixuan found manually.
         isinterp1020 = !NoBioRelation & numMemb <= 20) %>%
  ungroup()

grps12 <- grps12 %>% group_by(group_sym, group_id) %>%
  mutate(numMemb = n_distinct(c(gene_id1, gene_id2)),
         NoBioRelation = F %in% GMConnection,
         isinterp12 = !NoBioRelation & numMemb <= 20) %>%
  ungroup()

join_interpretable <- function(coef_tbl, grps, interpcol){
  interpcol <- enquo(interpcol)
  print(interpcol)
  tmp <- left_join(coef_tbl, select(grps, group_id, !! interpcol) %>% distinct(), by = c('group_id'))
  left_join(tmp, select(grps, group_id, !! interpcol) %>% distinct(), 
            suffix = c('.1', '.2'), by = c('group_id2' = 'group_id'))
}

Coefs <- join_interpretable(Coefs, grps1020, isinterp1020)
Coefs <- join_interpretable(Coefs, grps12, isinterp12)

# saves the Coefficients for each model as a separate table
# quick rules for applying interp at this stage:
# NA means it's a gene that is in a group of size 1 thus it is interpretable
# For pairs iff both groups in the pair are interpretable then the pair is interpretable
# O.W. it takes on the value of isinterp as defined by the grouping section above.
test2 <- Coefs %>%
  mutate_at(vars(starts_with('isinterp')), function(x) {if_else(is.na(x), T, x)}) %>%
  mutate(isinterp1020 = isinterp1020.1 & isinterp1020.2,
         isinterp12 = isinterp12.1 & isinterp12.2) %>%
  select(-matches('\\.[12]')) %>%
  nest(-Model) %>%
  pwalk(function(Model, data){
    fname <- paste0(str_split(Model, '.csv', 2, T)[1], '_Interp.csv')
    if (grepl('CellNum', fname)){
      write_csv(data %>% select(-isinterp12), fname)
    } else {
      write_csv(data %>% select(-isinterp1020), fname)
    }
  })
