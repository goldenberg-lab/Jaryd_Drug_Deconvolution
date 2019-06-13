library(tidyverse)
# Create table which identifies the primary target for each drug.

add_groups_to_primary_targs <- function(PlateBook, groups, file_indicator){
  # This function will take the original mapping of drugs and their targets with
  # each drug target pair's corresponding "IC50", and the groups used to
  # generate the coefficients. All tables assumed to use the standard naming 
  # scheme for this project. File indicator is a string to distinguish separate 
  # callings of this function. Returns 0 if successful.
  
  # Get mapping of group members to the representative
  Groupmem2Rep <- groups %>% select(gene_id, gene_sym, group_id, group_sym)
  
  # calculate the Primary Target for each compound not considering the groups.
  PrimaryTargets <- PlateBook %>% 
    group_by(inchi_key) %>%
    mutate(numGenesTargeted = length(unique(ortho_id))) %>%
    mutate(primaryTarg = if_else(numGenesTargeted == 1, T, 
                                 micromolar_activity == min(micromolar_activity, na.rm = T))) %>%
    mutate(primaryTarg = if_else(is.na(primaryTarg), F, primaryTarg)) %>%
    select(-homol_id, -homol_sym, -Plate, -Well, -main, -origin, -activity_moa, -gene_id, -gene_sym, -taxID) %>%
    distinct()
  
}

# calculate the Primary Target for each compound not considering the groups.
PrimaryTargets <- PlateBook %>% 
  group_by(inchi_key) %>%
  mutate(numGenesTargeted = length(unique(ortho_id))) %>%
  mutate(primaryTarg = if_else(numGenesTargeted == 1, T, 
                               micromolar_activity == min(micromolar_activity, na.rm = T))) %>%
  mutate(primaryTarg = if_else(is.na(primaryTarg), F, primaryTarg)) %>%
  select(-homol_id, -homol_sym, -Plate, -Well, -main, -origin, -activity_moa, -gene_id, -gene_sym, -taxID) %>%
  distinct()

# Save primaryTargets just the ortholog ids, without the grouping being applied 
# post.
write_csv(PrimaryTargs, './post fit analysis/Evaluate_Primary_and_off_Target/PrimaryTargetTbl_orthos.csv')


