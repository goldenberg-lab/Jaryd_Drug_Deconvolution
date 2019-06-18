library(tidyverse)
# Create table of drug to primaryTarg ortho_id to send to Yuan.

tbl <- read_csv('./post fit analysis/Evaluate_Primary_and_off_Target/PrimaryTargetTbl_orthos.csv')

HavePTnoties <- tbl %>% select(inchi_key, ortho_id, primaryTarg) %>% 
  group_by(inchi_key) %>% filter(sum(primaryTarg) == 1) %>% filter(primaryTarg == T) %>%
  select(inchi_key, ortho_id)

HavePTtied <- tbl %>% select(inchi_key, micromolar_activity, ortho_id, primaryTarg) %>%
  group_by(inchi_key) %>% filter(sum(primaryTarg) > 1) %>% filter(primaryTarg == T) %>%
  select(inchi_key, ortho_id)

UnknownPrimary <- tbl %>% select(inchi_key, ortho_id, primaryTarg) %>%
  group_by(inchi_key) %>% filter(sum(primaryTarg) == 0) %>% 
  mutate(ortho_id = NA) %>% select(inchi_key, ortho_id) %>% distinct()

finalTbl <- rbind(HavePTnoties, HavePTtied, UnknownPrimary)
# join the original gene_ids in with the ortho ids.

PlateBook <- read_csv('./Changed_Data/PlateBookID_Orthologs.csv', col_types = 'ddccicccdiicic')
PB <- PlateBook %>%
  select(inchi_key, micromolar_activity, gene_id, ortho_id)
finalTbl <- finalTbl %>% left_join(., PB, by = c('inchi_key', 'ortho_id')) %>% distinct() %>%
  group_by(inchi_key, ortho_id) %>% 
  filter(micromolar_activity == min(micromolar_activity) | is.na(micromolar_activity)) %>%
  distinct() %>% ungroup() %>% select(-ortho_id)

write_csv(finalTbl, './post fit analysis/Evaluate_Primary_and_off_Target/Drug_PrimTarg_Tbl_geneid.csv')

