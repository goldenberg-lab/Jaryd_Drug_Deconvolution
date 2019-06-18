library(tidyverse)
# Create PrimaryTargetTables

Create_primary_target_table <- function(Platebook_path, groups_path){
  # Create and save a table of the primary target labelled for each drug, 
  # assumes the default naming scheme for each table
  # 
  # Platebook_path : String, of the file location to read in the PlatbookID file
  # groups_path : String, of the file location to read in the groups file
  # save_file : String, location to store the resulting csv table
  
  stopifnot(require(tidyverse))
  
  PlateBookID <- read_csv(Platebook_path)
  PlateBookID <- rename(PlateBookID, Well = well, Plate = plate)
  
  PlateBookID <- as.tibble(PlateBookID)
  
  # get mapping of all gene members to their group representative, to remove cases
  # when grouped genes have the same IC50 value for a given drug which could fix 
  # the "ties" problem.
  groups <- read_csv(groups_path)
  Groupmem2Rep <- groups %>% select(gene_id, gene_sym, group_id, group_sym)
  
  # Calculate the Primary Target for each compound
  PrimaryTargs <- PlateBookID %>% 
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
  
  # add the goups to the primary target table, if a member of a group is a primary
  # target then the group is a primary target. 
  PrimaryTargs <- PrimaryTargs %>%
    left_join(., Groupmem2Rep, by = c('ortho_id' = 'gene_id', 'ortho_sym' = 'gene_sym')) %>% 
    mutate(group_id = if_else(is.na(group_id), ortho_id, group_id),
           group_sym = if_else(is.na(group_sym), ortho_sym, group_sym)) %>%
    group_by(inchi_key, group_id) %>%
    # setting micromolar_activity to min and primaryTarg to max to represent if a 
    # group has the primary targeted included in it, then we are considering the 
    # group a primary target, and I will record the IC50 value for the gene inside
    # the group which was small enough to warrant that label.
    mutate(primaryTarg = as.logical(max(primaryTarg)), 
           micromolar_activity = min(micromolar_activity, na.rm = T)) %>%
    select(inchi_key, micromolar_activity, group_id, 
           group_sym, numGenesTargeted, primaryTarg) %>%
    # Need to remove duplicate entries of any group, otherwise it will mess up the
    # sums of coefficients later on.
    distinct()
  
  PrimaryTargs
}

PB_path <- './Changed_Data/PlateBookID_Orthologs.csv'
group10 <- './grouping/latestgroups10,20.csv'
group1 <- './grouping/latestgroups1,2.csv'
save_g10_file <- './post fit analysis/Evaluate_Primary_and_off_Target/PrimaryTargetTbl_10,20.csv'
save_g1_file <- './post fit analysis/Evaluate_Primary_and_off_Target/PrimaryTargetTbl_1,2.csv'
PT_table10 <- Create_primary_target_table(PB_path, group10)
PT_table1 <- Create_primary_target_table(PB_path, group1)
write_csv(PT_table10, save_g10_file)
write_csv(PT_table1, save_g1_file)

# Test out if breaking ties using the predicted affinities from Yuan will help.
Prim_targ_preds <- read_csv('../../CellSummerProjectData/Primary_Target_Novartis/Drug_PrimTarg_Tbl_geneid_predictions.csv',
                            col_types = 'cid')

PT_tbl_wPred <- left_join(PT_table10, Prim_targ_preds, by = c('group_id' = 'gene_id', 'inchi_key'))

# four drugs had no prediction for any of their targets, all four were drugs 
# where there is only 1 target in my data (which could be due to filtering 
# targets inhibited by at least 4 drugs etc...)
# Add column that indicates if our Primary Target matches the predicted one.
PT_tbl_Match_PTPred <- PT_tbl_wPred %>% group_by(inchi_key) %>% 
  mutate(Pred_Prim = if_else(prediction == max(prediction, na.rm = T), T, F))

# four drugs had no prediction for any of their targets
NoPredictions <- PT_tbl_Match_PTPred %>% filter(max(prediction, na.rm = T) == -Inf)
#A: IQFYYKKMVGJFEH-CSMHCCOUSA-N; IWCWQNVIUXZOMJ-MISYRCLQSA-N; LCQLHJZYVOQKHU-VKHMYHEASA-N; COFOABYDLYTAJX-UHFFFAOYSA-N
# All have only 1 target and each has a micromolar_activity of Nan

# Were any drugs not included in the list of drugs sent to Yuan
Sent_Yuan <- read_csv('./post fit analysis/Evaluate_Primary_and_off_Target/Drug_PrimTarg_Tbl_geneid.csv')
setdiff(PT_tbl_Match_PTPred$inchi_key, Sent_Yuan$inchi_key)
#A: No










