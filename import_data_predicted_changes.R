# Import and prep data
library(ArgumentCheck)

#Fixed the ordering of the data preprocessing.
import_grpd <- function(keep_all_genes = F, min_num_drugs = 4, 
                        rmv_cell_count_2 = T, colinearity_threshold = 0.5, 
                        two_conc = NA, average_replicates = F, file_suffix = 'Predicted',
                        create_group_files = F, add_pair_features = F,
                        pair_min_num_drugs = 4, only_primary_targets = F,
                        create_group_files_long = F){
  #####
  # keep_all_genes : boolean indicating about keeping genes affected by fewer 
  #               than the min number of drugs
  # min_num_drugs : Minimum number of drugs that must inhibit a gene for it to 
  #               remain a predictor.
  # rmv_cell_count_2 : boolean indicating if removing experiments where the 
  #               population was adversely affected.
  # colinearity_threshold : The value where each pair with a Jaccard Index 
  #               greater than that value is considered to be grouped together.
  # two_conc : NA or a vector of two values where the first value is one of
  #               (1, 3, 10) and the second is one of (2, 6, 20). If NA then all
  #               concentrations will be included, else only the two 
  #               concentrations given will be used.
  # average_replicates : Average the replicate experiments, i.e. if the drug 
  #               applied and the concentration is the same treat those two 
  #               points as a single point
  # file_suffix : string to add to the end of all written files
  # create_group_files : boolean indicating if writing all files in function 
  #               add groups
  # add_pair_features : add features to indicate pairs of genes being inhibited
  # pair_min_num_drugs : the minimum number of drugs the pair features need to 
  #               be inhibited by to be included.
  # only_primary_targets : remove all non-primary tagets for a given drug, if
  #               multiple targets (genes) are targeted are tied for minimum 
  #               threshold then consider them both primary targets.
  # create_group_files_long : boolean to indicate if creating the larger file 
  #               including the jaccard index for all pairs of genes should be 
  #               written, set to T with caution it will increase the run time 
  #               significantly.
  #  DocString
  #####
  
  # Argument Check: needs to be updated for new arguments, using a deprecated 
  # package; should replace with suggested "checkmate" package.
  #####
  Check <- ArgumentCheck::newArgCheck()
  
  if (colinearity_threshold <= 0)
    ArgumentCheck::addError(
      msg = "colinearity_threshold must be greater than 0.",
      argcheck = Check
    )
  if (!is.character(file_suffix))
    ArgumentCheck::addError(
      msg = "file_suffix must be of type character.",
      argcheck = Check
    )
  if (!is.na(two_conc[[1]]) & sum(two_conc %in% c(1,2,3,6,10,20)) != 2){
    ArgumentCheck::addError(
      msg = "Both Concentrations must be concentrations used in the dataset.",
      argcheck = Check
    )
  }
  
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)
  #####
  source('Predicted isoform grouping functions only.R')
  #Import Data, and convert to tibbles
  #load('myParTab.RData')
  load('myParTabedgecorrected.RData')
  PlateBookID <- read_csv('./Changed_Data/PlateBookID_Orthologs.csv')
  genes <- PlateBookID %>% select(gene_id, gene_sym) %>% distinct()
  check <- PlateBookID %>% select(plate, well, inchi_key) %>%
    distinct()
  
  # Label non-human genes with human ortholog.
  orthos <- read_tsv('../../CellSummerProjectData/NCBI/gene_orthologs', 
                     col_types = 'iicii') %>%
    filter(`#tax_id` == 9606) # 9606 == Human, four entries in Other_tax which are human but those genes are not in our data.
  
  Predicted_Targs <- read_csv("../../CellSummerProjectData/Primary_Target_Novartis/Drug_PrimTarg_Tbl_geneid_predictions.csv",
                              col_types = "cid") %>% 
    filter(prediction >= 15) %>% 
    left_join(., orthos, by = c('gene_id' = 'Other_GeneID')) %>% 
    mutate(gene_id = if_else(is.na(GeneID), gene_id, GeneID)) %>%
    group_by(inchi_key, gene_id) %>% 
    summarize(prediction = max(prediction))
  
  PlateBookID <- left_join(check, Predicted_Targs, by = "inchi_key") %>% 
    left_join(., genes, by = "gene_id")
  
  PlateBookID <- rename(PlateBookID, Well = well, Plate = plate)
  PlateBookID <- as.tibble(PlateBookID)
  
  PlateConcentrations <- read.csv("/data/Jaryd/CellSummerProjectData/OriginalData/PlateConcentration_NovartisScreen.csv",
                                  header = T, fill = T, sep = ',')
  PlateConcentrations <- rename(PlateConcentrations, Plate = Plate..)
  PlateConcentrations <- as.tibble(PlateConcentrations)
  
  #Change factor "Dose Response" in Cpds.conc.uM to "Dose_response"
  #gave me problems in with glm recognizing it as a different level.
  levels(PlateConcentrations$Cpds.conc.uM) <- c(levels(PlateConcentrations$Cpds.conc.uM), "Dose_response")
  PlateConcentrations$Cpds.conc.uM[PlateConcentrations$Cpds.conc.uM == "Dose response"] <- "Dose_response"
  PlateConcentrations$Cpds.conc.uM <- factor(PlateConcentrations$Cpds.conc.uM)
  
  PlateBookID <- ungroup(PlateBookID) #Groups carry through when Joining so this
  #is a sanity check to ensure I don't add groups earlier and then it screws me up later.
  
  #create a unique number for each Plate well combo.
  ParTab <- ParTab %>% mutate(Platekey = (Plate-1)*384 + (Well-1))
  
  onlyConcentrations <- select(PlateConcentrations, Plate, Cpds.conc.uM)
  JoinedTable <- full_join(PlateBookID, ParTab)
  JoinedTable <- full_join(JoinedTable, onlyConcentrations)
  
  # Remove experiments where the drug applied is unknown, and remove known 
  # pseudogenes and 'CDC14C' which is missing from the GeneMania database.
  JoinedTable <- JoinedTable %>% filter(!is.na(inchi_key))
  
  if (only_primary_targets){
    primaryTargs <- read_csv('./post fit analysis/Evaluate_Primary_and_off_Target/PrimaryTargetTbl.csv') %>%
      select(inchi_key, group_id, group_sym, primaryTarg)
    JoinedTable <- JoinedTable %>% left_join(., primaryTargs) %>% 
      filter(primaryTarg == T) %>%
      select(Plate, Well, inchi_key, group_id, group_sym, CellNum, CellSize_1, 
             CellCount_1, CellCount_2, Cpds.conc.uM, p1s1_1, p1s1_2, p1s1_3)
  }
  
  # Store the largest population count labeled as adversely affected to be used 
  # later to remove averaged wells which are below this threshold.
  PopThresh <- max((JoinedTable %>% filter(CellCount_2 == 1))$CellCount_1)
  
  # Averaging the Replicate trials to help improve the interpretability of the 
  # model, it should essentially boost all stats by the variance attributable 
  # to the replicates. Do not remove CellCount_2 rows if averaging replicates. 
  # Ignoring the p1s1 outcomes for now because not sure what a mean of those 
  # measurements would tell us. 
  if(average_replicates){
    Platekey <- JoinedTable %>% group_by(inchi_key, Cpds.conc.uM) %>% 
      group_indices()
    
    JoinedTable$Platekey <- Platekey
    
    JoinedTable <- JoinedTable %>% group_by(inchi_key, Cpds.conc.uM) %>%
      mutate_at(vars(starts_with("CellSize"), starts_with("SizeMAD"),
                     starts_with("Fraction"), starts_with('Props'),
                     CellCount_1, CellNum), mean, na.rm = TRUE) %>%
      select(-Plate, -Well, -CellCount_2, -p1s1_1, -p1s1_2, -p1s1_3) %>%
      distinct()
  }
  
  # calculate num_drugs per ortho_id for grouping purposes.
  JoinedTable <- JoinedTable %>% group_by(gene_id) %>% 
    mutate(num_drugs = length(unique(inchi_key))) %>% ungroup()
  
  # call Add_groups function to JoinedTable to merge isoforms into representative groups.
  # Note: This should take place after all experiments that will be removed are 
  # removed, because that could change how separable different features are.
  if (!only_primary_targets){ # primary_targets_tbl already grouped, grouping groups is not good.
    if (create_group_files_long){
      JoinedTable <- AddGroups(JoinedTable, colinearity_threshold, generate_files_long = T, file_suffix,
                               Rmv_CC2 = rmv_cell_count_2, Concentration = two_conc) 
    } else if (create_group_files & !create_group_files_long) {
      JoinedTable <- AddGroups(JoinedTable, colinearity_threshold, generate_files_short = T, file_suffix,
                               Rmv_CC2 = rmv_cell_count_2, Concentration = two_conc)
    } else {
      JoinedTable <- AddGroups(JoinedTable, Rmv_CC2 = rmv_cell_count_2, Concentration = two_conc)
    }
  }
  
  #Remove samples where the population was adversely affected by the drug.
  #In general should be the first thing removed. Do not remove when calculating 
  #scores for Proliferation, i.e. fraction. Threshold for which Wells had a 
  # dramatically adversely affected calculated by Shixuan, I am merely using the
  # largest population labeled as such as the threshold.
  if(rmv_cell_count_2){
    if(average_replicates){
      JoinedTable <- filter(JoinedTable, CellCount_1 >= PopThresh)
    } else{
      JoinedTable <- filter(JoinedTable, CellCount_2 != 1)
    }
  }
  
  # Add paired groups as features
  if (add_pair_features){
    JoinedTable <- JoinedTable %>% 
      select(-starts_with('activity'), 
             -num_drugs, -gene_id, -gene_sym, -prediction) %>% 
      distinct()
    JoinedTable <- create_pair_features(JoinedTable)
  }
  
  # Select a single concentration for each drug
  if (!is.na(two_conc[[1]])){
    JoinedTable <- JoinedTable %>% filter(Cpds.conc.uM %in% two_conc)
  }
  
  # Calculate the number of drugs per group, this should replace previous counts
  # per gene.
  if (add_pair_features){
    JoinedTable <- JoinedTable %>% group_by(group_id, group_id2) %>% 
      mutate(num_drugs = length(unique(inchi_key))) %>% ungroup()
    
    TBP <- JoinedTable %>% filter(group_id != group_id2)
    histogram <- ggplot(TBP) + geom_histogram(mapping = aes(x = num_drugs), binwidth = 1)
    pdf("pairedFeatureHistogram_predicted.pdf")
    print(histogram)
    dev.off()
  } else{
    JoinedTable <- JoinedTable %>% group_by(group_id) %>% 
      mutate(num_drugs = length(unique(inchi_key))) %>% ungroup()
  }
  
  #remove groups only affected by a few drugs.
  if (!keep_all_genes){
    
    # Create histogram plot of num_drugs distribution, for presentation, 
    # does not work if fitting without pairs leave commented, keep to trace 
    # where plot came from for Medicine by design presentation
    # TBP <- JoinedTable %>% select(group_id, group_id2, num_drugs) %>% distinct()
    # 
    # TBPIndiv <- JoinedTable %>% filter(group_id == group_id2, num_drugs != 1) %>% 
    #   select(group_id, group_id2, num_drugs) %>% distinct()
    # 
    # count <- JoinedTable %>% filter(group_id == group_id2 & num_drugs == 1) %>% distinct()
    # 
    # ggplot(TBPIndiv, aes(x = num_drugs)) + geom_histogram(binwidth = 1)+theme(axis.text=element_text(size=16),
    #                                                                      axis.title=element_text(size=20))
    
    if (add_pair_features){
      JoinedTable <- JoinedTable %>% 
                      filter((num_drugs >= min_num_drugs & group_id == group_id2)|
                            (num_drugs >= pair_min_num_drugs & group_id != group_id2))
    } else{
      JoinedTable <- JoinedTable %>% filter(num_drugs >= min_num_drugs)  
    }
  }
  
  lst <- lst(JoinedTable, ParTab)
  
  return(lst)
}

create_pair_features <- function(Table) {
  # create pairs of groups as a new feature, also calculate the number of drugs 
  # that target these new pairs to remove them if they are not inhibited enough 
  # times.
  
  Table <- Table %>% group_by(Platekey) %>% 
                    select(Platekey, group_id, group_sym) %>% 
                    rename(group_id2 = group_id, group_sym2 = group_sym) %>% 
                    full_join(Table, ., by = c('Platekey')) %>% 
                    mutate(group_idcpy = group_id,
                           group_symcpy = group_sym,
                           group_id2cpy = group_id2,
                           group_sym2cpy = group_sym2,
                           group_id = if_else(group_idcpy < group_id2cpy, group_idcpy, group_id2cpy,
                                              missing = group_idcpy),
                           group_sym = if_else(group_idcpy < group_id2cpy, group_symcpy, group_sym2cpy,
                                               missing = group_symcpy),
                           group_id2 = if_else(group_id < group_id2cpy, group_id2cpy, group_idcpy),
                           group_sym2 = if_else(group_id < group_id2cpy, group_sym2cpy, group_symcpy)) %>%
                    select(-group_idcpy,-group_symcpy,-group_id2cpy, -group_sym2cpy,
                           -group_num) %>%
                    distinct()
  
  Table
  
}

