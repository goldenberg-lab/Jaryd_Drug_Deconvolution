library(tidyverse)
library(Matrix)
library(glmnet)

ModelResults <- tibble(alpha = NA, pairs = NA, intercept = NA, penalty = NA, 
                       conc = NA, `%dev.min` = NA, `HOdev.min` = NA, cor_with_actuals.min = NA,
                       `%dev.1se` = NA, `HOdev.1se` = NA, cor_with_actuals.1se = NA)
# Constants...
min_drugs <- 4
pairs <- T
avg_reps <- T
intercept <- F
conc <- c(10,20)
alpha <- 0
file_suffix <- ''
for(conc in list(c(10,20))){#list(NA, c(1,2), c(3,6), c(10,20))){ #list(c(10,20))){
  source('import_data.R')
  source('glmnetHelperFuns.R')
  lst <- import_grpd(average_replicates = avg_reps, rmv_cell_count_2 = F, 
                     create_group_files = T, two_conc = conc,
                     min_num_drugs = min_drugs, add_pair_features = pairs,
                     file_suffix = file_suffix)
  #save(file = 'glmnetdataMar21st.RData', list = c('lst'))
  #load('glmnetdataMar21st.RData')
  JoinedTable <- lst[[1]]
  
  # initial solution to trial if leaving these drugs in as duplicated 
  # experiments makes a difference. For clarity I should implement this into 
  # import_grpd if we choose to keep it.
  # Remove the "duplicate" experiments, i.e. the same drug applied at two 
  # concentrations in the same bracket i.e. same drug at 10 and 20 or 3 
  # and 6 etc... 
  if (!is.na(conc)){
    JoinedTable <- JoinedTable %>% group_by(inchi_key) %>% filter(Cpds.conc.uM == max(as.integer(as.character(Cpds.conc.uM)))) %>% ungroup()
  }
  
  # Remove all experiments where the phenotype of interest is unknown.
  JoinedTable <- filter(JoinedTable, !is.na(CellNum))
  ################################################################################
  # This function will split the data frames into train and test sets. 
  # 5 bins total, it will separate the drugs into different bins but spread 
  # the genes out so the domain is as even as possible across the bins.
  
  # At this point the function assumes the existence of group_sym, in future
  # should change to take vector of columns to be spread and separated.
  numbins <- 5
  JoinedTable <- BinData(JoinedTable, inchi_key, group_id, num_bins = numbins)
  
  # Rekey the plates after removing missing observations.
  # If statement to handle interaction features or if interactions are present.
  if(names(JoinedTable)[[1]] == 'Plate'){
    Platekey <- JoinedTable %>% group_by(Plate, Well) %>% group_indices()
  } else {
    Platekey <- JoinedTable %>% group_by(inchi_key, Cpds.conc.uM) %>% 
      group_indices()
  }
  
  JoinedTable$Platekey <- Platekey
  ################################################################################
  #Make vector of binarised value for if there is an effect on the gene. Used if
  #including all predictors. 
  
  # using binary indicator for if the gene was inhibited by a drug at a 
  # concentration at least it's IC50 value.
  JoinedTable <- mutate(JoinedTable, value = 1)
  
  #Alternative adding of "values": Use known activities which we classified 
  #as positive of negative influence and set the value vector to negative 1 
  #or positive 1 respectively.
  
  ####
  # If the group has multiple members then it will be given a value of -1, 
  # otherwise it will get the value given by the known activity. Read in the
  # groups created by the import call above and apply the appropriate value
  # to all groups of more than 1 member.
  
  # groups <- read_csv(paste0('./grouping/latestgroups',
  #                           paste0(conc, collapse = ','), '.csv')) %>% select(-num_drugs, -group_num)
  # 
  # orthoGeneMap <- groups %>% select(gene_id, gene_sym, ortho_id, ortho_sym) %>% distinct()
  # 
  # ActDirs <- read_csv('./moa_activities/ActivitiesDirections.csv')
  # ActDirs <- ActDirs %>% unite(activity_moa, starts_with('activity'), sep = ';') %>%
  #   separate(col = activity_moa, into = c('activity_moa', 'extra'), sep = ';NA') %>%
  #   select(-extra, -count, -starts_with('value'), -X15, -comments) %>%
  #   rename(value = `Final Choice`)
  # 
  # PlateBookID <- read_csv('./Changed_Data/PlateBookID_Orthologs.csv')
  # PlateBookID <- rename(PlateBookID, Well = well, Plate = plate)
  # GeneActivityMapping <- as.tibble(PlateBookID) %>% select(inchi_key, gene_id, activity_moa) %>% distinct()
  # 
  # test <- JoinedTable %>% mutate(multimembs1 = group_id %in% groups$group_id, 
  #                                multimembs2 = group_id2 %in% groups$group_id)
  # 
  # test2 <- test %>% left_join(., orthoGeneMap, by = c('group_id' = 'ortho_id'))
  # 
  # # Ideally there would be the same number of rows for test as test3, 
  # # but there are additional rows created 
  # 
  # test3 <- test2 %>% left_join(., GeneActivityMapping, 
  #                              by = c('inchi_key', 'gene_id')) %>% 
  #   filter((!is.na(activity_moa) & multimembs1 == F) | (multimembs1 == T & gene_id == group_id))
  # 
  # test4 <- test3 %>% left_join(., orthoGeneMap, by = c('group_id2' = 'ortho_id'), suffix = c('', '2'))
  # 
  # test5 <- test4 %>% left_join(., GeneActivityMapping, 
  #                              by = c('inchi_key', 'gene_id2' = 'gene_id'), suffix = c('', '2')) %>% 
  #   filter((!is.na(activity_moa2) & multimembs2 == F) | (multimembs2 == T & gene_id2 == group_id2))
  # 
  # cornercase <- ActDirs %>% filter(!is.na(`Drugs if different`))
  # ActDirs <- ActDirs %>% filter(is.na(`Drugs if different`)) %>% select(-`Drugs if different`)
  # 
  # test6 <- test5 %>% left_join(., ActDirs, by = c('activity_moa'))
  # 
  # test7 <- test6 %>% mutate(value = if_else(multimembs1, -1L, value))
  # test8 <- test7 %>% mutate(value = if_else(group_id != group_id2, -1L, value))
  # 
  # JoinedTable <- test8
  ##########################################################################
  #create a key for the group_id.
  GeneKey <- JoinedTable %>% select(group_id, group_id2) %>% distinct()
  #sort by ascending order of group_id
  GeneKey <- arrange(GeneKey, group_id)
  #create key column
  GeneKey <- mutate(GeneKey, genekey = 1:length(GeneKey[[1]]))
  
  #incorporate genekey into JoinedTable
  JoinedTable <- full_join(JoinedTable, GeneKey, by = c("group_id", 'group_id2'))
  
  Predictors <- sparseMatrix(i = JoinedTable$Platekey, j = JoinedTable$genekey,
                             x = JoinedTable$value, use.last.ij = T)
  
  Outcomes <- JoinedTable %>% select(Platekey, CellNum, bin) %>% 
    arrange(Platekey) %>% distinct()
  
  #CellNum sequence
  sequence <- c(seq(1, 2, 0.1), seq(3, 9, 1), seq(10, 100, 10))
  
  #CellSize sequence
  #sequence <- c(seq(1, 2, 0.1), seq(3, 9, 1), seq(10, 90, 10), seq(100, 900, 100))#, seq(1000, 8000, 1000))#, 2^(13:30))
  
  Models <- vector('list', length(sequence))
  names(Models) <- as.character(sequence)
  
  trainingTable <- tibble(Penalty_Factor = NA, MSE = NA, `%Dev` = NA, myDev = NA, numCoef = NA)
  
  # Uncomment for running CellNum
  lambdaseq <- exp(seq(3, 10, length.out = 100))
  
  pdf(paste0('./ModelStats/CellNum/glmnet_Grps_Penalty', '_alpha:', alpha, 
             '_conc:', paste0(conc, collapse = ','), '_Int:', intercept, Sys.Date(),'.pdf'))
  for (p in sequence){
    # add penalty value to all pair features
    GeneKey <- GeneKey %>% mutate(penalty = if_else(group_id != group_id2, p, 1))  
    
    if (alpha == 0){
      fit <- cv.glmnet(Predictors, as.matrix(select(Outcomes, CellNum)),
                       alpha=alpha, family = "gaussian", foldid = Outcomes$bin, 
                       intercept = intercept, penalty.factor = GeneKey$penalty,
                       lambda = lambdaseq)
    } else{
      fit <- cv.glmnet(Predictors, as.matrix(select(Outcomes, CellNum)),
                       alpha=alpha, family = "gaussian", foldid = Outcomes$bin, 
                       intercept = intercept, penalty.factor = GeneKey$penalty)  
    }
    
    plot(fit, main = paste0('Pen:', p, ' alpha:', alpha, 
                            ' conc:', paste0(conc, collapse = ',')))
    
    Models[[as.character(p)]] <- fit
    
    index.min <- which(fit$glmnet.fit$lambda == fit$lambda.min, arr.ind = T)
    
    dev.min <- fit$glmnet.fit$dev.ratio[index.min]
    
    numCoef <- fit$nzero[fit$lambda == fit$lambda.min]
    
    trainingTable <- trainingTable %>% add_row(Penalty_Factor = p, MSE = min(fit$cvm),
                                               `%Dev` = dev.min, numCoef = numCoef)
  }
  
  dev.off()
  trainingTable <- trainingTable %>% filter(!is.na(Penalty_Factor))
  
  PenaltySpace <- matrix(nrow = 0, ncol = 3)
  
  MinCoefFullfit <- matrix(nrow = 0, ncol = 4)
  
  for (penalty in sequence){
    tmp <- cbind(rep(penalty, length(Models[[as.character(penalty)]]$lambda)), 
                 Models[[as.character(penalty)]]$lambda, Models[[as.character(penalty)]]$cvm)
    
    PenaltySpace <- rbind(PenaltySpace, tmp)
    
    tmp2 <- cbind(rep(penalty, length(coef(Models[[as.character(penalty)]], 
                                           'lambda.min')[-1])), 
                  coef(Models[[as.character(penalty)]])[-1],
                  GeneKey$group_id, GeneKey$group_id2)
    
    MinCoefFullfit <- rbind(MinCoefFullfit, tmp2)
  }
  
  PenaltySpace <- as.tibble(PenaltySpace)
  names(PenaltySpace) <- c('penalty', 'lambda', 'MSE')
  
  MinCoefFullfit <- as.tibble(MinCoefFullfit)
  names(MinCoefFullfit) <- c('penalty', 'Coefs', 'group_id', 'group_id2')
  CoefSumMag <- MinCoefFullfit %>% mutate(Pair = group_id != group_id2) %>% 
    group_by(penalty, Pair) %>% summarize(CoefMagnitudeAbs = sum(abs(Coefs)),
                                          CoefMagnitudeSqr = sum(Coefs**2))
  
  
  pdf(paste0('./ModelStats/CellNum/lambdaError', '_alpha:', alpha, '_Int:', intercept,
             '_conc:', paste0(conc, collapse = ','), '_rmvdups', Sys.Date(), '.pdf'))
  print(ggplot(trainingTable) + geom_point(mapping = aes(x = Penalty_Factor, y = MSE)))
  
  print(ggplot(PenaltySpace) + geom_point(mapping = aes(x = penalty, y = lambda, colour = MSE)))
  
  print(ggplot(CoefSumMag) + geom_point(mapping = aes(x = penalty, y = CoefMagnitudeAbs, colour = Pair)))
  
  print(ggplot(CoefSumMag) + geom_point(mapping = aes(x = penalty, y = CoefMagnitudeSqr, colour = Pair)))
  
  dev.off()
  
  minPenalty <- trainingTable[which(trainingTable$MSE == min(trainingTable$MSE)), 'Penalty_Factor']
  OptimalModel <- Models[[as.character(min(minPenalty[[1]]))]]
  
  ##########################################################################
  # checking cv.glmnet i.e. calculating held out deviance explained ("R^2")
  minlambda <- trainingTable %>% filter(MSE == min(trainingTable$MSE))
  
  if (alpha == 0){
    predictions.min <- vector('list', 5)
    predictions.1se <- vector('list', 5)
    residuals.min <- vector('list', 5)
    residuals.1se <- vector('list', 5)
    CoefDet.min <- vector('numeric', 0)
    CoefDet.1se <- vector('numeric', 0)
    for (i in 1:5){
      tmptrain <- JoinedTable %>% filter(bin != i) %>%
        mutate(penalty = if_else(group_id != group_id2, minPenalty[[1]], 1))
      
      newPlatekey <- tmptrain %>% select(Platekey) %>% distinct() %>% arrange(Platekey) %>% mutate(tmpPlatekey = 1:length(.[[1]]))
      
      tmptrain <- tmptrain %>% full_join(., newPlatekey, by = c('Platekey'))
      
      tmpHO <- JoinedTable %>% filter(bin == i)
      
      newPlatekey <- tmpHO %>% select(Platekey) %>% distinct() %>% arrange(Platekey) %>% mutate(tmpPlatekey = 1:length(.[[1]]))
      
      tmpHO <- tmpHO %>% full_join(., newPlatekey, by = c('Platekey'))
      
      tmpPred <- sparseMatrix(i = tmptrain$tmpPlatekey, j = tmptrain$genekey,
                              x = tmptrain$value, use.last.ij = T,
                              dims = c(max(tmptrain$tmpPlatekey), max(JoinedTable$genekey)))
      tmpOutcomes <- tmptrain %>% select(Platekey, CellNum, bin) %>%
        arrange(Platekey) %>% distinct()
      
      tmpPenalty <- tmptrain %>% select(genekey, penalty) %>% arrange(genekey) %>% distinct()
      
      tmpfit <- glmnet(tmpPred, tmpOutcomes$CellNum, family = 'gaussian', 
                       penalty.factor = tmpPenalty$penalty,
                       intercept = intercept, 
                       lambda = exp(seq(log(Models[[as.character(minPenalty[[1]])]]$lambda.1se), 
                                        log(Models[[as.character(minPenalty[[1]])]]$lambda.min), 
                                        length.out = 100)), alpha = alpha)
      
      newx <- sparseMatrix(i = tmpHO$tmpPlatekey, j = tmpHO$genekey,
                           x = tmpHO$value, use.last.ij = T,
                           dims = c(max(tmpHO$tmpPlatekey), max(JoinedTable$genekey)))
      
      predictions.min[[i]] <- predict(tmpfit, newx = newx, s = Models[[as.character(minPenalty[[1]])]]$lambda.min)
      predictions.1se[[i]] <- predict(tmpfit, newx = newx, s = Models[[as.character(minPenalty[[1]])]]$lambda.1se)
      
      Actuals <- tmpHO %>% select(tmpPlatekey, CellNum) %>% arrange(tmpPlatekey) %>% distinct()
      
      residuals.min[[i]] <- Actuals$CellNum - predictions.min[[i]]
      residuals.1se[[i]] <- Actuals$CellNum - predictions.1se[[i]]
      
      CoefDet.min <- c(CoefDet.min, 1 - (sum(residuals.min[[i]]**2)/sum(Actuals$CellNum**2)))
      CoefDet.1se <- c(CoefDet.1se, 1 - (sum(residuals.1se[[i]]**2)/sum(Actuals$CellNum**2)))
    }
  } else {
    CoefDet.min <- ''
    CoefDet.1se <- ''
  }
  
  
  
  ##########################################################################
  
  #recreate dev plot for new Allfit model
  #plot(OptimalModel$glmnet.fit, xvar = "dev", label = FALSE) + title("Final Fraction", line = 3)
  
  #get coefficients for model on All data.
  Allcoefs.min<- coef(OptimalModel$glmnet.fit, s = OptimalModel$lambda.min)
  
  #extract tibbles from coef data
  AllCoefs.min <- as.tibble(summary(Allcoefs.min))
  
  Allcoefs.1se<- coef(OptimalModel$glmnet.fit, s = OptimalModel$lambda.1se)
  
  #extract tibbles from coef data
  AllCoefs.1se <- as.tibble(summary(Allcoefs.1se))
  # summary will always return the coefficients starting from i = 2 because even 
  # when fitting with no intercept there is a placeholder row for the intercept value of 0.
  
  AllCoefs.min <- AllCoefs.min %>% rename(genekey = i) %>% mutate(genekey = genekey - 1) %>% 
    left_join(GeneKey, by = 'genekey') %>%
    left_join(unique(select(JoinedTable, group_id, group_sym)), by = c('group_id')) %>%
    left_join(unique(select(JoinedTable, group_id, group_sym)), by = c('group_id2' = 'group_id'),
              suffix = c('','2')) %>%
    select(-j, -penalty) %>% arrange(desc(x)) %>% filter(genekey != 0) %>%
    mutate(genekey = as.integer(genekey), group_id = as.integer(group_id),
           group_id2 = as.integer(group_id2))
  
  AllCoefs.1se <- AllCoefs.1se %>% rename(genekey = i) %>% mutate(genekey = genekey - 1) %>% 
    left_join(GeneKey, by = 'genekey') %>%
    left_join(unique(select(JoinedTable, group_id, group_sym)), by = c('group_id')) %>%
    left_join(unique(select(JoinedTable, group_id, group_sym)), by = c('group_id2' = 'group_id'),
              suffix = c('','2')) %>%
    select(-j, -penalty) %>% arrange(desc(x)) %>% filter(genekey != 0) %>%
    mutate(genekey = as.integer(genekey), group_id = as.integer(group_id),
           group_id2 = as.integer(group_id2))
  
  
  write_csv(AllCoefs.min, 
            path = paste0('/data/Jaryd/R/CellSummerProject/Coefs/Penalty/CellNum/', 
                          'CellNum','MinDrgs:', min_drugs ,
                          '_alpha:', alpha,'_pairs:',pairs,'_Inter:',
                          intercept,'_avgreps:',avg_reps, '_concs:',
                          paste(conc, collapse = ','),file_suffix,'lambdaMin_rmvdups', Sys.Date(), '.csv'))
  write_csv(AllCoefs.1se, 
            path = paste0('/data/Jaryd/R/CellSummerProject/Coefs/Penalty/CellNum/',
                          'CellNum','MinDrgs:', min_drugs,
                          '_alpha:', alpha,'_pairs:',pairs,'_Inter:',
                          intercept,'_avgreps:',avg_reps, '_concs:',
                          paste(conc, collapse = ','),file_suffix,'lambda1se_rmvdups', Sys.Date(), '.csv'))
  
  #Store the %deviance explained for the model.
  index.min <- which(OptimalModel$glmnet.fit$lambda == OptimalModel$lambda.min, arr.ind = T)
  deviance.min <- OptimalModel$glmnet.fit$dev.ratio[index.min]
  MSEFnl.min <- OptimalModel$cvm[OptimalModel$lambda == OptimalModel$lambda.min]
  
  index.1se <- which(OptimalModel$glmnet.fit$lambda == OptimalModel$lambda.1se, arr.ind = T)
  deviance.1se <- OptimalModel$glmnet.fit$dev.ratio[index.1se]
  MSEFnl.1se <- OptimalModel$cvm[OptimalModel$lambda == OptimalModel$lambda.1se]
  
  Prediction.min <- predict(OptimalModel, Predictors, s = 'lambda.min') %>% 
    as.tibble(.) %>%
    mutate(Platekey = row_number())
  
  PredvsAct.min <- full_join(Prediction.min, Outcomes)
  
  Cor.min <- cor(PredvsAct.min$`1`,PredvsAct.min$CellNum, method = 'spearman')
  
  Prediction.1se <- predict(OptimalModel, Predictors, s = 'lambda.1se') %>% 
    as.tibble(.) %>%
    mutate(Platekey = row_number())
  
  PredvsAct.1se <- full_join(Prediction.1se, Outcomes)
  
  Cor.1se <- cor(PredvsAct.1se$`1`,PredvsAct.1se$CellNum, method = 'spearman')
  
  ModelResults <- ModelResults %>% 
    add_row(alpha = alpha, pairs = pairs, intercept = intercept, penalty = minPenalty[[1]],
            conc = paste(conc, collapse = ', '), `%dev.min` = deviance.min, 
            HOdev.min = paste0(CoefDet.min, collapse = ','),
            cor_with_actuals.min = Cor.min, `%dev.1se` = deviance.1se,
            HOdev.1se = paste0(CoefDet.1se, collapse = ','),
            cor_with_actuals.1se = Cor.1se)
}

# }
#
#   Compare these results before and after adding the activities to see if it 
#   makes a significant difference.
#      write_csv(ModelResults, 'ModelResultsNoActtest.csv')
#      
#      
write_csv(ModelResults, paste0('./ModelStats/CellNum/ModelResultsGrpdfixedPenaltyCellNum_rmvdups', Sys.Date(), file_suffix, '.csv'))
