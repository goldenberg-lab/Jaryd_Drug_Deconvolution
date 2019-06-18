library(tidyverse)
library(Matrix)
library(glmnet)

ModelResults <- tibble(alpha = NA, pairs = NA, intercept = NA, penalty = NA, 
                       conc = NA, `%dev.min` = NA, `HOdev.min` = NA, cor_with_actuals.min = NA,
                       `%dev.1se` = NA, `HOdev.1se` = NA, cor_with_actuals.1se = NA)
# Constants...
min_drugs <- 4
# alpha <- 0
pairs <- F
# # Grouped <- T # set up by different script
avg_reps <- T
intercept <- T
conc <- c(1,2)
alpha <- 0
new_value <- ''
intercept <- F
for(conc in list(c(1,2), c(3,6), c(10,20))){
  source('import_data_predicted_changes.R')
  source('glmnetHelperFuns.R')
  # added CellSize_1 suffix to double check the groupings are always being 
  # calculated the same if doubted.
  lst <- import_grpd(average_replicates = avg_reps, rmv_cell_count_2 = T, 
                     create_group_files = T, two_conc = conc,
                     min_num_drugs = min_drugs, add_pair_features = pairs,
                     file_suffix = paste0(paste0(conc, collapse = ','), 'CellSize_1'))
  #save(file = 'glmnetdataMar21st.RData', list = c('lst'))
  #load('glmnetdataMar21st.RData')
  JoinedTable <- lst[[1]]
  
  # Remove all experiments where the phenotype of interest is unknown.
  JoinedTable <- filter(JoinedTable, !is.na(CellSize_1))
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
  
  ##########################################################################
  #create a key for the group_id.
  GeneKey <- JoinedTable %>% select(group_id) %>% distinct()
  #sort by ascending order of group_id
  GeneKey <- arrange(GeneKey, group_id)
  #create key column
  GeneKey <- mutate(GeneKey, genekey = 1:length(GeneKey[[1]]))
  
  #incorporate genekey into JoinedTable
  JoinedTable <- full_join(JoinedTable, GeneKey, by = c("group_id"))
  
  Predictors <- sparseMatrix(i = JoinedTable$Platekey, j = JoinedTable$genekey,
                             x = JoinedTable$value, use.last.ij = T)
  
  Outcomes <- JoinedTable %>% select(Platekey, CellSize_1, bin) %>% 
    arrange(Platekey) %>% distinct()
  
  Model <- cv.glmnet(x = Predictors, y = Outcomes$CellSize_1,
                     alpha = alpha, family = "gaussian", foldid = Outcomes$bin,
                     intercept = intercept)
  
  pdf(paste0('./ModelStats/CellSize/glmnetplot_Conc:', 
             paste0(conc, collapse = ''), '_Int:', intercept, '.pdf'))
  plot(Model)
  dev.off()
  
  ##########################################################################
  # checking cv.glmnet i.e. calculating held out deviance explained (R^2)
  minlambda <- Model$lambda.min
  
  if (alpha == 0){
    predictions.min <- vector('list', 5)
    predictions.1se <- vector('list', 5)
    residuals.min <- vector('list', 5)
    residuals.1se <- vector('list', 5)
    CoefDet.min <- vector('numeric', 0)
    CoefDet.1se <- vector('numeric', 0)
    for (i in 1:5){
      tmptrain <- JoinedTable %>% filter(bin != i)
      
      newPlatekey <- tmptrain %>% select(Platekey) %>% distinct() %>% arrange(Platekey) %>% mutate(tmpPlatekey = 1:length(.[[1]]))
      
      tmptrain <- tmptrain %>% full_join(., newPlatekey, by = c('Platekey'))
      
      tmpHO <- JoinedTable %>% filter(bin == i)
      
      newPlatekey <- tmpHO %>% select(Platekey) %>% distinct() %>% arrange(Platekey) %>% mutate(tmpPlatekey = 1:length(.[[1]]))
      
      tmpHO <- tmpHO %>% full_join(., newPlatekey, by = c('Platekey'))
      
      tmpPred <- sparseMatrix(i = tmptrain$tmpPlatekey, j = tmptrain$genekey,
                              x = tmptrain$value, use.last.ij = T,
                              dims = c(max(tmptrain$tmpPlatekey), max(JoinedTable$genekey)))
      tmpOutcomes <- tmptrain %>% select(Platekey, CellSize_1, bin) %>%
        arrange(Platekey) %>% distinct()
      
      tmpPenalty <- tmptrain %>% select(genekey) %>% arrange(genekey) %>% distinct()
      
      tmpfit <- glmnet(tmpPred, tmpOutcomes$CellSize_1, family = 'gaussian', 
                       intercept = intercept, 
                       lambda = exp(seq(log(Model$lambda.1se), 
                                        log(Model$lambda.min), 
                                        length.out = 100)), alpha = alpha)
      
      newx <- sparseMatrix(i = tmpHO$tmpPlatekey, j = tmpHO$genekey,
                           x = tmpHO$value, use.last.ij = T,
                           dims = c(max(tmpHO$tmpPlatekey), max(JoinedTable$genekey)))
      
      predictions.min[[i]] <- predict(tmpfit, newx = newx, s = Model$lambda.min)
      predictions.1se[[i]] <- predict(tmpfit, newx = newx, s = Model$lambda.1se)
      
      Actuals <- tmpHO %>% select(tmpPlatekey, CellSize_1) %>% arrange(tmpPlatekey) %>% distinct()
      
      residuals.min[[i]] <- Actuals$CellSize_1 - predictions.min[[i]]
      residuals.1se[[i]] <- Actuals$CellSize_1 - predictions.1se[[i]]
      
      CoefDet.min <- c(CoefDet.min, 1 - (sum(residuals.min[[i]]**2)/sum(Actuals$CellSize_1**2)))
      CoefDet.1se <- c(CoefDet.1se, 1 - (sum(residuals.1se[[i]]**2)/sum(Actuals$CellSize_1**2)))
    }
  } else {
    CoefDet.min <- ''
    CoefDet.1se <- ''
  }
  
  ##########################################################################
  
  #recreate dev plot for new Allfit model
  #plot(Model$glmnet.fit, xvar = "dev", label = FALSE) + title("Final Fraction", line = 3)
  
  #get coefficients for model on All data.
  Allcoefs.min<- coef(Model$glmnet.fit, s = Model$lambda.min)
  
  #extract tibbles from coef data
  AllCoefs.min <- as.tibble(summary(Allcoefs.min))
  
  Allcoefs.1se<- coef(Model$glmnet.fit, s = Model$lambda.1se)
  
  #extract tibbles from coef data
  AllCoefs.1se <- as.tibble(summary(Allcoefs.1se))
  # summary will always return the coefficients starting from i = 2 because even 
  # when fitting with no intercept there is a placeholder row for the intercept value of 0.
  
  AllCoefs.min <- AllCoefs.min %>% rename(genekey = i) %>% mutate(genekey = genekey - 1) %>% 
    left_join(GeneKey, by = 'genekey') %>%
    left_join(unique(select(JoinedTable, group_id, group_sym)), by = c('group_id')) %>%
    select(-j) %>% arrange(desc(x)) %>% filter(genekey != 0) %>%
    mutate(genekey = as.integer(genekey), group_id = as.integer(group_id))
  
  AllCoefs.1se <- AllCoefs.1se %>% rename(genekey = i) %>% mutate(genekey = genekey - 1) %>% 
    left_join(GeneKey, by = 'genekey') %>%
    left_join(unique(select(JoinedTable, group_id, group_sym)), by = c('group_id')) %>%
    select(-j) %>% arrange(desc(x)) %>% filter(genekey != 0) %>%
    mutate(genekey = as.integer(genekey), group_id = as.integer(group_id))
  
  write_csv(AllCoefs.min, 
            path = paste0('/data/Jaryd/R/CellSummerProject/Coefs/Penalty/CellSize/Predicted/', 
                          'CellSize_1','MinDrgs:', min_drugs ,
                          '_alpha:', alpha,'_pairs:',pairs,'_Inter:',
                          intercept,'_avgreps:',avg_reps, '_concs:',
                          paste(conc, collapse = ','),new_value,'lambdaMin',
                          Sys.Date(),'Predicted.csv'))
  write_csv(AllCoefs.1se, 
            path = paste0('/data/Jaryd/R/CellSummerProject/Coefs/Penalty/CellSize/Predicted/',
                          'CellSize_1','MinDrgs:', min_drugs,
                          '_alpha:', alpha,'_pairs:',pairs,'_Inter:',
                          intercept,'_avgreps:',avg_reps, '_concs:',
                          paste(conc, collapse = ','),new_value,'lambda1se',
                          Sys.Date(),'Predicted.csv'))
  
  #Store the %deviance explained for the model.
  index.min <- which(Model$glmnet.fit$lambda == Model$lambda.min, arr.ind = T)
  deviance.min <- Model$glmnet.fit$dev.ratio[index.min]
  MSEFnl.min <- Model$cvm[Model$lambda == Model$lambda.min]
  
  index.1se <- which(Model$glmnet.fit$lambda == Model$lambda.1se, arr.ind = T)
  deviance.1se <- Model$glmnet.fit$dev.ratio[index.1se]
  MSEFnl.1se <- Model$cvm[Model$lambda == Model$lambda.1se]
  
  Prediction.min <- predict(Model, Predictors, s = 'lambda.min') %>% 
    as.tibble(.) %>%
    mutate(Platekey = row_number())
  
  PredvsAct.min <- full_join(Prediction.min, Outcomes)
  
  Cor.min <- cor(PredvsAct.min$`1`,PredvsAct.min$CellSize_1, method = 'spearman')
  
  Prediction.1se <- predict(Model, Predictors, s = 'lambda.1se') %>% 
    as.tibble(.) %>%
    mutate(Platekey = row_number())
  
  PredvsAct.1se <- full_join(Prediction.1se, Outcomes)
  
  Cor.1se <- cor(PredvsAct.1se$`1`,PredvsAct.1se$CellSize_1, method = 'spearman')
  
  ModelResults <- ModelResults %>% 
    add_row(alpha = alpha, pairs = pairs, intercept = intercept,
            conc = paste(conc, collapse = ', '), `%dev.min` = deviance.min, 
            HOdev.min = paste0(CoefDet.min, collapse = ','),
            cor_with_actuals.min = Cor.min, `%dev.1se` = deviance.1se,
            HOdev.1se = paste0(CoefDet.1se, collapse = ','),
            cor_with_actuals.1se = Cor.1se)
}
#}
# }
#
#   Compare these results before and after adding the activities to see if it 
#   makes a significant difference.
#      write_csv(ModelResults, 'ModelResultsNoActtest.csv')
#      
#      
write_csv(ModelResults, paste0('./ModelStats/CellSize/ModelResultsGrpdfixedPenaltyCellSize_Predicted', Sys.Date(),'.csv'))
