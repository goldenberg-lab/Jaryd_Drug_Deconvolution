library(tidyverse)

# Helper function for running glmnet

BinData <- function(Data, Separate, Spread, num_bins = 5){
  # This function will take two variables as input and try to separate one 
  # variable as much as possible between the bins and spread the other variable 
  # out as much as possible across num_bins bins. Plan to accomodate multiple columns passed 
  # in for Separate or Spread, but they must be redundant columns, i.e. for all 
  # unique elements in the first column the second column does not distinguish 
  # between them though it can group them. (In my case this refers to how 
  # gene_id is unique but gene_sym is a more human friendly name, but for every 
  # unique gene_id there is a unique gene_symbol).
  
  Separate_Q <- enquo(Separate)
  Spread_Q <- enquo(Spread)
  #step one : Get number of times a "gene_group" was inhibited by any drug in the data.
  tmp <- Data %>% select(!!Separate_Q, !!Spread_Q) %>% distinct()
  GeneInhibitedcounts <- tmp %>% count(!!Spread_Q)
  
  #Join that information into a temporary object
  Test <- full_join(Data, GeneInhibitedcounts, by = c(quo_name(Spread_Q)))
  #Set a representative Gene affected by each drug as the one inhibited the least.
  Test <- Test %>% select(Platekey, !!Separate_Q, !!Spread_Q, group_sym, n) %>%
    group_by(!!Separate_Q) %>%
    filter(n == min(n)) %>%
    mutate(id = row_number()) %>%
    filter(id == min(id))
  # Set internal id to count up for each group, i.e. the first drug that targets geneA etc.
  Test <- Test %>% group_by(!!Spread_Q) %>% mutate(internal_id = row_number()%%num_bins)
  # external id gives each group a unique number.
  Test <- Test %>% ungroup() %>% mutate(., external_id = group_indices(., !!Spread_Q))
  # summing these two ids will spread out the data to be evenly dispersed amongst 
  # the bins and ensures that for a given gene it is represented in the most bins. 
  # Modulo 5 is used to force them into the bins
  Test <- Test %>% mutate(bin = (internal_id + external_id)%%num_bins + 1)
  #ungroup Test so joining can drop grouping variables
  Test <- Test %>% ungroup()
  #plot to see distribution in bins, not perfectly even but very close.
  ggplot(Test) + geom_histogram(mapping = aes(x = bin), binwidth = 1)
  #Join temp object back to JoinedTable to have those bin numbers
  #included per line in original data.
  Data <- select(Test, !!Separate_Q, bin) %>%
    full_join(Data, ., by = c(quo_name(Separate_Q)))
  
  Data
}