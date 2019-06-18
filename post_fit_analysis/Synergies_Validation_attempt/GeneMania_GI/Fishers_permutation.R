library(parallel)
library(tidyverse)

# perform enrichment, and correlation analysis on the GeneMania known genetic 
# interactions with our set of synergies.

# call script to load the GM_GI data.
source("./post_fit_analysis/Synergies_Validation_attempt/GeneMania_GI/Load_GM_GI_data.R")

# Prior notes on how I sped up naive implementation of Thresh_Fish
#####
# # There might be potential for speedups since the fisher test is only a couple 
# # milliseconds to run, but I don't think there is any way faster than tabulate 
# # to get all the counts out of the data frame.
# # A sample of tests done to create Faster_TF
# library(microbenchmark)
# tmpdat <- tibble(rnorm(100000))
# tmpdat <- tmpdat %>% mutate(in_set_col = as.logical(runif(nrow(.))))
# microbenchmark::microbenchmark(tmpdat %>% mutate(pos = c(rep(2, 7), rep(4, dim(tmpdat)[[1]] - 7)),
#                                                  pos = pos - in_set_col) %>%
#                                  .$pos %>%
#                                  tabulate(., nbins = 4), 
#                                tmpdat %>% mutate(rnk = seq.int(nrow(.)),
#                                                  pos = if_else(rnk <= 100, 
#                                                                if_else(in_set_col == T, 1, 2),
#                                                                if_else(in_set_col == T, 3, 4)),
#                                                  pos = factor(pos, levels = c(1,2,3,4))) %>%
#                                  .$pos %>%
#                                  table(.))
# 
# # This didn't make any difference to the speed.
# xint <- 1L:100000L
# yint <- rep(1L, 100000)
# ylog <- rep(T, 100000)
# microbenchmark(xint - yint, xint - ylog)
# 
# # check that for a set number of thresholds they produce the same pvalues
# thresholds <- 1:100
# fishers1 <- lapply(thresholds, function(x){print(x);Faster_TF(Synergies_Tbl, x, Syn_coef, known)})
# fishers2 <- lapply(thresholds, function(x){print(x);Thresholded_fishers(Synergies_Tbl, x, Syn_coef, known)})
# 
# pvals1 <- sapply(fishers1, function(x){x$p.value})
# pvals2 <- sapply(fishers2, function(x){x$p.value})
# # the two vectors are identical looks like faster is as correct as the slow 
# # implementation
#####

# Fast implementation of Thresholded Fishers, assumes in_set_vec is presorted.
# 1 if <T & known, 2 if <T & !known, 3 if > T & known, 4 if > T & !known
Thresh_Fish <- function(threshold, in_set_vec){
  {c(rep(2L, threshold), rep(4L, length(in_set_vec) - threshold)) - in_set_vec} %>%
    tabulate(., nbins = 4) %>%
    matrix(., ncol = 2) %>%
    fisher.test(., alternative = 'greater')
}

# Calculates the pvalue for each permutation given.
Perm_Fish <- function(threshold, in_set_vecs){
  pvals <- numeric(length = length(in_set_vecs))
  for (i in 1:length(in_set_vecs)){
    pvals[[i]] <- Thresh_Fish(threshold, in_set_vecs[[i]])$p.value
  }
  pvals
}

thresholds <- 1:dim(Synergies_Tbl)[[1]]
nperm <- 1000
permutations <- vector("list", length = nperm)
permutations[[1]] <- arrange(Synergies_Tbl, desc(abs(Syn_coef))) %>% .$known
for (i in 2:nperm){
  # Permutation method which doesn't preserve the fact that groups get the same 
  # coefficient.
  #permutations[[i]] <- sample(permutations[[1]], length(permutations[[1]]))
  # Permutation method which will preserve the groups.
  permutations[[i]] <- Synergies_Tbl %>% select(Syn_coef) %>%
    distinct() %>% mutate(new_order = sample(Syn_coef, nrow(.))) %>%
    left_join(Synergies_Tbl, ., by = "Syn_coef") %>% 
    arrange(desc(abs(new_order))) %>% 
    .$new_order
}
names(permutations) <- paste0('Starting Permutation: ', 1:nperm)
fisher_pvals <- mapply(permutations, names(permutations), FUN = 
                           function(p, id){print(id); 
                             sapply(thresholds, 
                                    function(t){if(t%%1000 == 0){print(t)}; 
                                      Thresh_Fish(t, p)$estimate})})#$p.value})})#,
                         #mc.cores = 7)

# get the pvalue given the list of fisher.test results, across thresholds, 
# across permutations.
best_pvals <- sapply(as.data.frame(fisher_pvals), min)

pval <- (1 + sum(best_pvals < best_pvals[[1]]))/(nperm)

save.image('./post fit analysis/Synergies_Validation_attempt/GeneMania_GI/Fishers_results_test_greater_fixed_permutation.RData')

# plot a random ten permutations, pvalues over thresholds.
set.seed(7)
as_tibble(fisher_pvals) %>% select(c(1, sample(1:(ncol(.)-1), 8))) %>%
  mutate(., threshold = 1:nrow(.)) %>%
  gather(Permutation, p_value, starts_with('Starting')) %>%
  ggplot(., aes(x = threshold, y = p_value)) + geom_point() + 
  facet_wrap('~Permutation') +
  ggtitle('Sample of fisher p_values for every threshold on true ordering and random permutations')
ggsave('./post fit analysis/Synergies_Validation_attempt/GeneMania_GI/Fisher_thresholds_fixed_permutation.pdf')
# function to calculate a ratio for each index in a logical vector how many True
# values are before that point of the total number of True values.
#   x : logical vector
# roll_ratio <- function(x){
#   after <- sum(x) # inital values, for the zeroth index.
#   before <- 0
#   sapply(x, function(val){before <<- before + val; before/after})
# }
# 
# Synergies_Tbl %>% arrange(desc(abs(Syn_coef))) %>%
#   mutate(rank = 1:nrow(.), known_ratio = roll_ratio(known)) %>%
#   ggplot(., aes(x = rank, y = known_ratio, colour = known)) +
#   geom_point() + geom_abline(intercept = 0, slope = 1/nrow(Synergies_Tbl))

