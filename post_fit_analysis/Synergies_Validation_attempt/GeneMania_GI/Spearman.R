library(tidyverse)
library(cowplot)

source('./post fit analysis/Synergies_Validation_attempt/GeneMania_GI/Load_GM_GI_data.R')

known_Synergies <- Synergies_Tbl %>% filter(known == T)

tmp <- cor.test(log10(known_Synergies$Weight), abs(known_Synergies$Syn_coef), method = 'spearman')

nperm <- 1000
permutations <- vector('list', length = nperm)
permutations[[1]] <- known_Synergies %>% select(Weight, Syn_coef)
for(i in 2:nperm){
  permutations[[i]] <- known_Synergies %>% mutate(Weight = sample(Weight, dim(known_Synergies)[[1]])) %>%
    select(Weight, Syn_coef)
}

pvals <- lapply(permutations, function(x){cor.test(x$Weight, abs(x$Syn_coef), method = 'spearman')$p.value})

pval <- (1 + sum(pvals < pvals[[1]]))/(nperm)

TBP <- known_Synergies %>% rename(known_Synergy_confidence = Weight, Synergy_Beta = Syn_coef)
ggplot(TBP, aes(x = abs(Synergy_Beta), y = log10(known_Synergy_confidence))) +
  geom_point() +
  ggtitle("rho: -0.12, Spearman pval: ~8.06E-09, 1000 permutation pval: 0.001")
ggsave('./post fit analysis/Synergies_Validation_attempt/GeneMania_GI/WeightvCoef.pdf')

TBP2 <- Synergies_Tbl %>% 
  mutate(known_Synergy_confidence = if_else(is.na(log10(Weight)), -5, log10(Weight))) %>%
  rename(Synergy_Beta = Syn_coef)
ggplot(TBP2, aes(x = abs(Synergy_Beta), y = known_Synergy_confidence)) +
  geom_point()

# plot some of the permuted weight v coefs.
ggplot(permutations[[2]], aes(x = abs(Syn_coef), y = log10(Weight))) +
  geom_point() +
  ggtitle("permutation_1")
ggplot(permutations[[3]], aes(x = abs(Syn_coef), y = log10(Weight))) +
  geom_point() +
  ggtitle("permutation_2")
ggplot(permutations[[4]], aes(x = abs(Syn_coef), y = log10(Weight))) +
  geom_point() +
  ggtitle("permutation_3")

##### Attempt 2
# Try again removing the less certain coefficients, i.e. keep only the largest 
# 10% of coefficients, and the smallest 10% of coefficients.

thresh <- known_Synergies %>% 
  arrange(desc(abs(Syn_coef))) %>%
  .$Syn_coef %>%
  .[c(0.1*length(.), 0.9*length(.))]

Conf_known_Synergies <- known_Synergies %>%
  filter(abs(Syn_coef) >= abs(thresh[[1]]) | abs(Syn_coef) <= abs(thresh[[2]]))

tmp <- cor.test(Conf_known_Synergies$Weight, Conf_known_Synergies$Syn_coef, method = 'spearman')

ggplot(Conf_known_Synergies, aes(x = abs(Syn_coef), y = log10(Weight))) + geom_point()

# A: doesn't help, the spearman correlation actually went down there are plenty 
# of high confidence synergies which have very small coefficients.

##### Attempt 3
# Only keep known genes with high confidence. Then run permutation test to get a
# p-value, since artificially selecting cut off.
# I feel like this is starting to get sort of data hacky...
lrg_conf_syn <- known_Synergies %>% filter(Weight >= 0.01)

cor <- cor.test(lrg_conf_syn$Weight, abs(lrg_conf_syn$Syn_coef), method = "spearman")
ggplot(lrg_conf_syn, aes(x = abs(Syn_coef), y = Weight)) + geom_point() +
  ggtitle("rho: 0.637, Spearman pval: 0.0045, 1000 permutation pval: 0.001")
ggsave('./post fit analysis/Synergies_Validation_attempt/GeneMania_GI/WeightvCoef_highConfOnly.pdf')

nperm <- 1000
permutations <- vector('list', length = nperm)
permutations[[1]] <- lrg_conf_syn %>% select(Weight, Syn_coef)
for(i in 2:nperm){
  # a permutation is a reordering of the weights for all coefficients and then 
  # filter to the most confident pairing of synergy and Weight.
  permutations[[i]] <- known_Synergies %>% mutate(Weight = sample(Weight, dim(known_Synergies)[[1]])) %>%
    filter(Weight >= 0.01) %>% select(Weight, Syn_coef)
}

rhos <- lapply(permutations, function(x){cor.test(x$Weight, abs(x$Syn_coef), method = 'spearman')$estimate})

pval <- (1 + sum(rhos > rhos[[1]]))/(nperm)

# plot some of the permutations.
ggplot(permutations[[sample(2:nperm, 1)]], aes(x = abs(Syn_coef), y = Weight)) +
  geom_point() +
  ggtitle("permutation")
ggplot(permutations[[sample(2:nperm, 1)]], aes(x = abs(Syn_coef), y = Weight)) +
  geom_point() +
  ggtitle("permutation")
ggplot(permutations[[sample(2:nperm, 1)]], aes(x = abs(Syn_coef), y = Weight)) +
  geom_point() +
  ggtitle("permutation")
















