library(tidyverse)
library(ggrepel)

# compare predictive performance between our model and baseline model.

Baseline <- read_csv(paste0('./Coefs/Penalty/CellNum/PTBaseline/CellNumMinDrgs:4_alpha:0_',
                            'pairs:TRUE_Inter:FALSE_avgreps:TRUE_concs:10,20',
                            'lambda1se_rmvdups2019-01-09PTBaseline.csv')) %>%
  rename(PT_coef = x) %>% select(-genekey)

#Updated script on June 14th to read the latest coefficient list, which is 
# identical to the previous one imported here, but now it is consistent with 
# other scripts.
Original <- read_csv(paste0('./Coefs/Penalty/CellNum/CellNumMinDrgs:4_alpha:0_',
                            'pairs:TRUE_Inter:FALSE_avgreps:TRUE_concs:10,20',
                            'lambda1se_rmvdups2019-03-25.csv')) %>%
  rename(Orig_coef = x) %>% select(-genekey)

source('./import_data.R')
lst <- import_grpd(rmv_cell_count_2 = F, two_conc = c(10,20), 
                   average_replicates = T, file_suffix = 'delete_me',
                   add_pair_features = T)
JoinedTable <- lst[[1]]

JoinedTable <- JoinedTable %>% group_by(inchi_key) %>%
  # removes duplicate drug trials once from MOA and once from Kinome
  filter(Cpds.conc.uM == max(as.integer(as.character(Cpds.conc.uM)))) %>% 
  ungroup()

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

tmp <- JoinedTable %>% left_join(., Baseline) %>% left_join(., Original) %>%
  group_by(inchi_key) %>% 
  mutate(PT_pred = sum(PT_coef, na.rm = T), Orig_pred = sum(Orig_coef),
         PT_res = CellNum - PT_pred, Orig_res = CellNum - Orig_pred) %>%
  filter(!is.na(PT_coef))

plots <- list()
i <- 1
TBP <- tmp %>% gather(key = "Model", value = "Res", PT_res, Orig_res) %>% select(inchi_key, Model, Res) %>%
  distinct()
plots[[i]] <- ggplot(TBP, aes(x = abs(Res))) + geom_histogram() + facet_wrap(~Model) + xlim(0,5)

TBP <- tmp  %>% select(inchi_key, PT_res, Orig_res) %>% filter(!is.na(PT_res), !is.na(Orig_res)) %>%
  distinct()
# dens would change the colour of points if the density of points was high, a 
# handy function I think, found on Stack exchange.
TBP$dens <- get_density(abs(TBP$PT_res), abs(TBP$Orig_res), n = 200)

# leave the colour out for now, used kernel density estimation to guess the 
# scale, roughly fit with a hyperparameter n = 200, as seen in the call to 
# get_density. Not sure how/if can interpret the density values given but it 
# looks cool, if not super interpretable it does convey the fact that there are 
# more points in a given place than another. just add the colour aesthetic back in, if wanted back.

i <- i + 1
plots[[i]] <- tmp %>% select(inchi_key, CellNum, Orig_res, PT_res) %>% ungroup() %>%
  mutate(tag = if_else(abs(PT_res) > 150, substr(inchi_key, 1, 5), '')) %>%
  distinct() %>%
  ggplot(., aes(x = abs(Orig_res), y = abs(PT_res),
                label = tag)) + geom_point() + 
  geom_text_repel() + geom_abline(slope = 1, intercept = 0)

pdf('./Baseline/Comparison_plots.pdf')
plots
dev.off()


