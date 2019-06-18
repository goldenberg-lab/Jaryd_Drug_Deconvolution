library(tidyverse)
library(cowplot)
# Create waterfall plots of the pIC50 = -log_10(IC50) against the coefficient for each drug that
# is classified as well predicted.

base_dir <- './post_fit_analysis/Evaluate_Primary_and_off_Target/'
fs <- list.files(base_dir)
fs <- fs[grepl("^FeaturePerturbance", fs)]
# extract the date from each file, known beforehand the date is the last part 
# before the file extension separated by underscores.
dates <- fs %>% 
  str_split(., '_') %>% 
  sapply(., function(x) {tail(x, 1)}) %>%
  str_split(., '[.]') %>%
  sapply(., function(x) {head(x, 1)}) %>%
  as.Date.character(., "%b-%d-%Y")

f <- fs[which(dates == max(dates))]
if (length(f) != 1){
  stop("more than one file had the latest date, uncertain which to read in.")
}

drugTargTbl <- read_csv(paste0(base_dir, f), col_types = 'cicicddiidldddddi')

# remove any drug where at least one targets IC50 is unknown (assumed small), and
# if it was a drug that was not "PredictedWell"
PredWellDrug <- drugTargTbl %>% group_by(inchi_key) %>%
  filter(sum(micromolar_activity) != Inf & !is.na(sum(PredWell)) & group_id == group_id2)

# Get chemical name from pubchem database.
ChemNames <- read_csv2('../../CellSummerProjectData/PubChem/out.csv') # some chem names have commas so using ; as sep
tmp <- function(x) {
  logic <- grepl('^.*[[:digit:][:punct:][:blank:]].*$', x)
  x[logic] <- ''
  x
}
ChemNames <- ChemNames %>% mutate_at(vars(starts_with('X')), tmp)

# Replace all strings that start with "old_prefix" with the 
# `new_prefix`separator`idx`, where idx is how many new_prefixes have been added
# so far.
prefix_replacement <- function(prev_name, cur_name, replace_prefix, separator){
  if (!grepl(paste0("^", replace_prefix), cur_name)){
    return(c(cur_name, prev_name))
  }
  tmp <- str_split(prev_name, separator)
  rep(paste0(tmp[[1]][1], separator, as.numeric(tmp[[1]][2]) + 1), 2)
}
get_replace_fun <- function(replace_prefix, separator){
  replacement <- "name_0"
  function(prev_name, cur_name){
    state <- prefix_replacement(replacement, cur_name, replace_prefix, separator)
    replacement <<- state[[2]]
    state[[1]]  
  }
}
replace_X <- get_replace_fun("X", "_")
  
names(ChemNames) <- get_replace_fun("X", "_") %>% Reduce(names(ChemNames), accumulate = T)

tChemNames <- as_tibble(cbind(t(ChemNames %>% select(-inchi_key))))
names(tChemNames) <- ChemNames$inchi_key


get_name <- function(x){
  name <- x[grepl("^[A-Z][a-z]", x)]
  if (length(name) < 1){
    name <- x[grepl("^[A-Z]", x)]
  }
  if (length(name) < 1){
    name <- x[grepl("^[a-z]", x)]
  }
  if (length(name) <1){
    name <- NA
  }
  name[[1]]
}
Names <- tChemNames %>% mutate_all(get_name) %>% .[1,] %>% as_tibble(.) %>% 
  gather(key = "inchi_key", value = "drg_name")

PredWellDrug <- left_join(PredWellDrug, Names, by = "inchi_key") %>% 
  mutate(drg_name = if_else(is.na(drg_name), inchi_key, drg_name))

# Pick random 20 drugs for preliminary plots, of groups with at least 5 members
rdm20 <- sample(unique(PredWellDrug$inchi_key), 20) 
TBP <- PredWellDrug %>% filter(inchi_key %in% rdm20 & group_id == group_id2)

ggplot(TBP, aes(x = x, y = -log10(micromolar_activity), colour = primaryTarg,
                size = ProportionPrediction)) + 
  geom_point(shape = 13) + 
  facet_wrap(~drg_name) +
  scale_color_manual(values=c("#000000", "#ff0000"))
ggsave('./post fit analysis/Waterfall draft.png')

# Interesting inchi_keys (part of): BIIVYFLTOXDAOV-YVEFUNNKSA-N, RTIXKCRFFJGDFG-UHFFFAOYSA-N,
# AHJRHEGDXFFMBM-UHFFFAOYSA-N, GGPZCOONYBPZEW-UHFFFAOYSA-N, DKXHSOUZPMHNIZ-UHFFFAOYSA-N,
# LCNDUGHNYMJGIW-UHFFFAOYSA-N

# Manually selected Drugs
selected <- c('BIIVYFLTOXDAOV-YVEFUNNKSA-N', 'RTIXKCRFFJGDFG-UHFFFAOYSA-N',
              'AHJRHEGDXFFMBM-UHFFFAOYSA-N', 'GGPZCOONYBPZEW-UHFFFAOYSA-N', 'DKXHSOUZPMHNIZ-UHFFFAOYSA-N',
              'LCNDUGHNYMJGIW-UHFFFAOYSA-N')

TBP <- PredWellDrug %>% filter(inchi_key %in% selected & group_id == group_id2)
ggplot(TBP, aes(x = x, y = -log10(micromolar_activity), colour = primaryTarg,
                size = ProportionPrediction)) + 
  geom_point(shape = 13) + 
  facet_wrap(~drg_name) +
  scale_color_manual(values=c("#000000", "#ff0000"))
# Get the distribution of all correlations between the IC50 and the Coefficient
check <- PredWellDrug %>% 
  filter(!is.na(PredWell)) %>% group_by(inchi_key) %>% 
  mutate(pcor = cor(x, micromolar_activity), 
         scor = cor(x, micromolar_activity, method = 'spearman'), 
         numGenesTargeted = n()) %>% 
  select(inchi_key, pcor, scor, numGenesTargeted) %>% 
  distinct()
ggsave("./post_fit_analysis/Chosen_IC50vEff_plts.pdf")

