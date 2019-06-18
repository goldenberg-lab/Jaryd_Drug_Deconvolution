library(tidyverse)
# genomeRNAi database files exploration.

dir <- "../../CellSummerProjectData/genomernai/proliferation/"
fs <- paste0(dir, list.files(dir))
tbls <- lapply(fs, function(x) {read_tsv(x, skip = 19)})

# Note two of them don't have scores, one has the scores truncated to just the 
# ones with a score less than -2 and the other seemingly has all scores reported.
#####
# initially comparing the scores from fs[[2]] and the coefficients from Cell 
# Proliferation, to see if similar things were found. Since it has actual scores 
# we can compare the scores together.

coefs <- read_csv('./Coefs/Penalty/CellNum/CellNumMinDrgs:4_alpha:0_pairs:TRUE_Inter:FALSE_avgreps:TRUE_concs:10,20lambda1se_rmvdups2018-12-19_Interp.csv') %>%
  filter(isInterp & group_id == group_id2) %>% 
  arrange(desc(abs(x))) %>% 
  mutate(MWgroup = c(rep('A', 100), rep('B', length(x) - 100)))

study2 <- tbls[[2]]
study_scores <- study2 %>% #select(`Entrez ID`, Score) %>% 
  mutate(`Entrez ID` = as.numeric(`Entrez ID`)) %>% 
  filter(!is.na(`Entrez ID`))

test <- coefs %>% left_join(., study_scores, by = c('group_id' = 'Entrez ID'))

# cor as in Nish experiments
cor.test(test$x, test$Score, use = 'complete.obs' )
plot(test$x, test$Score, use = 'complete.obs')

# Mann whitney U test as planned for Novartis validation
MWtest <- test %>% filter(!is.na(Score)) %>% mutate(abs_x = abs(x), abs_score = abs(Score))

wilcox.test(formula = abs_score ~ MWgroup, data = MWtest, exact = FALSE, correct = FALSE)
wilcox.test(formula = abs_x ~ MWgroup, data = MWtest, exact = FALSE, correct = FALSE)

library(survey)
MWdesign <- svydesign(~0, data = MWtest)
svyranktest(abs_score ~ MWgroup, MWdesign, test = "vanderWaerden")

MWdesignweight <- svydesign(~0, data = MWtest, weights = ~abs_x)
svyranktest(abs_score ~ MWgroup, MWdesignweight)

# only far from zero coefficients
coefs_nz <- coefs %>% filter(abs(x) > 0.04) # also tried 0.05

testnz <- coefs_nz %>% left_join(., study_scores, by = c('group_id' = 'Entrez ID'))

cor.test(testnz$x, testnz$Score, use = 'complete.obs' )
plot(testnz$x, testnz$Score, use = 'complete.obs')

#####
# Tbls[[3]], only has certain genes classified as influencing a given phenotype,
# so let's try to compare how many of those hits are in the top coefficients in our model.
study3 <- tbls[[3]][grepl('cell death', tbls[[3]]$Phenotype)|grepl('proliferation', tbls[[3]]$Phenotype),]

most_neg_coefs <- coefs %>% arrange(x) %>% filter(row_number() <= 66)
most_pos_coefs <- coefs %>% arrange(desc(x)) %>% filter(row_number() <= 66)

study_hit_in_dat <- study3 %>% filter(`Entrez ID` %in% coefs$group_id)

TP_neg <- study_hit_in_dat %>% mutate(mostneg = study_hit_in_dat$`Entrez ID` %in% most_neg_coefs$group_id)
TP_pos <- study_hit_in_dat %>% mutate(mostpos = study_hit_in_dat$`Entrez ID` %in% most_pos_coefs$group_id)
# hits that are common are very very low, not really necessary to run the fishers exact test I know it will be bad...

# TODO: finish comparison
cont_tbl <- matrix(c(), nrow = 2, ncol = 2)

#####
# compare the two studies to see if they even found similar things

study_1 <- tbls[[2]] %>% select(`Entrez ID`, Score) %>% filter(!is.na(`Entrez ID`))
study_2 <- tbls[[3]] %>% select(`Entrez ID`, Phenotype) %>% filter(!is.na(`Entrez ID`))

studys <- full_join(study_1, study_2, by = c('Entrez ID')) %>%
  mutate(Neg_eff = grepl('cell death', Phenotype), 
         Pos_eff = grepl('proliferation', Phenotype)) %>%
  group_by(`Entrez ID`) %>%
  mutate(study2Res = sum(c(-1*Neg_eff, Pos_eff)))

ggplot(studys, aes(x = study2Res, y = Score)) + geom_point()

studys %>% group_by(study2Res) %>% summarise(mean(Score, na.rm = T))

