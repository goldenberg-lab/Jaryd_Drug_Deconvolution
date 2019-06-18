library(tidyverse)
# Read in validation data and analyze the controls to see if the two populations
# differ significantly.

folder <- '../../CellSummerProjectData/ValidNish/'
val_p1 <- read_csv(paste0(folder, '20180801_siRNA-reseed_SL.csv'), 
                   col_names = c('exp_num', 'well_group', 'gene_sym', 
                                 'num_cells_48h', 'num_cells_72h', 
                                 'med_diam_48h', 'med_diam_72h'), 
                   col_types = 'ccciidd', skip = 3) %>%
  mutate(exp_num = 1) %>% filter(!is.na(gene_sym))

val_p2 <- read_csv(paste0(folder, '20180809_siRNA-reseed_SL.csv'), 
                   col_names = c('exp_num', 'well_group', 'gene_sym', 
                                 'num_cells_48h', 'num_cells_72h', 
                                 'med_diam_48h', 'med_diam_72h'), 
                   col_types = 'ccciidd', skip = 3) %>%
  mutate(exp_num = 2) %>% filter(!is.na(gene_sym))

val <- rbind(val_p1, val_p2)

# the experiments are ordered in the original file so that the first control in 
# each group is the first control in the table, so I will number them to make it
#  explicit which well corresponds to which well.
ctrls <- val %>% filter(gene_sym == 'siCON') %>% mutate(well_group = c(1,2,3,4,5,6,1,2,3,4,5,6))

# likely there is a mean shift between the two experiments in regards to the controls
t.test(ctrls$num_cells_48h[1:6], ctrls$num_cells_48h[7:12], paired = T) # pval: 0.2405
t.test(ctrls$num_cells_72h[1:6], ctrls$num_cells_72h[7:12], paired = T) # pval: 0.07275
t.test(ctrls$med_diam_48h[1:6], ctrls$med_diam_48h[7:12], paired = T) # pval: 0.0001216
t.test(ctrls$med_diam_72h[1:6], ctrls$med_diam_72h[7:12], paired = T) #pval: 0.1128

# compare the correlation of the difference between the two time points to 
# ascertain if there is evidence the same process was ongoing and there is 
# simply a shift in mean between the two experiments. Using KS.test, cor.test, 
# and t.test on the differences between the two time points.
val <- val %>% mutate(num_diff = num_cells_72h - num_cells_48h, 
                      diam_diff = med_diam_72h - med_diam_48h,
                      num_fold = num_cells_72h/num_cells_48h, 
                      diam_fold = med_diam_72h/med_diam_48h)

# boxplot shows there is a mean shift, but the values across experiments look 
# quite correlated on the scatter plots, so I think the best course of action is
# to use the mean between the two experiments as the single value, when 
# correlating with the coefficients.
ggplot(val) + geom_boxplot(aes(group = exp_num, y = num_cells_48h))
plot(x = val$num_cells_48h[1:31], y = val$num_cells_48h[32:62])

ggplot(val) + geom_boxplot(aes(group = exp_num, y = num_cells_72h))
plot(x = val$num_cells_72h[1:31], y = val$num_cells_72h[32:62])

ggplot(val) + geom_boxplot(aes(group = exp_num, y = med_diam_48h))
plot(x = val$med_diam_48h[1:31], y = val$med_diam_48h[32:62])

ggplot(val) + geom_boxplot(aes(group = exp_num, y = med_diam_72h))
plot(x = val$med_diam_72h[1:31], y = val$med_diam_72h[32:62])


# Looking at the differences between the two time points in each experiment, so 
# it carries forward that the base values are correlated so are the differences 
# so I believe the general trajectory of the cells growth was similar and there 
# was simply a reasonable shift in mean between the two experiments.
ks.test(val$num_diff[1:31], val$num_diff[32:62], 't')
ks.test(val$diam_diff[1:31], val$diam_diff[32:62], 't')

plot(x = val$num_diff[1:31], y = val$num_diff[32:62])
plot(x = val$diam_diff[1:31], y = val$diam_diff[32:62])

t.test(val$num_diff[1:31], val$num_diff[32:62], paired = T)
t.test(val$diam_diff[1:31], val$diam_diff[32:62], paired = T)

# Given the correlation was significantly not zero and the KS test failed to 
# reject the null hypothesis, I believe the differences in the original values 
# are primarily differences in mean so taking the average between the two 
# experiments seems reasonable.