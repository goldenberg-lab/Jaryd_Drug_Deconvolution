# rewritten for cleanliness and legibility, and to remove unuseful sections.

library(tidyverse)

PlateBookID <- read_csv('../../CellSummerProjectData/OriginalData/Platebook_with_Activity.csv')
max_activities <- max(str_count(PlateBookID$activity_moa, ';') + 1)
PlateBookID <- PlateBookID %>% separate(activity_moa, paste0('activity_moa', 1:max_activities))
# build counts of types of activity on genes.
# code assumes activity_moa and PlateBookID have been joined previously.

all_activities <- gather(PlateBookID) %>% .$value %>% unique()

activity_counts <- PlateBookID %>% group_by_at(vars(starts_with('activity'))) %>%
  summarise(count = n())

write_csv(activity_counts, './moa_activities/AllActivities.csv')
