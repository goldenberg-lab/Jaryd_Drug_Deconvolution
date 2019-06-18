# Create Plot of num_targets per compound, pre, mid, and post data cleaning
library(tidyverse)
source('./import_data.R')


PlateBookOrig <- read_csv('../../CellSummerProjectData/OriginalData/PlateBookID2.txt') %>% 
  select(inchi_key, gene_id, gene_sym) %>% distinct()
PlateBookOrthos <- read_csv('./Changed_Data/PlateBookID_Orthologs.csv') %>%
  select(inchi_key, ortho_id, ortho_sym) %>% distinct()
lst <- import_grpd(average_replicates = T, rmv_cell_count_2 = F, 
                    create_group_files = F, two_conc = c(10,20),
                    min_num_drugs = 4, add_pair_features = T,
                    file_suffix = paste0(paste0(c(10,20), collapse = ','), 'deleteme'))
JoinedTable <- lst[[1]] %>% select(inchi_key, group_id, group_sym) %>% distinct()

PlateBookOrigCount <- PlateBookOrig %>% count(inchi_key) %>% rename(`#targets` = n)
PlateBookOrthosCount <- PlateBookOrthos %>% count(inchi_key) %>% rename(`#targets` = n)
JoinedTableCount <- JoinedTable %>% count(inchi_key) %>% rename(`#targets` = n)

pdf('./graphs/DrugTargetDistributions.pdf')
ggplot(PlateBookOrigCount, aes(x = `#targets`)) + geom_histogram(binwidth = 1) + ggtitle('Original Distribution')
ggplot(PlateBookOrthosCount, aes(x = `#targets`)) + geom_histogram(binwidth = 1) + ggtitle('Transformed to Human Orthologs Distribution')
ggplot(JoinedTableCount, aes(x = `#targets`)) + geom_histogram(binwidth = 1) + ggtitle('Post Grouping and filter out Genes with Fewer than 4 compounds targeting it')
dev.off()
