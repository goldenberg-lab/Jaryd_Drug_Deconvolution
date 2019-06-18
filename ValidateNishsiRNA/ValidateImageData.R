library(tidyverse)

# Looking at the image processed validation data to confirm similar results.
dir <- '../../CellSummerProjectData/ValidNish/HeLa_image/'

Well_gene_map <- read_csv(paste0(dir, 'Well2Gene.csv'))

i <- 0
HeLa_image <- list.files(dir) %>% 
  .[grep('statTab.*csv', .)] %>% 
  paste0(dir, .) %>%
  lapply(., read_csv) %>%
  lapply(., function(x){assign('i', get('i', globalenv()) + 1, globalenv()); mutate(x, Plate = i)}) %>%
  bind_rows() %>%
  full_join(., Well_gene_map, by = c('Plate', 'Row', 'Column')) %>%
  .[c('Plate', 'Row', 'Column', 'Label', 'gene_sym', 
      setdiff(names(.), c('Plate', 'Row', 'Column', 'Label', 'gene_sym')))] %>%
  group_by(gene_sym) %>%
  summarise_at(vars(CellNum_1:sizeMAD_7), mean)

# read in groups used for model fit and use the rep coefs for these genes. Add 
# the group_sym for the genes for the different concentrations... Then 
# recalculate the correlations and don't forget to calculate the pvals with 
# those and report together.
grpsNA <- read_csv('./grouping/latestgroupsNA.csv') %>% 
  select(gene_sym, group_sym) %>% dplyr::rename(group_symNA = group_sym)
grps1 <- read_csv('./grouping/latestgroups1,2.csv') %>% 
  select(gene_sym, group_sym) %>% dplyr::rename(group_sym1 = group_sym)
grps3 <- read_csv('./grouping/latestgroups3,6.csv') %>% 
  select(gene_sym, group_sym) %>% dplyr::rename(group_sym3 = group_sym)
grps10 <- read_csv('./grouping/latestgroups10,20.csv') %>% 
  select(gene_sym, group_sym) %>% dplyr::rename(group_sym10 = group_sym)

HeLa_image <- left_join(HeLa_image, grpsNA, by = 'gene_sym') %>% 
  left_join(., grps1, by = 'gene_sym') %>% 
  left_join(., grps3, by = 'gene_sym') %>%
  left_join(., grps10, by = 'gene_sym')

get_coef_files <- function(directory){
  files <- list.files(directory)
  files <- files[grep('old|test|Interp|.sh|CellSize.*pairs:T|Baseline|NoPairs', files, invert = T)]
  files <- paste0(directory, files)
  files
}

CellNum_files <- get_coef_files('./Coefs/Penalty/CellNum/')
CellSize_files <- get_coef_files('./Coefs/Penalty/CellSize/')

read_coefs <- function(files, split_by = '', keep = 1, suffix = ''){
  
  new_x <- str_split(files[[1]], split_by)[[1]][keep] %>% str_flatten(.) %>% paste0(suffix, .)
  coefs <- read_csv(files[[1]]) %>% select(-genekey)
  coefs <- cbind(coefs[-1], coefs[1]) # reorder columns
  
  # rename coefficient column to represent the file in it's name
  names(coefs) <- c(names(coefs)[-length(names(coefs))], new_x)
  
  for (file in files[-1]){
    tmp <- read_csv(file) %>% select(-genekey)
    new_x <- str_split(file, split_by)[[1]][keep] %>%  str_flatten(.) %>% paste0(suffix, .)
    names(tmp) <- c(new_x, names(tmp)[-1])
    coefs <- full_join(coefs, tmp)
  }
  
  coefs
}

CellNum_coefs <- read_coefs(CellNum_files, '_', c(4,6), suffix = 'CellNum_')

CellSize_coefs <- read_coefs(CellSize_files, '_', c(5,7), suffix = 'CellSize_')

# fill in the missing value if the gene is it's own group representative, then 
# loop over different concs and join by it's respective rep genes per concentration.
HeLa_image <- HeLa_image %>% mutate_at(vars(starts_with('group')), funs(if_else(is.na(.), gene_sym, .)))
all_data <- HeLa_image
for(con in c('NA', '1,', '3,', '10,')){
  sym <- paste0('group_sym', gsub('[[:punct:] ]+','',con))
  by <- set_names('group_sym', sym)
  all_data <- left_join(all_data, CellNum_coefs %>% 
                          filter(group_sym == group_sym2) %>% 
                          select(group_sym, contains(con)), 
                        by = by) %>% 
    left_join(., CellSize_coefs %>% 
                select(group_sym, contains(con)), by = by)
  print(con)
}

Calc_Cor <- function(data, column, meth){
  split_col <- str_split(column, '_')[[1]]
  names <- paste(c("cor", "pval"),
                 split_col[1],
                 toupper(substr(meth, 1, 1)), sep = '_')
  
  cor <- data %>% 
    summarise_at(vars(contains(':')), 
                 function(x){cor(x, all_data[[column]], 
                                 use = 'complete.obs', method = meth)}) %>% 
    gather("Model", !! names[[1]])
  
  pval <- data %>%
    summarise_at(vars(contains(':')), 
                 function(x){cor.test(x, all_data[[column]], 't', 
                                      method = meth)$p.value}) %>%
    gather("Model", !! names[[2]])
  
  full_join(cor, pval)
}

CellNum_S <- Calc_Cor(all_data, 'CellNum_2', 'spearman')
#CellNum_1S <- Calc_Cor(all_data, 'CellNum_1', 'spearman')
CellSize_S <- Calc_Cor(all_data, 'MedianSize_1', 'spearman')


CellNum_P <- Calc_Cor(all_data, 'CellNum_2', 'pearson')
#CellNum_1P <- Calc_Cor(all_data, 'CellNum_1', 'pearson')
CellSize_P <- Calc_Cor(all_data, 'MedianSize_1', 'pearson')

cors <- left_join(CellNum_S, CellNum_P, by = "Model") %>% 
  left_join(., CellSize_S, by = "Model") %>%
  left_join(., CellSize_P, by = "Model")

write_csv(cors, './ValidateNishsiRNA/Correlation_Table_HeLa_image.csv')

# Produce scatter plots of the models with intercept forced to 0.
labeller_CN <- function(chars){
  correlations <- (filter(cors, Model %in% chars) %>% arrange(Model))
  m <- regexpr("concs.*l", chars)
  cons <- regmatches(chars, m)
  paste(substr(cons, 1, nchar(cons)-1), 'cor=',
        round(correlations$cor_CellNum_S, 3), 'pval=',
        round(correlations$pval_CellNum_S, 3), sep = ' ')
}

labeller_CS <- function(chars){
  correlations <- (filter(cors, Model %in% chars) %>% arrange(Model))
  m <- regexpr("concs.*l", chars)
  cons <- regmatches(chars, m)
  paste(substr(cons, 1, nchar(cons)-1), 'cor=',
        round(correlations$cor_MedianSize_S, 3), 'pval=',
        round(correlations$pval_MedianSize_P, 3), sep = ' ')
}

# create object To Be Plotted (TBP)
TBP <- all_data %>% gather('Model', 'coef', contains(':'))

pdf('./ValidateNishsiRNA/correlation_plots_HeLa_Image.pdf')
ggplot(TBP[grep('CellNum.*Inter:F.*Min', TBP$Model),], 
       aes(x = CellNum_1, y = coef, label = gene_sym)) + 
  facet_wrap(~ Model, nrow = 2, labeller = as_labeller(labeller_CN)) + 
  geom_point() + ggtitle('CellNum 48h, Intercept=0, Min lambda') +
  geom_text(size = 2.5, vjust = -1)

ggplot(TBP[grep('CellNum.*Inter:F.*1se', TBP$Model),], 
       aes(x = CellNum_1, y = coef, label = gene_sym)) + 
  facet_wrap(~ Model, nrow = 2, labeller = as_labeller(labeller_CN)) + 
  geom_point() + ggtitle('CellNum 48h, Intercept=0, 1se lambda') +
  geom_text(size = 2.5, vjust = -1)

ggplot(TBP[grep('CellSize.*Inter:F.*Min', TBP$Model),], 
       aes(x = MedianSize_1, y = coef, label = gene_sym)) + 
  facet_wrap(~ Model, nrow = 2, labeller = as_labeller(labeller_CS)) + 
  geom_point() + ggtitle('CellSize 48h, Intercept=0, Min lambda') +
  geom_text(size = 2.5, vjust = -1)

ggplot(TBP[grep('CellSize.*Inter:F.*1se', TBP$Model),], 
       aes(x = MedianSize_1, y = coef, label = gene_sym)) + 
  facet_wrap(~ Model, nrow = 2, labeller = as_labeller(labeller_CS)) + 
  geom_point() + ggtitle('CellSize 48h, Intercept=0, 1se lambda') +
  geom_text(size = 2.5, vjust = -1)
dev.off()

