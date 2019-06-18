library(tidyverse)

# Calculate the correlation between the siRNA results and the magnitudes of the
#coefficients found by our linear regression.

# Read in the siRNA dataset from Nish
folder <- '../../CellSummerProjectData/ValidNish/Coulter/'
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

val <- rbind(val_p1, val_p2) %>% group_by(gene_sym) %>%
  summarise_at(vars(num_cells_48h:med_diam_72h), mean)


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

val <- left_join(val, grpsNA, by = 'gene_sym') %>%
  left_join(., grps1, by = 'gene_sym') %>%
  left_join(., grps3, by = 'gene_sym') %>%
  left_join(., grps10, by = 'gene_sym')

# read in the coefficients from the models
# Helper fun
get_coef_files <- function(directory){
  files <- list.files(directory)
  files <- files[grep('.csv', files)]
  files <- files[grep('old|test|Interp|.sh|CellSize.*pairs:T|PTBaseline', files, invert = T)]
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
val <- val %>% mutate_at(vars(starts_with('group')), funs(if_else(is.na(.), gene_sym, .)))
all_data <- val
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
                 substr(split_col[3], 1, 2),
                 toupper(substr(meth, 1, 1)), sep = '_')
  
  cor <- data %>%
    summarise_at(vars(contains(':')),
                 function(x){cor(x, data[[column]],
                                 use = 'complete.obs', method = meth)}) %>%
    gather("Model", !! names[[1]])
  
  pval <- data %>%
    summarise_at(vars(contains(':')),
                 function(x){cor.test(x, data[[column]], 't',
                                      method = meth)$p.value}) %>%
    gather("Model", !! names[[2]])
  
  full_join(cor, pval)
}


Num48hS <- Calc_Cor(all_data, 'num_cells_48h', 'spearman')
#Num72hS <- Calc_Cor(all_data, 'num_cells_72h', 'spearman')

Num48hP <- Calc_Cor(all_data, 'num_cells_48h', 'pearson')
#Num72hP <- Calc_Cor(all_data, 'num_cells_72h', 'pearson')

Size48hS <- Calc_Cor(all_data, 'med_diam_48h', 'spearman')
#Size72hS <- Calc_Cor(all_data, 'med_diam_72h', 'spearman')

Size48hP <- Calc_Cor(all_data, 'med_diam_48h', 'pearson')
#Size72hP <- Calc_Cor(all_data, 'med_diam_72h', 'pearson')

cors <- Num48hS %>% 
  #left_join(., Num72hS, by = "Model") %>%
  left_join(., Num48hP, by = "Model") %>%
  #left_join(., Num72hP, by = "Model") %>%
  left_join(., Size48hS, by = "Model") %>%
  #left_join(., Size72hS, by = "Model") %>%
  left_join(., Size48hP, by = "Model") #%>%
  #left_join(., Size72hP, by = "Model")

write_csv(cors, './ValidateNishsiRNA/Correlation_Table_test.csv')

# Produce scatter plots of the models with intercept forced to 0.
labeller_CN <- function(chars){
  correlations <- (filter(cors, Model %in% chars) %>% arrange(Model))
  m <- regexpr("concs.*l", chars)
  cons <- regmatches(chars, m)
  paste(substr(cons, 1, nchar(cons)-1), 'cor=',
                   round(correlations$cor_num_48_S, 3), 'pval=',
                   round(correlations$pval_num_48_S, 3), sep = ' ')
}

labeller_CS <- function(chars){
  correlations <- (filter(cors, Model %in% chars) %>% arrange(Model))
  m <- regexpr("concs.*l", chars)
  cons <- regmatches(chars, m)
  paste(substr(cons, 1, nchar(cons)-1), 'cor=',
        round(correlations$cor_med_48_S, 3), 'pval=',
        round(correlations$pval_med_48_S, 3), sep = ' ')
}

# create object To Be Plotted (TBP)
TBP <- all_data %>% gather('Model', 'coef', starts_with('Cell'))
# two plots used in presentation of best correlation
ggplot(TBP[grep('CellNum.*Inter:F.*cs:10.*1se', TBP$Model),],
       aes(x = log(num_cells_48h), y = coef, label = gene_sym)) +
  geom_point() + ggrepel::geom_text_repel(size = 4, vjust = -1) +
  ggtitle('Cell Proliferation, Validated Genes')
ggsave('./graphs/For Paper/PDFs/CellNumValidationCorrelationtest.pdf')
ggsave('./graphs/For Paper/PNGs/CellNumValidationCorrelationtest.png')

ggplot(TBP[grep('CellSize.*Inter:F.*cs:1,.*Min', TBP$Model),],
       aes(x = med_diam_48h, y = coef, label = gene_sym)) +
  geom_point() + ggrepel::geom_text_repel(size = 4, vjust = -1) +
  ggtitle('Cell Size, Validated Genes')
ggsave('./graphs/For Paper/PDFs/CellSizeValidationCorrelationtest.pdf')
ggsave('./graphs/For Paper/PNGs/CellSizeValidationCorrelationtest.png')

pdf('./ValidateNishsiRNA/correlation_plotstest.pdf')
ggplot(TBP[grep('CellNum.*Inter:F.*Min', TBP$Model),],
       aes(x = num_cells_48h, y = coef, label = gene_sym)) +
  facet_wrap(~ Model, nrow = 2, labeller = as_labeller(labeller_CN)) +
  geom_point() + ggtitle('CellNum 48h, Intercept=0, Min lambda') +
  geom_text(size = 2.5, vjust = -1)

ggplot(TBP[grep('CellNum.*Inter:F.*1se', TBP$Model),],
       aes(x = num_cells_48h, y = coef, label = gene_sym)) +
  facet_wrap(~ Model, nrow = 2, labeller = as_labeller(labeller_CN)) +
  geom_point() + ggtitle('CellNum 48h, Intercept=0, 1se lambda') +
  geom_text(size = 2.5, vjust = -1)

ggplot(TBP[grep('CellSize.*Inter:F.*Min', TBP$Model),],
       aes(x = med_diam_48h, y = coef, label = gene_sym)) +
  facet_wrap(~ Model, nrow = 2, labeller = as_labeller(labeller_CS)) +
  geom_point() + ggtitle('CellSize 48h, Intercept=0, Min lambda') +
  geom_text(size = 2.5, vjust = -1)

ggplot(TBP[grep('CellSize.*Inter:F.*1se', TBP$Model),],
       aes(x = med_diam_48h, y = coef, label = gene_sym)) +
  facet_wrap(~ Model, nrow = 2, labeller = as_labeller(labeller_CS)) +
  geom_point() + ggtitle('CellSize 48h, Intercept=0, 1se lambda') +
  geom_text(size = 2.5, vjust = -1)
dev.off()






