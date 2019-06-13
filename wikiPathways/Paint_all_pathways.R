library(rWikiPathways)
library(tidyverse)
library(RCy3)
library(RColorBrewer)

# General notes to self:
# - Need version 3.3.1 of wikipathways app in cytoscape newer versions fail to 
#   load gpml files.
# - Cytoscape must be running this script will check that.
# - all directories are hardcoded, future project to fix that if helpful.
# - assumes Human pathways already downloaded from wikipathways, code is 
#   commented if needed later.
# 

# Check cytoscape is online
stopifnot(cytoscapePing() == c("You are connected to Cytoscape!"))

# if need to download pathways from wikipathways
#downloadPathwayArchive(organism="Homo sapiens", format='gpml', destpath = './wikiPathways/HumanTest/')

dir <- '/data/Jaryd/R/CellSummerProject/wikiPathways/Human/'

Pathways <- list.files('./wikiPathways/Human/') %>%
  .[grep('\\.gpml', .)] %>%
  paste0(dir, .)

# Read in database and ID info from pathway files.
source("wikiPathways/GPML_Helpers.R")
all_nodes <- lapply(Pathways, function(path){get_Nodes(path) %>% mutate(file = path)})

all_nodes_Tbl <- Reduce(rbind, all_nodes) %>%
  mutate(pathway = basename(file))

NCBI_Info <- read_tsv('../../CellSummerProjectData/NCBI/Homo_sapiens.gene_info') 


#TODO: 
# 1. There seems to be genes which when mapped with this NCBI synonyms data 
# there are still those that are missing. My idea is to then incorporate a mapping
# from the Ensembl id's for some of the nodes (those with Ensembl listed as the
# database) to the Entrez_ids. I think I might have the file already in my NCBI 
# folder. Should try this if synonyms is not enough.
# 
# 2. Some network's images look like the image was saved prior to completion of 
# the painting process. I think it could be more complicated networks are sending 
# a complete signal back to the script before it is fully finished (could try 
# adding a wait to solve this in that case)? My other theory about why some 
# failed was that since it is running through the gui my using the computer at 
# the same time caused some interruptions so it messed up the execution.
#  sol'n attempts: 1. tried adding sys.sleep(0.5) between each call, didn't affect 
#                  rate of problems on sample size of 13.
#  sol'n: loading the network twice if it is the fist network loaded after 
#  clearing fixes the issue in most cases there are still some instances where 
#  the background shapes tend to be misplaced but they are few.
#                  
max_Syns <- NCBI_Info$Synonyms %>% str_count('\\|') %>% max(.)

NCBI_Syns <- NCBI_Info %>%
  select(Symbol, Synonyms) %>%
  distinct() %>%
  separate(Synonyms, paste0('Synonym', 1:max_Syns), sep = '\\|', 
           extra = 'merge', fill = 'right') %>%
  mutate(Synonym0 = Symbol) %>%
  gather("Syn_num", "Synonym", -Symbol) %>%
  filter(!is.na(Synonym))

grps <- read_csv('./grouping/latestgroups10,20CellNum.csv') %>% 
  select(group_id, group_sym, ortho_id, ortho_sym) %>% 
  mutate(group_id = if_else(is.na(group_id), ortho_id, group_id), 
         group_sym = if_else(is.na(group_sym), ortho_sym, group_sym)) %>%
  distinct()

CNSingleGene <- read_csv(paste0('./Coefs/Penalty/CellNum/CellNumMinDrgs:4_alpha:0_',
                      'pairs:TRUE_Inter:FALSE_avgreps:TRUE_concs:10,20lambda1se',
                      '_rmvdups2019-03-25_Interp.csv')) %>% filter(group_id == group_id2, isinterp1020) %>%
  select(-group_id2, -group_sym2, -genekey, -isinterp1020) %>%
  left_join(., grps, by = c('group_id', 'group_sym'))

Node_w_Coefs <- left_join(all_nodes_Tbl, NCBI_Syns, by = c('name' = 'Synonym')) %>%
  mutate(Symbol = if_else(database == 'Entrez Gene' & is.na(Symbol), name, Symbol)) %>%
  left_join(., CNSingleGene, by = c('Symbol' = 'ortho_sym'))

# The missing values in the x column prevents the colour mapping from being 
# continuos. Store missing in separate column.
Name2coef <- Node_w_Coefs %>% select(name, x, group_sym) %>% as.data.frame(.) %>%
  mutate(missing = if_else(is.na(x), T, F), x = if_else(is.na(x), 0, x),
         group_sym = if_else(is.na(group_sym), 'def', group_sym)) %>%
  distinct()

#Load, and "Paint" a pathway from a gpml with the relevant coefficients.
Paint_Pathway <- function(gpmlfile, out_folder = '/data/Jaryd/R/CellSummerProject/wikiPathways/Cytoscape_Out/',
                          group_colours = c('#000000')){
  if(length(getCollectionList()) %% 10 == 0){
    deleteAllNetworks()
    # There is a high probability of the network being loaded with errors if it 
    # is the first network after deleting all others hence the double load.
    try(importNetworkFromFile(gpmlfile)) 
  }
  try(importNetworkFromFile(gpmlfile))
  loadTableData(Name2coef, data.key.column = 'name',
                table.key.column = 'name')
  setVisualStyle('default')
  setNodeColorMapping('x', c(min(Node_w_Coefs$x, na.rm = T),0,max(Node_w_Coefs$x, na.rm = T)), 
                      c('#FF0000', '#AAAAAA', '#0000FF'), 
                      default.color = '#FFFFFF', mapping.type = 'c')
  #setNodeBorderColorDefault('#000000')
  ### added these lines
  # TODO: finish changing the border colour to correspond with which group the gene belongs to.
  # Create a mapping of a "random" colour for each group with more than one member.
  
  browser()
  grouptocolour <- Node_w_Coefs %>% filter(pathway == basename(gpmlfile) & !is.na(x)) %>%
    select(Symbol, group_sym) %>% distinct() %>% group_by(group_sym) %>% add_count() %>%
    filter(n > 1) %>% select(group_sym) %>% distinct() %>% 
    ungroup() %>%
    mutate(colour = sample(c(brewer.pal(9, "Set1"), brewer.pal(12,'Paired'))[c(-1,-2,-5,-7,-11,-15)], 
                           length(.[[1]])))
  
  if(length(grouptocolour[[1]] > 0)){
    setNodeBorderColorMapping('group_sym', 
                              table.column.values = grouptocolour$group_sym, 
                              colors = grouptocolour$colour, 
                              mapping.type = 'd', default.color = '#000000')  
    setNodeBorderWidthMapping('missing', c(T,F), c(0.5,3), mapping.type = 'd')
  }
  
  ### 
  # setNodeBoderWidthDefault(0.5)
  missing <- selectNodes(T, by.col = 'missing')
  setNodeColorBypass(missing$nodes, new.colors = c('#FFFFFF'))
  clearSelection()
  fitContent()
  exportImage(paste0(out_folder, basename(gpmlfile)), type = 'SVG')
}
##### Interactive code which helped me find bugs and make choices.
### Figure out which colours to use as a palette to draw randomly from for groups.
# source('../Helpers/Colour_helper.R')
# pot_cols <- c(brewer.pal(9, "Set1"), brewer.pal(12,'Paired'))
# disp_cols(pot_cols)
# 
# rmvd_cols <- pot_cols[c(-1,-2,-5,-7,-11,-15)]
# disp_cols(rmvd_cols)
# # Leave me with exactly 15 colours, which is the largest number of groups in a 
# # single pathway so far.
###

###
# Figure out if there are more genes which are both a gene and a synonym for a different gene.
All_Synonyms <- NCBI_Info %>% filter(Synonyms != '-')
tmp <- paste0(All_Synonyms$Synonyms, collapse = '|')
All_Syms_which_are_Syns <- NCBI_Info$Symbol[sapply(NCBI_Info$Symbol, function(Sym){grepl(paste0('\\|', Sym, '\\|'), tmp)})]

MyGenes_which_are_Syns <- CNSingleGene$ortho_sym[CNSingleGene$ortho_sym %in% NCBI_Info$Symbol[All_Syms_which_are_Syns]]

###
#####
dir.create('/data/Jaryd/R/CellSummerProject/wikiPathways/Cytoscape_Out/SessionFiles/', showWarnings = F)
# file.remove(dir('/data/Jaryd/R/CellSummerProject/wikiPathways/Out_test/',
#                 pattern = "\\.svg$", full.names = TRUE))

### for testing
# Pathways <- list.files('./wikiPathways/HumanTest/') %>%
#   .[grep('\\.gpml', .)] %>%
#   paste0(dir, .)
###

results <- lapply(Pathways, 
                  function(x, out_dir){
                    if(file.exists(paste0(out_dir, sub('gpml$', 'svg', basename(x))))){
                      return("File already exists")
                    } else {
                      Paint_Pathway(x, out_folder = out_dir)
                    }
                  }, out_dir = '/data/Jaryd/R/CellSummerProject/wikiPathways/Cytoscape_Out/')
#lapply(Pathways[219:length(Pathways)], Paint_Pathway)

results

#TODO: split the groups into a row per member, and figure out a way to visually 
# represent multiple members of the same group on the pathway. Run this script 
# over the weekend and save these ones in a separate folder.
# Meeting with Ran week of June 10th-14th, talk to him about adding a visual for
# the groups then. 
# My first idea is to add a colour to the border so that all same coloured 
# borders are members of the same group?






