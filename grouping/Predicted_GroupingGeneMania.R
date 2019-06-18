library(tidyverse)
library(hashmap)
library(igraph)
library(ArgumentCheck)


GeneManiaGrouping <- function(Network, GeneManiaDistance = 1, 
                              AppliedAfterNamingConvention = F,
                              generateTables = T, FileSuffix = 'DefaultSuf',
                              import_threshold = 0.5){
  # GMEdge refers to edge in the GeneMania network, CLEdges is an edge in the 
  # Collinearity network.
  
  # Group genes together using known gene interactions i.e. Physical 
  # interaction and protein domain information from GeneMania database
  
  # Argument check:
  #####
  Check <- ArgumentCheck::newArgCheck()
  
  if (GeneManiaDistance < 0) 
    ArgumentCheck::addError(
      msg = "'Distance' must be >= 0",
      argcheck = Check
    )
  
  if (typeof(AppliedAfterNamingConvention) != 'logical')
    ArgumentCheck::addError(
      msg = "AfterNamingConbention must be a logical value",
      argcheck = Check
    )
  
  if (typeof(generateTables) != 'logical')
    ArgumentCheck::addError(
      msg = "generateTables must be a logical value",
      argcheck = Check
    )
  
  if (!(typeof(FileSuffix) %in% c('character', 'NULL'))) 
    ArgumentCheck::addError(
      msg = "'FilSuffix' must be a string or NULL",
      argcheck = Check
    )
  
  #* Return errors and warnings (if any)
  ArgumentCheck::finishArgCheck(Check)
  #####
  # all edges filtered for jaccard index above 0.5 also has row that labels if the
  # connection is kept due to the naming scheme or not, i.e. column `"isoform"?`  
  GeneNetwork <- Network#read_csv(paste('./grouping/Gene_Collinearity', import_threshold, 
                          #'.csv', sep = ''), col_types = cols())
  
  Gene_ids <- unique(c(GeneNetwork[['gene_id1']], GeneNetwork[['gene_id2']]))
  
  # Tested this and any gene in my data not in this mapping is not in the 
  # GeneMania database by entrez id. Found out later this is a problem of just 
  # missing the Entrez id and the gene is still in GeneMania only not listed 
  # under that identification, using secondary ensembl ids from NCBI to denote 
  # those that are missing
  GMNameMappings <- read_tsv('../../CellSummerProjectData/GeneMania/identifier_mappings.txt',
                              col_types = cols())
  Ensemble2Entrez <- GMNameMappings %>% filter(Source == 'Entrez Gene ID') %>%
    select(-Source) %>% filter(Name %in% Gene_ids)
  
  Physical_Interactions <- read_tsv('../../CellSummerProjectData/GeneMania/Physical_interactions.tsv', 
                                    col_types = cols())
  # Remove extra headers from concatenating the entire downloads with wget.
  Physical_Interactions <- filter(Physical_Interactions, Gene_A != "Gene_A") 
  
  Protein_Domains <- read_tsv('../../CellSummerProjectData/GeneMania/Protein_Domains.tsv', 
                              col_types = cols())
  Protein_Domains <- Protein_Domains %>% filter(Gene_A != 'Gene_A')
  
  GMEdges <- rbind(Physical_Interactions, Protein_Domains)
  
  NCBIEnsb2ID <- read_tsv('../../CellSummerProjectData/NCBI/gene2ensembl.tsv') %>%
    filter(GeneID %in% Gene_ids) %>% select(GeneID, Ensembl_gene_identifier) %>%
    distinct()
  
  NCBIEnsb2ID <- NCBIEnsb2ID %>% left_join(select(GMNameMappings, -Source), 
                                    by = c( 'Ensembl_gene_identifier' = 'Name')) %>% 
    mutate(GeneID = as.character(GeneID))
  
  GMEdges <- GMEdges %>% left_join(select(NCBIEnsb2ID, -Ensembl_gene_identifier), by = c('Gene_A' = 'Preferred_Name')) %>%
    rename('NCBIid_A' = 'GeneID')
  GMEdges <- GMEdges %>% left_join(select(NCBIEnsb2ID, -Ensembl_gene_identifier), by = c('Gene_B' = 'Preferred_Name')) %>%
    rename('NCBIid_B' = 'GeneID')
  
  GMEdges <- left_join(GMEdges, Ensemble2Entrez, by = c('Gene_A' = 'Preferred_Name')) %>%
    rename('GM_A' = 'Name')
  GMEdges <- left_join(GMEdges, Ensemble2Entrez, by = c('Gene_B' = 'Preferred_Name')) %>%
    rename('GM_B' = 'Name')
  
  GMEdges <- GMEdges %>% mutate(Name_A = if_else(!is.na(GM_A), GM_A, NCBIid_A),
                                Name_B = if_else(!is.na(GM_B), GM_B, NCBIid_B))
  
  #Ensembl Id's in the physical Interaction and Protein Domain connections but not
  # in the Ensembl to Entrez mapping provided by GeneMania.
  MissingEntrez_A <- GMEdges %>% filter(is.na(Name_A)) %>% 
    select(Gene_A, Name_A) %>% 
    rename(Gene = Gene_A, Name = Name_A) %>%
    distinct()
  
  MissingEntrez_B <- GMEdges %>% filter(is.na(Name_B)) %>% 
    select(Gene_B, Name_B) %>%
    rename(Gene = Gene_B, Name = Name_B) %>%
    distinct()
  
  MissingEntrez <- rbind(MissingEntrez_A, MissingEntrez_B) %>% distinct()
  
  # build graph from GMEdges
  GMEdges <- GMEdges %>% mutate(Name_1 = ifelse(Name_A < Name_B, Name_A, Name_B),
                            Name_2 = ifelse(Name_A < Name_B, Name_B, Name_A)) %>%
    select(Name_1, Name_2) %>%
    filter(!is.na(Name_1) & !is.na(Name_2)) %>%
    distinct()
  
  UniqueIds <- union(GMEdges$Name_1, GMEdges$Name_2)
  
  graph <- graph_from_data_frame(GMEdges, directed = F)
  
  neighbours <- ego(graph = graph, nodes = UniqueIds, order = 1)
  names(neighbours) <- UniqueIds
  
  
  # Plan is to find neighbours for a given gene of a specific order and see if the
  # target node in our network is in the neighbours of a given order in the 
  # genemania data.
  
  if (AppliedAfterNamingConvention) {
    UnclassifiedCLEdges <- filter(GeneNetwork, `"isoform"?` == F)
  } else {
    UnclassifiedCLEdges <- GeneNetwork
  }
  
  GMClassifiedCLEdges <- UnclassifiedCLEdges %>% rowwise() %>%
    mutate(gene_id1 = as.character(gene_id1), gene_id2 = as.character(gene_id2)) %>%
    mutate(neighbour = list(names(neighbours[[gene_id1]]))) %>%
    mutate(GMConnection = gene_id2 %in% neighbour)
  
  
  # Generate tables to import into cytoscape to see the networks.
  checkConnections <- GMClassifiedCLEdges %>% select(-neighbour)
  # this tables information is contained in the following written file All Collinearity, use that instead.
  #if(generateTables){
  #  write_csv(checkConnections, paste('./grouping/GeneManiaGroupedDistance_', 
  #                                    GeneManiaDistance, FileSuffix, '.csv', sep = ''))  
  #}
  
  
  tmp <- GeneNetwork %>% mutate(gene_id1 = as.character(gene_id1), gene_id2 = as.character(gene_id2))
  test <- full_join(GMClassifiedCLEdges, tmp) %>% 
    mutate(GMConnection = if_else(is.na(GMConnection), F, GMConnection)) %>%
    select(-neighbour)
  
  # test <- test %>% mutate(edgecolour = if_else(value == 1, 4, 
  #                                              if_else(`"isoform"?` & GMConnection, 3, 
  #                                                      if_else(GMConnection, 2, 
  #                                                              if_else(`"isoform"?`, 1, 0)))))
  
  if(generateTables){
    write_csv(test, paste('./Cytoscape/AllCollinearity', import_threshold,
                          'GMGrpdDistance_', GeneManiaDistance, FileSuffix, 
                          '.csv', sep = ''))  
  }
  # Removing column neighbours because it's a list which is tough to deal with 
  # and the information is no longer needed.
  GMClassifiedCLEdges <- GMClassifiedCLEdges %>% select(-neighbour)
  
  GMClassifiedCLEdges
}

#test <- GeneManiaGrouping(generateTables = T, FileSuffix = 'PreNameApplication_NoCC2_AvgReps', import_threshold = 0.5)
#test2 <- GeneManiaGrouping(generateTables = T, FileSuffix = 'PreNameApplication', import_threshold = 0.5)




