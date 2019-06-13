library(tidyverse)

# Chosen isoform function
if_isoform3 <- function(genesym_1, genesym_2) {
  start1 <- strsplit(genesym_1, '[0123456789]')[[1]][1]
  start2 <- strsplit(genesym_2, '[0123456789]')[[1]][1]
  
  min_len <- max(min(str_length(start1), str_length(start2)), 2)
  
  if (tolower(substr(start1, 1, min_len)) == tolower(substr(start2, 1, min_len))){
    return(T)
  } else if(tolower(substr(genesym_1, 1, 3)) == tolower(substr(genesym_2, 1, 3))){
    return(T)
  } else{
    return(F)
  }
}

# function to label if there is a connection in GeneMania or not.

#PlateBook must include column num_drugs. Function called in import_data script.
AddGroups <- function(PlateBook, threshold = 0.5, generate_files_short = F, Concentration = NA,
                      file_suffix = '', Rmv_CC2 = F, generate_files_long = F){
  if(generate_files_long){
    generate_files_short <- T
  }

  if (is.na(Concentration[[1]])){
    Allgenes <- read_csv('./Linearity/AllLinearityConc:NA.csv')
  } else if(length(Concentration)== 2){
    if(1 %in% Concentration & 2 %in% Concentration){
      Allgenes <- read_csv('./Linearity/AllLinearityConc:1,2.csv')
    } else if (3 %in% Concentration & 6 %in% Concentration){
      Allgenes <- read_csv('./Linearity/AllLinearityConc:3,6.csv')
    } else if (10 %in% Concentration & 20 %in% Concentration){
      Allgenes <- read_csv('./Linearity/AllLinearityConc:10,20.csv')
    } else {
      stop('AllLinearity has not been calculated for this pair of concentrations.')
    } 
  } else{
    stop('AllLinearity has not been calculated for this set of concentrations.')
  }
  
  # Remove any genes in Allgenes (the colinearity table) which have been removed 
  # elsewhere in the pipeline post colinearity calculation. 
  # ex: the pseudogenes removed in import_data.
  GenesinPB <- (PlateBook %>% mutate(gene = paste(ortho_sym, '(', ortho_id, ')', sep = ' ')))$gene
  
  AllgenesGrtr0.5 <- filter(Allgenes, value > threshold & Var1 %in% GenesinPB & Var2 %in% GenesinPB)
  
  if (generate_files_short){
    write_csv(AllgenesGrtr0.5, paste0('./Linearity/AllLinearityAbove0.5', file_suffix,'.csv'))
  }
  
  # Add counts
  DrugCounts1 <- PlateBook %>% mutate(Var1 = paste(ortho_sym, '(', ortho_id, ')'),
                                        Var2 = paste(ortho_sym, '(', ortho_id, ')')) %>%
    rename(ortho_id1 = ortho_id) %>%
    ungroup() %>% select(Var1, ortho_id1, num_drugs) %>% distinct()
  
  DrugCounts2 <- PlateBook %>% mutate(Var1 = paste(ortho_sym, '(', ortho_id, ')'),
                                        Var2 = paste(ortho_sym, '(', ortho_id, ')')) %>%
    rename(ortho_id2 = ortho_id) %>%
    ungroup() %>% select(Var2, ortho_id2, num_drugs) %>% distinct()
  
  tmp <- full_join(AllgenesGrtr0.5, DrugCounts1 , by = c('Var1'), suffix = c('1', '2')) %>% filter(!is.na(Var2))
  AllgenesGrtr0.5 <- full_join(tmp, DrugCounts2 , by = c('Var2'), suffix = c('1', '2')) %>% filter(!is.na(Var1))
  
  # Check if the two genes are Isoforms, then select the gene with the most 
  # drugs targeting to be the representative of the two. Essentially, building 
  # the edges of the graph which will be traversed later to construct final groupings.
  AllgenesGrtr0.5 <- AllgenesGrtr0.5 %>% rowwise() %>% mutate(`"isoform"?` = if_isoform3(Var1, Var2))
  
  # browser()
  source('GroupingGeneMania.R')
  AllgenesGrtr0.5 <- GeneManiaGrouping(AllgenesGrtr0.5, generateTables = generate_files_short,
                                       FileSuffix = file_suffix)
  if(generate_files_long){
    tmp <- full_join(Allgenes, DrugCounts1 , by = c('Var1'), suffix = c('1', '2')) %>% filter(!is.na(Var2))
    Allgenes <- full_join(tmp, DrugCounts2 , by = c('Var2'), suffix = c('1', '2')) %>% filter(!is.na(Var1))
    AllgenesPresent <- Allgenes %>% filter(value > 0)
    AllgenesPresent <- GeneManiaGrouping(AllgenesPresent, generateTables = generate_files_long,
                                  AppliedAfterNamingConvention = F)
  }
  
  # Saving csv to import in cytoscape to see unmolested network of drug 
  # collinearities. Use function at bottom to save this table.
  UnFilteredNetwork <- AllgenesGrtr0.5
  
  # Old code when using naming scheme to group drugs instead of GeneMania data.
  # changes <- AllgenesGrtr0.5 %>% rowwise() %>% mutate(merge = if_isoform3(Var1, Var2)) %>%
  #   filter(merge == T) %>% mutate(old = if_else(num_drugs1 > num_drugs2, Var2, 
  #                                               ifelse(num_drugs1 == num_drugs2, 
  #                                                      ifelse(ortho_id1 < ortho_id2, Var2, Var1), Var1)),
  #                                 new = if_else(num_drugs1 > num_drugs2, Var1, 
  #                                               ifelse(num_drugs1 == num_drugs2, 
  #                                                      ifelse(ortho_id1 < ortho_id2, Var1, Var2), Var2))) %>% select(old, new)
  # 
  # Sets Representative gene for each connection i.e. gives directionality to 
  # the edges.
  
  changes <- AllgenesGrtr0.5 %>% filter(GMConnection == T | value == 1) %>% 
    mutate(old = if_else(num_drugs1 > num_drugs2, Var2, 
                                                if_else(num_drugs1 == num_drugs2, 
                                                       if_else(ortho_id1 < ortho_id2, Var2, Var1), Var1)),
                                  new = if_else(num_drugs1 > num_drugs2, Var1, 
                                                if_else(num_drugs1 == num_drugs2, 
                                                       if_else(ortho_id1 < ortho_id2, Var1, Var2), Var2))) %>% select(old, new)
  
  # Depth first search to set value to be the best name for the isoform group.
  lst <- changes$new
  names(lst) <- changes$old
  original <- lst # is for testing purposes.
  for (i in 1:length(lst)){
    gene <- lst[i]
    while (!is.na(lst[gene])){
      lst[i] <- lst[gene]
      gene <- lst[i]
    }
  }
  
  visited <- rep(0, length(unique(c(names(lst), lst))))
  names(visited) <- unique(c(names(lst), lst))
  groups <- vector("list", length(unique(c(names(lst), lst))))
  names(groups) <- unique(c(names(lst), lst))
  backedge <- names(lst)
  names(backedge) <- lst
  
  # traverse the network labeling each bunch of connected nodes as a group.
  groups <- Label_subnets(list_of_edges = lst)
  
  # function from online Hadley Wickham post, change list to df
  list_to_df <- function(listfordf){
    if(!is.list(listfordf)) stop("it should be a list")
       
    df <- list(list.element = listfordf)
    class(df) <- c("tbl_df", "data.frame")
    attr(df, "row.names") <- .set_row_names(length(listfordf))
       
    if (!is.null(names(listfordf))) {
      df$name <- names(listfordf)
    }
  
    df
  }
  
  groups_df <- list_to_df(groups) %>% rename(group_num = list.element) %>% 
    mutate(group_num = unlist(group_num)) %>% separate(name, c('ortho_sym', 'ortho_id', 'extra'), sep = '[ ][()]') %>%
    mutate(ortho_id = as.integer(ortho_id))
  
  groups_df <- full_join(groups_df, select(PlateBook, ortho_sym, ortho_id, num_drugs) %>% distinct(), 
                         by = c('ortho_sym', 'ortho_id'))
  
  group_sym <- list()
  group_id <- list()
  for (i in 1:max(groups_df$group_num, na.rm = T)){
    tmp <- filter(groups_df, group_num == i) %>% arrange(ortho_id)
    row <- which.max(tmp$num_drugs)
    if (!is_empty(row)){
      group_id[i] <- tmp$ortho_id[[row]]
      group_sym[i] <- tmp$ortho_sym[[row]]
    }
  }
  
  #function copied from a Jon M on stack exchange
  nullToNA <- function(x) {
    x[sapply(x, is.null)] <- NA
    return(x)
  }
  
  group_sym <- nullToNA(group_sym)
  group_id <- nullToNA(group_id)
  group_sym[311] <- NA
  group_id[311] <- NA
  
  groups_df <- groups_df %>% mutate(group_sym = as.character(group_sym[group_num]),  
                                    group_id = as.integer(as.character(group_id[group_num]))) %>% 
    select(group_num, group_sym, group_id, ortho_sym, ortho_id, num_drugs)
  
  # Fix how as.character set some NA's to "NA"
  groups_df$group_sym[groups_df$group_sym == 'NULL'] <- NA 
  groups_df$group_sym[groups_df$group_sym == 'NA'] <- NA
  
  # Write groups_df to csv after fixing multiple leaf problem. Shows mapping 
  # from gene through ortholog till group representative. 
  check <- select(PlateBook, gene_id, gene_sym, ortho_id, ortho_sym) %>% 
    full_join(., groups_df, by = c('ortho_id', 'ortho_sym')) %>%
    filter(!is.na(gene_id)) %>%
    distinct()
  
  if(generate_files_short){
    write_csv(check, paste0('./grouping/latestgroups',file_suffix,'.csv'))
  }
  
  #used to store group names for later merging to produce summary table of edges.
  group_names <- check %>% select(-ortho_id, -ortho_sym, -num_drugs, -gene_sym, -group_num)

  # export new groups for cytoscape visualization.
  tmp <- groups_df %>% filter(!is.na(group_num)) %>% 
    mutate(CytoGroup = paste0(paste(group_sym, group_id, sep = " ( "), ' )'),
           CytoOrtho = paste0(paste(ortho_sym, ortho_id, sep = " ( "), ' )')) %>% 
    filter(!(CytoGroup == CytoOrtho))
  
  if(generate_files_short){
  write_csv(tmp, paste0('./Cytoscape/group_netfixednewest',file_suffix,'.csv'))
  }
  
  newPlateBook <- full_join(PlateBook, select(groups_df, -num_drugs), 
                            by = c('ortho_sym', 'ortho_id')) %>% 
    mutate(group_id = ifelse(is.na(group_id), ortho_id, group_id), 
           group_sym = ifelse(is.na(group_sym), ortho_sym, group_sym))
  
  # These lines just for visualization of the genes that will become 
  # Representatives in Cytoscape.
  GroupReps <- filter(newPlateBook, ortho_sym != group_sym) %>% 
    select(group_id, group_sym) %>% distinct()
  
  # Adds columns Rep and edgecolour for visualization purposes.
  UnFilteredNetwork <- IndicateRepresentativeGenes(UnFilteredNetwork, GroupReps) %>%
    mutate(edgecolour = if_else(value == 1, 4, 
                                if_else(`"isoform"?` & GMConnection, 3,
                                        if_else(GMConnection, 2, 
                                                if_else(`"isoform"?`, 1, 0)))))
  
  # create numbering for groups from the Unfiltered Network
  if(generate_files_short){
    lst <- UnFilteredNetwork$Var2
    names(lst) <- UnFilteredNetwork$Var1
    
    groups <- Label_subnets(list_of_edges = lst)
  
  # Create summary table of edges with final names and an identifier for the 
  # collinearity graph it is in. Also write to csv the edges filtered above 
  # threshold.
    
    gene_list <- UnFilteredNetwork %>% 
      select(-ortho_id1, -ortho_id2, -num_drugs1, -num_drugs2, -`"isoform"?`, -edgecolour) %>%
      mutate(group_num = groups[[Var1]]) %>% separate(Var1, c('gene_sym1', 'gene_id1', 'extra1'), sep = '[ ][()]') %>%
      separate(Var2, c('gene_sym2', 'gene_id2', 'extra2'), sep = '[ ][()]') %>% 
      mutate(gene_id1 = as.integer(gene_id1), gene_id2 = as.integer(gene_id2)) %>%
      select(-extra1, -extra2, -IsRep, -IsRep2)
    
    gene_list <- left_join(gene_list, group_names, by = c("gene_id1" = "gene_id")) %>%
      rename(group_sym1 = group_sym, group_id1 = group_id)
    gene_list <- left_join(gene_list, group_names, by = c("gene_id2" = "gene_id")) %>%
      rename(group_sym2 = group_sym, group_id2 = group_id)
    
    gene_list <- gene_list %>% rename(JaccardIdx = value)
    
    write_csv(gene_list, paste0('./grouping/GroupingsByEdgesComplete',file_suffix,'.csv'))
    
    gene_list <- gene_list %>% filter(JaccardIdx == 1 | GMConnection == T) %>% 
      rename(group_id = group_id1, group_sym = group_sym1) %>% 
      select(-group_id2, -group_sym2, -group_num)
    
    write_csv(gene_list, paste0('./grouping/GroupingsByEdges',file_suffix,'.csv'))
    
    gene_list <- gene_list %>% mutate(gene1 = paste0(gene_sym1, '(', gene_id1, ')'),
                                      gene2 = paste0(gene_sym2, '(', gene_id2, ')'))
    write_csv(gene_list, paste0('./Cytoscape/GroupingsGraph', file_suffix, '.csv'))
    
    #Copy in Cytoscape folder for plotting
    write_csv(UnFilteredNetwork, 
            paste0('./Cytoscape/Gene_Collinearitycomplete', threshold, file_suffix,'.csv'))
    
    # Write_csv listing the collinearity for all pairs of genes that appear 
    # together at least once in the dataset. Same columns as GroupingsByEdges and
    # added as the second sheet in the grouping xls sent to everyone.
    
    if (generate_files_long){
      GroupingPg2 <- AllgenesPresent %>% filter(value > 0) %>% 
        separate(Var1, c('gene_sym1', 'gene_id1', 'extra1'), sep = '[ ][()]') %>%
        separate(Var2, c('gene_sym3', 'gene_id2', 'extra2'), sep = '[ ][()]') %>%
        select(-ortho_id1, -ortho_id2, -num_drugs1, -num_drugs2, -extra1, -extra2) %>%
        rename(JaccardIdx = value)
      
      write_csv(GroupingPg2, './grouping/GroupingsByEdgesPg2.csv')  
    }
    
    # just the grouped edges
    # remove the EntrezIDs from the gene names to help fit the label inside the node.
    GraphGroups <- UnFilteredNetwork %>% filter(edgecolour > 1) %>% 
      mutate(Var1 = str_split(Var1, ' ', simplify = T)[[1]],
             Var2 = str_split(Var2, ' ', simplify = T)[[1]])

    write_csv(GraphGroups, paste0('./Cytoscape/Gene_CollinearityGraph', threshold, file_suffix, '.csv'))
    
    #Copy in grouping folder
    write_csv(UnFilteredNetwork, 
              paste0('./grouping/Gene_Collinearity', threshold, file_suffix,'.csv'))
  }
  
  newPlateBook
}

IndicateRepresentativeGenes <- function(UnFilteredNetwork, GroupReps){
  
  GroupReps <- mutate(GroupReps, 
                      Var2 = paste(paste(group_sym, group_id, sep = " ( "), ' )', sep = ""),
                      IsRep2 = T) %>% select(-group_sym, -group_id)
  UnFilteredNetwork <- full_join(UnFilteredNetwork, GroupReps, by = c('Var2')) %>%
    mutate(IsRep2 = !is.na(IsRep2))
  
  UnFilteredNetwork <- mutate(UnFilteredNetwork,
                               CpyVar1 = Var1,
                               Var1 = if_else(IsRep2, Var2, Var1), 
                               Var2 = if_else(IsRep2, CpyVar1, Var2)) %>%
                        select(-CpyVar1)
  
  GroupReps <- rename(GroupReps, Var1 = Var2, IsRep = IsRep2)
  
  UnFilteredNetwork <- full_join(UnFilteredNetwork, GroupReps, by = c('Var1')) %>%
    mutate(IsRep = !is.na(IsRep)) %>% filter(!is.na(Var1) & !is.na(Var2))
  
  UnFilteredNetwork
  
}

Label_subnets <- function(list_of_edges){
  # Using DFS it returns a list with the names for each node corresponds to an
  # integer where all nodes with the same integer are connected.
  # Assumes no missing values in the list_of_edges
  visited <- rep(0, length(unique(c(names(list_of_edges), list_of_edges))))
  names(visited) <- unique(c(names(list_of_edges), list_of_edges))
  groups <- vector("list", length(unique(c(names(list_of_edges), list_of_edges))))
  names(groups) <- unique(c(names(list_of_edges), list_of_edges))
  backedge <- names(list_of_edges)
  names(backedge) <- list_of_edges
  
  # traverse the network labeling each bunch of connected nodes as a group.
  group <- 1
  for (name in unique(names(list_of_edges))){
    if (visited[name] == 1) {
      next
    }
    visited[name] <- 1
    groups[name] <- group
    children <- unique(c(list_of_edges[names(list_of_edges) == name], backedge[names(backedge) == name]))
    while (length(children) > 0){
      child <- children[length(children)]
      children <- children[1:length(children)-1]
      visited[child] <- 1
      groups[child] <- group
      new_children <- vector()
      for (new_child in unique(c(children, unique(c(list_of_edges[names(list_of_edges) == child], backedge[names(backedge) == child]))))){
        if (visited[new_child] != 1){
          new_children <- c(new_children, new_child)
        }
      }
      children <- c(children, new_children)
    }
    group <- group + 1
  }
  
  groups
}
