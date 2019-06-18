# trying out gProfiler R package
library(tidyverse)
library(gProfileR)

my_gprofiler <- function(file, num_sig, seed = NULL, block_size = NULL, col_types = NULL, correction_method = "bonferroni"){
  # Function to read in an ordered list of genes and run the gprofiler algorithm 
  # to find any enrichment of the found genes.
  # Args:
  # file is an ordered list of Entrez_ids separated by a newline character where the 
  # first line is group_id
  # num_sig is the number of genes from the top of the list taken to have been 
  # found significant, as a percentage or absolute number of genes.
  # Set seed to not null to take a random ordering of the genes with the given 
  # seed, if seed is null it will not reorder the input.
  # col_types is the same arg passed to read_csv.
  gene_list <- read_csv(file, col_types = col_types)
  
  if (num_sig < 0){
    stop("num_sig can't be negative")
  } else if(num_sig < 1){
    num_sig <- length(gene_list[[1]])*num_sig
  }
  
  if (is_null(num_sig)){
    num_sig <- length(gene_list[[1]])
  }
  
  if(!is_null(seed)){
    if (is_null(block_size)){
      block_size == length(gene_list[[1]])
    }
    set.seed(seed)
    gene_list <- gene_list %>% mutate(block = row_number() %% block_size) %>% 
      group_by(block) %>% group_by(member = row_number()) %>% sample_frac(1) %>%
      ungroup() %>% select(-member, -block)
  }
  
  #if (!names(gene_list) == 'group_id'){
  #  stop('header not group_id, file may not be formatted correctly')
  #}
  
  gprofiler(gene_list[[1]][1:num_sig], ordered_query = T, 
            custom_bg = gene_list[[1]], src_filter = c('KEGG', "REAC"),
            correction_method = correction_method, numeric_ns = 'ENTREZGENE_ACC')

}

# Removed the non-interpretable groups before writing this list of genes, 
# hence the name OrderedInterpretable...
dir1 <- './Coefs/Penalty/CellNum/OrderedInterpretableMainEffects/Syms/'
dir2 <- './Coefs/Penalty/CellSize/OrderedInterpretableMainEffects/Syms/'
files <- paste0(dir1, list.files(dir1))
files <- c(files, paste0(dir2, list.files(dir2)))

files <- files[grep('Inter:FALSE.*[10,20|1,2|3,6]lambda', files)]

# Test different thresholds with this for loop
Visualize_changing_paths_over_thresholds <- function(files, num_trials){
  # starts from keeping 5 genes as important, and incrementing by 5 as many 
  # times as indicated by num_trials.
  for (i in 1:length(files)){
    f <- files[[i]]
    if (length(grep('Syms/', f)) == 1){
      c_t = 'c'
    } else {
      c_t = NULL
    }
    results <- vector('list', num_trials)
    for(i in seq(1, num_trials)){
      t <- i*5
      results[[i]] <- as.tibble(my_gprofiler(f, t, col_types = c_t))$term.name
      names(results)[[i]] <- t
      print(t)
    }
    
    test <- lapply(results, function(x){length(x) <- max(lengths(results)); x})
    paths <- unlist(results) %>% unique(.)
    
    # Visualize the pathways for a given threshold applied to the ordered list of 
    # genes found to be significant to see how this choice affects the pathways 
    # found to be enriched.
    
    ncols <- length(test)
    nrows <- length(paths)
    mat <- matrix(nrow = nrows, ncol = ncols)
    colnames(mat) <- names(test)
    rownames(mat) <- paths
    for (thresh in names(test)){
      for(path in paths){
        mat[path, thresh] <- if_else(path %in% test[[thresh]], 1, 0)
      }
    }
    plot_name <- gsub('.*Syms_|MinDrgs.*concs:|_rmv.*|_1|_In.*', '', f)
    test <- as.tibble(mat, rownames = 'pathway') %>% 
      gather(key = threshold, value = present, -pathway) %>% 
      mutate(threshold = as.integer(threshold), present = as.logical(present))
    plt <- ggplot(test, aes(threshold, pathway)) + geom_tile(aes(fill = present), colour = 'black')+
      scale_fill_manual(values = c("white", "dodgerblue2"), breaks = c(T, F)) + 
      ggtitle(plot_name)
    
    plot(plt)
  }
}

Visualize_changing_paths_over_thresholds(files, 30)

# maybe don't need this whole function...
Visualize_enrich_permute_local_blocks <- function(files, t, num_trials, num_seeds){
  # for each file and at a given threshold, it will run gprofiler where the 
  # order is randomly permuted within local blocks, for a given num_trials it 
  # will increase the size of the blocks by 1 that many times. num_trials cannot 
  # be larger than t.
  for (i in 1:length(files)){
    f <- files[[i]]
    if (length(grep('Syms/', f)) == 1){
      c_t = 'c'
    } else {
      c_t = NULL
    }
    
    #browser()
    results <- vector('list', num_trials)
    for (j in seq(1,num_trials)) {
      for(k in seq(1, num_seeds)){
        results[[j]] <- c(results[[j]], as.tibble(my_gprofiler(f, t, seed = 7+k, block_size = j, 
                                               col_types = c_t))$term.name)
        names(results)[[j]] <- j
      }
    }
    test <- lapply(results, function(x){length(x) <- max(lengths(results)); x})
    paths <- unlist(results) %>% unique(.)
    
    # Visualize the pathways for a given threshold applied to the ordered list of 
    # genes found to be significant to see how this choice affects the pathways 
    # found to be enriched.
     
    #browser()
    
    ncols <- length(test)
    nrows <- length(paths)
    mat <- matrix(nrow = nrows, ncol = ncols)
    colnames(mat) <- names(test)
    rownames(mat) <- paths
    for (block_size in names(test)){
      for(path in paths){
        mat[path, block_size] <- sum(path == test[[block_size]], na.rm = T) #if_else(path %in% test[[block_size]], mat[path, block_size], 0)
      }
    }
    
    
    plot_name <- gsub('.*Syms_|MinDrgs.*concs:|_rmv.*|_1|_In.*', '', f)
    if (length(mat) != 0){
      test <- as.tibble(mat, rownames = 'pathway') %>% 
        gather(key = block_size, value = present, -pathway) %>% 
        mutate(block_size = as.integer(block_size))
    
      plt <- ggplot(test, aes(block_size, pathway)) + geom_tile(aes(fill = present), colour = 'black') + 
        scale_fill_gradient2(low = "White", mid= "blue", high = "purple", midpoint = 5) + ggtitle(plot_name)
      
      plot(plt)
    } else {
      print(paste0(plot_name,': Has no enriched pathways at the given threshold'))
    }
  }
}

Visualize_enrich_permute_local_blocks(files, 0.1, 10, 10)

Get_gprofile_results <- function(files, t){
  # Writes a csv in the gprofiler folder for each Ordered Syms files in the list
  #  files, using t as a threshold for how many in the list are deemed 
  #  important.
  results <- vector('list', length(files))
  for (i in 1:length(files)){
    f <- files[[i]]
    if (length(grep('Syms/', f)) == 1){
      c_t = 'c'
    } else {
      c_t = NULL
    }
    results[[i]] <- as.tibble(my_gprofiler(f, t, col_types = c_t, correction_method = 'analytical'))
    # when including 1se syms change this gsub pattern, so that the name reflects 
    # the different lambda choices also.
    names(results)[[i]] <- gsub('.*Syms_|MinDrgs.*concs:|_rmv.*|_1|_In.*', '', f)
    
    file_name <- sub('.*?_', '', f) %>% gsub('rmvdups.*.txt|_Interp.txt', '', .)
    write_csv(results[[i]], paste0('./gprofiler/', file_name, Sys.Date(), 'NoCorrection.csv'))
  }
  
  
}

Get_gprofile_results(files, 0.1)
