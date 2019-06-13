# ".gpml" file helpers


#Extract all unique DB entries for a gpml file.
# filepath : string which is the path to the file to read from
# pattern : string, regex pattern to match on each line.
# group : which group to return, 1 is the whole line, 2 is the first group etc.
get_unique_matches <- function(filepath, pattern = 'Database=\\"(.*?)\\"', group){
  stopifnot(require(readr))
  c <- file(filepath, "r")
  on.exit(close(c))
  fline <- readLines(con = c, n = 1)
  dbs <- vector("character", length = 0)
  while(length(fline) != 0){
    db <- str_match(fline, pattern)[1,group]
    dbs <- c(dbs,db)
    fline <- readLines(con = c, n = 1)
  }
  unique(dbs)
}

#Reading a GPML file, stores node name, database, and id information in a tibble
# filepath : string, path to read gpml file
get_Nodes <- function(filepath){
  stopifnot(require(readr) && require(tidyverse))
  c <- file(filepath, "r")
  on.exit(close(c))
  fline <- readLines(con = c, n = 1)
  nodes <- character() #list()
  dbs <- character() #list()
  ids <- character() #list()
  while(length(fline) != 0){
    if (grepl('<DataNode ', fline)){
      node <- str_match(fline, 'TextLabel=\\"(.*?)\\"')[1,2]
      fline <- readLines(con = c, n = 1)
      while(!grepl('</DataNode>', fline)){
        db <- str_match(fline, 'Database=\\"(.*?)\\"')
        id <- str_match(fline, 'ID=\\"(.*?)\\"')
        if (!is.na(db[1,1])){
          dbs <- c(dbs, db[1,2])
          ids <- c(ids, id[1,2])
        }
        fline <- readLines(con = c, n = 1)
      }
      nodes <- c(nodes, node)
      # Sanity check, if length of nodes is greater than dbs add an empty string
      # to dbs, i.e. assuming this means that node was missing a database 
      # attribute.
      if (length(nodes) != length(dbs)){
        dbs <- c(dbs, '')
        ids <- c(ids, '')
      }
    }
    fline <- readLines(con = c, n = 1)
  }
  tibble(name = nodes, database = dbs, id = ids)
}
