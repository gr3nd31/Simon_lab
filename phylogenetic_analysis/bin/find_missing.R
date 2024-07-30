#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
find_missing <- function(db_file = "nucleotide_list.csv",
                    list_file = "list.csv"){
  the_db <- read_csv(db_file, show_col_types = F)
  the_list <- read.csv(list_file, header = F)
  pull_list <- ""
  for (i in unique(unlist(the_list[,1]))){
    if (!i %in% unique(the_db$Accession_id)){
      pull_list <- paste0(pull_list, i, "\n")
    }
  }
  if (pull_list == ""){
    cat("All sequences present.\n")
  } else {
    cat("Found new sequences to download.\n")
    cat(pull_list, file = "pull_list.txt")
  }
}
find_missing(list_file = args[1], db_file = args[2])
