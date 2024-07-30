#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
find_names_and_frames <- function(list_file = "list.csv",
                                  db_file = "nucleotide_list.csv"){
  the_db <- read_csv(db_file, show_col_types = F)
  the_list <- read.csv(list_file, header=F)
  seqs <- list.files("nucleotide_seqs/")
  tree_file <- ""
  
  for (i in unique(seqs)){
    if (!i %in% unique(unlist(the_db$Accession_id))){
      cat(paste0("Found sequence not int database: ",i, ".\n"))
      the_fasta <- read.table(paste0("nucleotide_seqs/",i),
                                header = F,
                                sep = "\n")
      the_name <- substr(the_fasta[1,], 2, nchar(the_fasta[1,]))
      chicken <- tibble("Accession_id" = i,
                          "Name" = the_name)
      the_db <- rbind(the_db, chicken)
    }
  }
  for (i in unique(unlist(the_db$Accession_id))){
    #print(i)
    if (is.na(the_db[the_db$Accession_id == i,]$Name)){
      cat(paste0("Name for ", i, " not found. Appending from fasta.\n"))
      the_fasta <- read.table(paste0("nucleotide_seqs/",i),
                              header = F,
                              sep = "\n")
      the_name <- substr(the_fasta[1,], 2, nchar(the_fasta[1,]))
      the_db[the_db$Accession_id == i,]$Name <- the_name
    }
  }
  if (ncol(the_list) > 1){
    cat("Frames detected, generating tree file from subsetted sequences.\n")
    names(the_list) <- c("Acc", "Start", "Stop")
    frames <- T
  } else {
    cat("Generating tree file using the full sequences.\n")
    names(the_list) <- c("Acc")
    the_list$Start <- 0
    the_list$Stop <- 0
    frames <- F
  }
  if (TRUE %in% unique(duplicated(the_list))){
    cat("Duplicates detected. Removing...\n")
    the_list <- the_list[!duplicated(the_list),]
  }
  for (i in unique(unlist(the_list$Acc))){
    #print(the_db[the_db$Accession_id == i,]$Name)
    tree_file <- paste0(tree_file, ">", the_db[the_db$Accession_id == i,]$Name, "\n")
    the_read <- ""
    the_start <- 0
    the_fasta <- read.table(paste0("nucleotide_seqs/",i),
                            header = F,
                            sep = "\n")
    for (j in c(2:nrow(the_fasta))){
      the_read <- paste0(the_read, the_fasta[j,1])
    }
    the_end <- nchar(the_read)
    if (frames == T & !is.na(the_list[the_list$Acc == i,]$Start)){
      the_start <- the_list[the_list$Acc == i,]$Start
    } 
    
    if (frames == T & !is.na(the_list[the_list$Acc == i,]$Stop)){
      the_end <- the_list[the_list$Acc == i,]$Stop
    }
    tree_file <- paste0(tree_file, substr(the_read, the_start, the_end), "\n")
  }
  cat(tree_file, file="sequences.fasta")
  write_csv(the_db, db_file)
}

find_names_and_frames(list_file = args[1], db_file = args[2])
