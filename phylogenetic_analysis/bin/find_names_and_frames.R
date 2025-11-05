#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
find_names_and_frames <- function(list_file = "list.csv",
                                  db_file = "nucleotide_list.csv",
                                  orf_name = "none",
                                  translations = "none"){
  the_db <- read_csv(args[2], show_col_types = F)
  the_list <- read.csv(args[1], header=F)
  if (db_file == "nucleotide_list.csv"){
    seqs_type <- "nucleotide_seqs/"
  } else {
    seqs_type <- "aa_seqs/"
  }
  seqs <- list.files(seqs_type)
  tree_file <- ""

  if (length(args) > 2) {
    orf_name <- args[3]
  }
  if (length(args) > 3) {
    translations <- args[4]
  }
  
  # Adds new rows to the annotation file
  for (i in unique(seqs)){
    if (!i %in% unique(unlist(the_db$Accession_id))){
      cat(paste0("Found sequence not in database: ",i, ".\n"))
      the_fasta <- suppressWarnings(read.table(paste0(seqs_type,i),
                                header = F,
                                sep = "\n"))
      the_name <- substr(the_fasta[1,], 2, nchar(the_fasta[1,]))
      chicken <- tibble("Accession_id" = i,
                        "Name" = the_name,
                        "ORFs" = "Unknown",
                        "Tag" = "None")
      the_db <- rbind(the_db, chicken)
    }
  }
  
  # Adds name from fasta file if missing
  for (i in unique(unlist(the_db$Accession_id))){
    if (is.na(the_db[the_db$Accession_id == i,]$Name)){
      cat(paste0("Name for ", i, " not found. Appending from fasta.\n"))
      the_fasta <- suppressWarnings(read.table(paste0(seqs_type,i),
                              header = F,
                              sep = "\n"))
      the_name <- substr(the_fasta[1,], 2, nchar(the_fasta[1,]))
      the_db[the_db$Accession_id == i,]$Name <- the_name
    }
  }

  if (ncol(the_list) > 1 & (is.na(orf_name) | orf_name == "none")){
    cat("Frames detected, generating tree file from subsetted sequences.\n")
    names(the_list) <- c("Acc", "Start", "Stop")
    frames <- T
  } else if (!is.na(orf_name) & orf_name != "none"){
    frames <- T
    # set default start and stops
    the_list$Start <- 0
    the_list$Stop <- 0
    the_list$FS <- 0
    the_list$FS_position <- 0
    names(the_list) <- c("Acc", "Start", "Stop", "FS", "FS_position")
    cat(paste0("Pulling sequences from orf: ", orf_name,".\n"))
    #iterates through the list
    for (i in unique(unlist(the_list$Acc))){
      orf_found <- F
      #splits the ORFs to find the correct one
      orf_string <- strsplit(the_db[the_db$Accession_id == i,]$ORFs, ",")[[1]]
      for (j in orf_string){
        j <- trimws(j)
        j_split <- strsplit(j, "\\[")[[1]][1]
        #If the orf name matches, then its split to find things
        if (j_split == orf_name){
          orf_range <- strsplit(trimws(j), "\\[")[[1]][2]
          orf_range <- strsplit(orf_range, "\\]")[[1]][1]
          orf_range <- strsplit(orf_range, "\\:")[[1]]
          the_list[the_list$Acc == i,]$Start <- as.numeric(orf_range[1])
          the_list[the_list$Acc == i,]$Stop <- as.numeric(orf_range[length(orf_range)])
          if (length(orf_range) == 3){
            #print("Found a frameshift")
            fs_number <- as.numeric(strsplit(orf_range[2], "_")[[1]][1])
            fs_position <- as.numeric(strsplit(orf_range[2], "_")[[1]][2])
            the_list[the_list$Acc == i,]$FS <- fs_number
            the_list[the_list$Acc == i,]$FS_position <- fs_position
          }
          orf_found <- T
        }
      }
      if (!orf_found){
        cat(paste0("Orf ", orf_name, " not found in ", i,". Make sure annotation is formatted correctly or remove from list.\n"))
        the_list <- the_list[the_list$Acc != i,]
      }
    }
  } else {
    cat("Generating tree file using the full sequences.\n")
    names(the_list) <- c("Acc")
    the_list$Start <- 0
    the_list$Stop <- 0
    the_list$FS <- 0
    the_list$FS_position <- 0
    frames <- F
  }
  # Removes duplicates from the list
  if (TRUE %in% unique(duplicated(the_list))){
    cat("Duplicates detected. Removing...\n")
    the_list <- the_list[!duplicated(the_list),]
  }
  # iterates through the list to create the subset
  for (i in unique(unlist(the_list$Acc))){
    #print(the_db[the_db$Accession_id == i,]$Name)
    tree_file <- paste0(tree_file, ">", the_db[the_db$Accession_id == i,]$Name, "_", i, "\n")
    
    # Pulls in the read sequence
    the_read <- ""
    the_fasta <- suppressWarnings(read.table(paste0(seqs_type,i),
                            header = F,
                            sep = "\n"))
    for (j in c(2:nrow(the_fasta))){
      the_read <- paste0(the_read, the_fasta[j,1])
    }
    
    # sets the start and stop defaults
    the_start <- 1
    the_end <- nchar(the_read)
    #print(the_end)
    
    # Stop and start sites are taken from the file
    if (frames == T & 
        !is.na(the_list[the_list$Acc == i,]$Start) & 
        the_list[the_list$Acc == i,]$Start > 0){
      the_start <- the_list[the_list$Acc == i,]$Start
    } 
    
    if (frames == T & !is.na(the_list[the_list$Acc == i,]$Stop) & the_list[the_list$Acc == i,]$Stop > 0){
      the_end <- the_list[the_list$Acc == i,]$Stop
    }
    
    if (the_list[the_list$Acc == i,]$FS != 0 & translations != "none"){
      the_read <- substr(the_read, the_start, the_end)
      cat(paste("Found a frameshift in", i, "of", the_list[the_list$Acc == i,]$FS, "at", the_list[the_list$Acc == i,]$FS_position, "\n"))
      shifted <- 0
      for (j in c(1:nchar(the_read))){
        if (j == the_list[the_list$Acc == i,]$FS_position-the_start+1){
          #print("ding")
          #print(j)
          shifted <- the_list[the_list$Acc == i,]$FS
          #print(j+shifted)
          #print(substr(the_read, j+shifted, j+shifted))
        }
        tree_file <- paste0(tree_file, substr(the_read, j+shifted, j+shifted))
      }
      tree_file <- paste0(tree_file, "\n")
    } else {
      tree_file <- paste0(tree_file, substr(the_read, the_start, the_end), "\n")
    }
  }
  #print(tree_file)
  cat(tree_file, file="sequences.fasta")
  write_csv(the_db, db_file)
}

find_names_and_frames(list_file = args[1], db_file = args[2], orf_name = args[3])
