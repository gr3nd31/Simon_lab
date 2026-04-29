#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
MASSive <- function(alignmentFile="db_aligned.tsv",
                    genomeFile="sequence.fasta",
                    sirnaFile="sirna.csv",
                    threshold=50,
                    outFile="FilteredsiRNA.csv"){
  siRNA <- read_csv(sirnaFile, show_col_types = F)
  alignmentFile <- read_tsv(alignmentFile, show_col_types = F, col_names = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'))
  z <- readLines(genomeFile)
  gCount=sum(str_count(z, ">"))
  
  siRNA$tab <- paste0(str_replace_all(siRNA$Name, " ", "_"), "_", siRNA$Position)
  siRNA$TotalHit <- 0
  siRNA$HundredHit <- 0
  for(i in unique(siRNA$tab)){
    siRNA[siRNA$tab == i,]$TotalHit <- length(unique(alignmentFile[alignmentFile$qseqid ==i,]$sseqid))
    siRNA[siRNA$tab == i,]$HundredHit <- length(unique(alignmentFile[alignmentFile$qseqid ==i & alignmentFile$pident == 100,]$sseqid))
  }
  siRNA$TotalHit <- round(100*(siRNA$TotalHit/gCount), 2)
  siRNA$HundredHit <- round(100*(siRNA$HundredHit/gCount), 2)
  siRNA <- siRNA[siRNA$HundredHit > threshold,]
  write_csv(siRNA[-11], outFile)
}
MASSive(alignmentFile = args[1], genomeFile = args[2], sirnaFile = args[3])