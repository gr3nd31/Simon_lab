suppressPackageStartupMessages(library(tidyverse))

pull_and_thull <- function(dataFile = "data.csv"){
  if (!"./subs" %in% list.dirs()){
    dir.create("subs")
  }
  datum <- read_csv(dataFile, show_col_types = F)
  setwd("subs")
  for (i in c(1:nrow(datum))){
    nameo <- strsplit(datum[i,]$Name, " ")[[1]][1]
    nameo <- str_replace(nameo, ">", "")
    cat(datum[i,]$Sequence, file = paste0(nameo, ".seq"))
    cat(datum[i,]$Structure, file = paste0(nameo, ".str"))
  }
  setwd("../")
}
pull_and_thull()
