suppressPackageStartupMessages(library(tidyverse))
datum <- read_tsv("db_aligned.tsv", show_col_types = F, col_names = c("read_id", "ref", "percIdent",
                                                                      "length", "X", "Y", "q_start",
                                                                      "q_end", "h_start", "h_end",
                                                                      "expect", "score"))
read_db <- tibble("ref" = unique(datum$ref))
read_db$reads <- ""

if (!file.exists("reads")){
  dir.create("reads")
}
setwd("reads")
for (i in unique(datum$read_id)){
  interim <- datum[datum$read_id == i,]
  if (length(unique(interim$ref)) == 1){
    the_ref <- unique(interim$ref)[1]
  } else {
    if (max(table(interim$ref)) == 1 & nrow(interim[interim$expect == min(interim$expect),]) < 2){
      the_ref <- interim[interim$expect == min(interim$expect),]$ref
    } else if (length(table(interim$ref)[table(interim$ref) == max(table(interim$ref))]) == 1) {
      the_ref <- names(sort(table(interim$ref), decreasing = T))[1]
    } else {
      the_ref <- "multiple"
      if (!file.exists("multiple")){
        dir.create("multiple")
      }
      write(i, file = "multiple/multiple.txt", append = TRUE)
      #j <- interim
      the_ref <- interim[interim$expect == min(interim$expect),]$ref
      if (length(the_ref) > 1){
        the_ref <- the_ref[1]
        j <- interim
      }
    }
  }
  if (!file.exists(the_ref)){
    dir.create(the_ref)
  }
  write(i, file = paste0(the_ref, "/", the_ref, ".txt"), append = TRUE)
}
setwd("../")
