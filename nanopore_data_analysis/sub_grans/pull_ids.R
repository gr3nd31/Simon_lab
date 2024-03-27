pull_ids <- function(df,
                     file_name = "reads.txt"){
  read_list = ""
  for (unique_read in unique(df$read_id)){
    check_id <- strsplit(unique_read, " ")[[1]]
    read_id <- check_id[1]
    read_list <- paste0(read_list, read_id, "\n")
  }
  cat(read_list, file = file_name)
}
