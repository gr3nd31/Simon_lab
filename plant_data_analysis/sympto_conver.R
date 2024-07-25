library(tidyverse)

sympto_convert <- function(input,
                           t_date = "", 
                           i_date = "",
                           s_date = "", 
                           f_date = "") {
  if (!input %in% list.files()){
    print("File not found. Try again")
  } else {
    interim <- read_csv(input)
    #print(head(interim))
  for (i in c(1:nrow(interim))){
    plant_num <- interim[i,"Number"]
    b_rep <- unique(interim[i,]$Replicate)
    #print(plant_num)
    for (j in c(1:as.numeric(plant_num))){
      sympto_interim <- tibble("Infection" = as.character(interim[i,1]),
                               "Plant" = as.character(interim[i,"Plant"]),
                               "bio_replicate" = b_rep,
                               "Replicate" = j,
                               "T_date" = t_date, 
                               "I_date" = i_date, 
                               "S_date" = s_date, 
                               "F_date" = f_date, 
                               "Notes" = "")
      if (!exists("sympto_log")){
        sympto_log <- sympto_interim
      } else {
        sympto_log <- rbind(sympto_log, sympto_interim)
      }
    }
  }
  new_input <- paste0(strsplit(input, ".csv")[[1]][1], "_symptoLog.csv")
  #print(new_input)
  #print(sympto_log)
  write_csv(sympto_log, as.character(new_input))
  }
}

sympto_convert("infection_plan.csv", t_date = "2024.06.01", i_date = "2024.07.04", f_date = "2024.08.15")
