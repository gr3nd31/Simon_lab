library(tidyverse)
library(ggpubr)

symptograph <- function(file_name = "test.csv",
                       end_date = T,
                       format = "png"){
  datum <- read_csv(file_name)
  datum$T_date <- as.Date(datum$T_date, "%m/%d/%y")
  datum$I_date <- as.Date(datum$I_date, "%m/%d/%y")
  datum$S_date <- as.Date(datum$S_date, "%m/%d/%y")
  if (end_date){
    final_date <- max(datum[!is.na(datum$S_date),]$S_date)
  } else {
    final_date <- as.Date(today())
  }
  datum$age_infiltrated <- datum$I_date - datum$T_date
  datum$age <- as.Date(final_date) - datum$T_date
  datum$dpi <- as.Date(final_date) - datum$I_date
  datum$dus <- as.Date(final_date) - datum$S_date
  datum[is.na(datum$S_date),]$dus <- max(datum$age)+10
  
  for (i in c(0:as.numeric(max(datum$age)))){
    #print(paste0("Day: ", i))
    for (j in unique(datum$Infection)){
      interim <- datum[datum$Infection ==j,]
      sym_count <- nrow(interim[interim$S_date <= min(interim$I_date)+i & !is.na(interim$S_date),])
      sym_percent <- 100*sym_count/nrow(interim)
      #print(paste0(j, ": ", round(sym_percent, 2), "%"))
      interim_times <- tibble("Day" = i,
                              "Final_age" = unique(interim$age),
                              "Infiltration_age" = unique(interim$age_infiltrated),
                              "Sample" = j,
                              "Percent" = sym_percent,
                              "Total" = nrow(interim),
                              "Subset_total"= sym_count)
      if (!exists("all_times")){
        all_times <- interim_times
      } else {
        all_times <- rbind(all_times, interim_times)
      }
    }
  }
  trick <<- all_times
  write_csv(all_times, paste0("timed_", file_name))
  draft <- ggplot(data = trick, aes(x = Day,
                                        y = Percent,
                                        color = Sample))+
    geom_step(linewidth=1.5, alpha=0.6)+
    labs(x = "Day Post-Infiltration",
         y = "Percent Symptomatic")+
    theme_bw()+
    theme(legend.position = "top",
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 16,
                                   color = "black"),
          axis.line = element_line(linewidth = 1),
          axis.ticks = element_line(linewidth = 1,
                                    color = "black"),
          panel.border = element_rect(fill = NA, 
                                      color = "black", 
                                      linewidth = 2),
          legend.text = element_text(size = 14),
          legend.title = element_blank())
  print(draft)
  ggsave(paste0("timed_", str_replace(file_name, "csv", format)), width = 6, height = 5)
}
