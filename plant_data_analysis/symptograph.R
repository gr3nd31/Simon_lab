library(tidyverse)
library(ggpubr)

symptograph <- function(file_name = "test.csv",
                        graph_set = c(),
                        age_range = c(),
                        min_age = T,
                        percent_at_100 = T,
                        limit = 50,
                        line_at_limit = T,
                        alph = 0.6,
                        prefix = "",
                        format = "png"){
  datum <- read_csv(file_name)
  datum$T_date <- as.Date(datum$T_date, "%m/%d/%y")
  datum$I_date <- as.Date(datum$I_date, "%m/%d/%y")
  datum$S_date <- as.Date(datum$S_date, "%m/%d/%y")
  datum$F_date <- as.Date(datum$F_date, "%m/%d/%y")
  datum$age_infiltrated <- 0
  datum$age <- 0
  datum$dpi <- 0
  datum$dus <- 0
  
  for (i in unique(datum$I_date)){
    if (length(unique(datum[datum$I_date == i,]$F_date)) < 2 & 
        !is.na(unique(datum[datum$I_date == i,]$F_date))){
      final_date <- max(datum[datum$I_date ==i,]$F_date)
    } else {
      final_date <- as.Date(today())
    }
    datum[datum$I_date ==i,]$age_infiltrated <- as.numeric(datum[datum$I_date ==i,]$I_date - datum[datum$I_date ==i,]$T_date)
    datum[datum$I_date ==i,]$age <- as.numeric(final_date - datum[datum$I_date ==i,]$T_date)
    datum[datum$I_date ==i,]$dpi <- as.numeric(final_date - datum[datum$I_date ==i,]$I_date)
    datum[datum$I_date ==i,]$dus <- as.numeric(datum[datum$I_date ==i,]$S_date - datum[datum$I_date ==i,]$I_date)
    datum[is.na(datum$S_date) & datum$I_date ==i,]$dus <- max(datum[datum$I_date ==i,]$age)+10
  }
  write_csv(datum, file_name)
  
  for (j in unique(datum$Infection)){
    trigger_date <- 0
    if (min_age){
      day_set <- min(datum$dpi)
    } else {
      day_set <- max(datum$dpi)
    }
    for (i in c(0:day_set)){
      interim <- datum[datum$Infection == j,]
      gin <- nrow(interim)
      
      sym_count <- nrow(interim[interim$dus <= i & !is.na(interim$S_date),])
      
      sym_percent <- 100*sym_count/nrow(interim)
      if (sym_percent >= limit & trigger_date == 0){
        trigger_date <- i 
      }
      
      interim_times <- tibble("Day" = i,
                              #"Final_age" = unique(interim$age),
                              #"Infiltration_age" = unique(interim$age_infiltrated),
                              "Sample" = j,
                              "Percent" = sym_percent,
                              "Total" = gin,
                              "Subset_total"= sym_count)
      if (!exists("all_times")){
        all_times <- interim_times
      } else {
        all_times <- rbind(all_times, interim_times)
      }
      if (i == day_set){
        max_percent <- sym_percent
      }
    }
    id50 <- tibble("Infection" = j,
                   "Total" = gin,
                   "Percent" = max_percent,
                   "days_to_50" = trigger_date)
    if (!exists("trigger_days")){
      trigger_days <- id50
    } else {
      trigger_days <- rbind(trigger_days, id50)
    }
  }
  
  write_csv(all_times, paste0("timed_", file_name))

  if (length(graph_set)> 0){
    trick <- all_times[all_times$Sample %in% graph_set,]
  } else {
    trick <- all_times
  }
  draft <- ggplot(data = trick, aes(x = Day,
                                    y = Percent,
                                    color = Sample))+
    geom_step(linewidth=1.5, alpha=alph)+
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
  if (line_at_limit){
    draft <- draft +
      geom_hline(yintercept = limit, color = "red", alpha = 0.4, linetype = 2, linewidth = 1)
  }
  if (length(age_range) == 2){
    draft <- draft+xlim(age_range[1], age_range[2])
  }
  if (percent_at_100){
    draft <- draft+ylim(0,100)
  }
  print(draft)
  ggsave(paste0(prefix, str_replace(file_name, "csv", format)), width = 6, height = 5)
  print(trigger_days)
  write_csv(trigger_days, paste0(prefix, "id", limit,".csv"))
}
