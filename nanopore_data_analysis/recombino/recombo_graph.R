library(tidyverse)
library(ggpubr)
library(plotly)
library(htmlwidgets)

recombo_graph <- function(path="recombo",
                          save_it = T,
                          graph_it = T,
                          save_widget = T,
                          prefix = ""){
  read_list <- list.dirs(path = path)
  if(exists("datum")){
    rm(datum)
  }
  print("Generating read graphs...")
  for (i in read_list){
    #print(i)
    #print(list.files(path = i))
    if (i != "recombo" & "windowed.csv" %in% list.files(path = i)){
      interim <- read_csv(paste0(i, "/windowed.csv"), show_col_types = F)
      interim$rel_score <- 0
      for (j in unique(interim$window)){
        interim[interim$window ==j,]$rel_score <- 100*interim[interim$window ==j,]$score/max(interim[interim$window ==j,]$score)
      }
      interim$window_start <- interim$window_start-min(interim$window_start)+1
      if (!exists("datum")){
        datum <- interim
      } else {
        datum <- rbind(datum, interim)
      }
      if (graph_it){
        draft <- ggplot(data = interim, aes(x = window_start, y = rel_score, color = reference))+
          geom_line(linewidth = 2)+
          theme_bw()+
          labs(x = "Read position", y = "Relative alignment score")+
          scale_color_viridis_d()+
          ylim(0,100)+
          theme(legend.position = "top",
                axis.title = element_text(size = 16),
                axis.text = element_text(size = 14))
        print(draft)
        ggsave(paste0(i, "/window.png"), dpi=300, width = 15, height = 4, units = "in")
      }
    }
  }
  if (save_it){
    print("Saving total read file...")
    print(head(datum))
    write_csv(datum, paste0(prefix, "windows_df.csv"))
  }
  if(exists("cheeky")){
    rm(cheeky)
  }
  if (exists("summary_df")){
    rm(summary_df)
  }
  print("Generating summary file...")
  for (i in unique(datum$reference)){
    for (j in c(min(datum[datum$reference == i,]$hit_start):max(datum[datum$reference == i,]$hit_start))){
      cheeky <- tibble("Reference" = i,
                        "Position" = j,
                        "Mean_Score" = mean(datum[datum$reference ==i & datum$hit_start==j,]$rel_score),
                        "SD_Score" = sd(datum[datum$reference ==i & datum$hit_start==j,]$rel_score))
      if(!exists("summary_df")){
        summary_df <- cheeky
      } else {
        summary_df <- rbind(summary_df, cheeky)
      }
    }
  }
  write_csv(summary_df, "summary_df.csv")
  if (graph_it){
    draft <- ggplot(data = datum, aes(x = hit_start, y = rel_score, color = reference))+
      geom_line(linewidth = 1)+
      theme_bw()+
      labs(x = "Reference position", y = "Relative alignment score")+
      scale_color_viridis_d()+
      facet_wrap(~reference, ncol = 1)+
      ylim(0,100)+
      theme(legend.position = "none",
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
    print(draft)
    ggsave("windows.png", dpi=300, width = 15, height = 14, units = "in")
    if (save_widget){
      saveWidget(ggplotly(draft), "windows.html", selfcontained = T)
    }
    draft <- ggplot(data = summary_df, aes(x = Position, y = Mean_Score, color = Reference))+
      geom_bar(stat = "identity")+
      theme_bw()+
      labs(x = "Reference position", y = "Relative alignment score")+
      scale_color_viridis_d()+
      facet_wrap(~Reference, ncol = 1)+
      ylim(0,100)+
      theme(legend.position = "none",
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
    draft
  }
}

recombo_graph()
