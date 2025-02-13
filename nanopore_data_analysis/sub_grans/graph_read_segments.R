suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plotly))

graph_read_segments <- function(file_name = "reads_df.csv",
                                as_percent = F,
                                color_by = "sense",
                                save_it = T,
                                save_as = "read_hsps.png",
                                subset_data = F,
                                subset_number = 100,
                                skips = T,
                                plotly_it = F){
  if (typeof(file_name)=="character"){
    datum <- read_csv(file_name, show_col_types = F)
  } else if (typeof(file_name)=="list"){
    datum <- file_name
  }
  
  if (subset_data & typeof(subset_number) == "double"){
    print("Sampling data...")
    read_ids <- unique(datum$read_id)
    read_ids <- sample(read_ids, subset_number)
    datum <- datum[datum$read_id %in% read_ids,]
  }
  if (as_percent){
    print("Adjusting reads as percent.")
  }
  datum <- datum[order(-datum$aligned_length),]
  y_tick <- 1
  datum$y1 <- 0
  for (i in unique(datum$read_id)){
    datum[datum$read_id == i,]$y1 <- y_tick
    if(as_percent){
      datum[datum$read_id == i,]$q_start <- 1+as.integer(100*datum[datum$read_id == i,]$q_start/max(datum[datum$read_id ==i,]$read_length))
      datum[datum$read_id == i,]$q_end <- as.integer(100*datum[datum$read_id == i,]$q_end/max(datum[datum$read_id ==i,]$read_length))
    }
    chicken <- tibble("read_id" = i,
                      "start" = 1,
                      "end" = max(datum[datum$read_id == i,]$read_length),
                      "y1" = y_tick)
    if (!exists("interim")){
      interim <- chicken
    } else {
      interim <- rbind(interim, chicken)
    }
    y_tick <- y_tick + 1
  }
  if (as_percent){
    interim$end <- 100
  }
  #interim$y1 <- as.character(interim$y1)
  #datum$y1 <- as.character(datum$y1)
  
  draft <- ggplot()+
    geom_segment(data = interim, aes(x = start, y = y1, xend = end, yend = y1),
                 color = "darkgray", linewidth = 5, alpha = 1)
  if (color_by == "sense"){
    draft <- draft+
      geom_segment(data = datum, aes(x = q_start, y = y1, xend = q_end, yend = y1, color = h_strand),
                   linewidth = 3, alpha = 1)
  } else if (color_by == "start"){
    draft <- draft+
      geom_segment(data = datum, aes(x = q_start, y = y1, xend = q_end, yend = y1, color = h_start),
                   linewidth = 3, alpha = 1)+
      scale_color_viridis_b(option = "inferno")
  } else if (color_by == "end"){
    draft <- draft+
      geom_segment(data = datum, aes(x = q_start, y = y1, xend = q_end, yend = y1, color = h_end),
                   linewidth = 3, alpha = 1)+
      scale_color_viridis_b(option = "inferno")
  }
  draft <- draft+
    labs(x = "Read Position",
         y = "Read Number")+
    theme_bw()+
    theme(axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(), 
          axis.ticks.y = element_blank(),
          legend.position = "top")
  print(draft)
  if (save_it & nchar(save_as) > 0){
    ggsave(save_as, dpi = 300)
  }
}

graph_read_segments(file_name = "reads_df.csv",
                    subset_data = T, 
                    subset_number = 50, 
                    as_percent = F,
                    color_by = "start")

