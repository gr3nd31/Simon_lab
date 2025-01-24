suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))

genome_view <- function(data = "data.csv",
                        type = "APE",
                        infer_name = T,
                        single_fold = F,
                        save_it = F){
  if (typeof(data) == "character"){
    datum <- read_csv(data, show_col_types = F)
  } else {
    datum <- data
  }
  datum$real_name <- "HOLD"
  datum$start <- 0
  datum$end <- 0
  datum$position <- 0
  for (i in unique(datum$Name)){
    datum[datum$Name == i,]$real_name <- str_split(i, "_")[[1]][1]
    datum[datum$Name == i,]$start <- as.numeric(str_split(i, "_")[[1]][2])
    datum[datum$Name == i,]$end <- as.numeric(str_split(i, "_")[[1]][3])
  }
  if (!single_fold){
    for (i in unique(datum$real_name)){
      if (!infer_name){
        iter <- readline(prompt = paste0("What should the name for '", i, "' be: "))
        datum[datum$real_name == i,]$real_name <- iter
      } else {
        iter <- i
      }
      datum[datum$real_name == iter,]$position <- 1000*round(datum[datum$real_name == iter,]$end/max(datum[datum$real_name == iter,]$end),4)
      draft <- ggbarplot(data = datum[grepl("\\(", datum$Structure),], x = "position", 
                         y = "Length", fill = type, color = type)+
        scale_color_viridis_b()+
        scale_fill_viridis_b()+
        labs(x = "Relative Position",
             y = "Hairpin Length",
             title = iter)+
        theme_bw()+
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid = element_blank(),
              legend.position = "top")
      print(draft)
      if(save_it){
        ggsave(paste0(type),"_",iter,".png", dpi = 300)
      }
    }
  } else {
    if (type != "APE"){
      print("Reminder, single folds can only graph APE. Graphing APE now...")
    }
    for (i in unique(datum$real_name)){
      if (!infer_name){
        iter <- readline(prompt = paste0("What should the name for '", i, "' be: "))
        datum[datum$real_name == i,]$real_name <- iter
      } else {
        iter <- i
      }
      the_apes <- str_replace_all(datum[datum$real_name == iter,]$PE, "\\(", "")
      the_apes <- str_split(str_replace_all(the_apes, "\\)", "")," ")[[1]]
      the_apes <- round(as.numeric(the_apes[2:length(the_apes)]),4)
      if (!exists("all_datum")){
        all_datum <- tibble("Name" = iter,
                            "Position" = c(1:length(the_apes)),
                            "PE" = the_apes)
      } else {
        interim <- tibble("Name" = iter,
                          "Position" = c(1:length(the_apes)),
                          "PE" = the_apes)
        all_datum <- bind_rows(all_datum, interim)
      }
    }
    print(unique(all_datum$Name))
    draft <- ggplot(data = all_datum, aes(x = Position, y = PE, color=PE, fill = PE))+
      geom_bar(stat = "identity")+
      scale_color_viridis_c()+
      scale_fill_viridis_c()+
      theme_bw()+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid = element_blank())+
      facet_wrap(~Name, ncol=1)
    print(draft)
  }
}
