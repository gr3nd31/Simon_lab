library(ggpubr)
library(tidyverse)

fault_line <- function(genome = "genome_df.csv",
                       reads = "reads_df.csv",
                       graph_images = T,
                       save_images = T,
                       max_percent = 10,
                       max_size = F){
  genome_df <- read_csv(genome)
  reads_df <- read_csv(reads)
  
  if (max_size){
    maxie <- readline(prompt = paste0("What max size should be applied? (genome is ", nrow(genome_df), "): "))
  }
  reads_df <- reads_df[reads_df$aligned_length <= as.numeric(maxie),]
  
  genome_df$plus_depth <- 0
  genome_df$minus_depth <- 0
  genome_df$plus_starts <- 0
  genome_df$plus_ends <- 0
  genome_df$minus_starts <- 0
  genome_df$minus_ends <- 0
  
  for (i in c(1:nrow(genome_df))){
    for (j in unique(reads_df$h_strand)){
      crossers <- nrow(reads_df[reads_df$h_start <= i &
                                  reads_df$h_end >= i &
                                  reads_df$h_strand == j,])
      starters <- nrow(reads_df[reads_df$h_start ==i &
                                  reads_df$h_strand == j,])
      enders <- nrow(reads_df[reads_df$h_end ==i &
                                reads_df$h_strand == j,])
      if (j == "Plus"){
        genome_df[genome_df$Position_num == i,]$plus_depth <- crossers
        genome_df[genome_df$Position_num == i,]$plus_starts <- starters
        genome_df[genome_df$Position_num == i,]$plus_ends <- enders
      } else {
        genome_df[genome_df$Position_num == i,]$minus_depth <- crossers
        genome_df[genome_df$Position_num == i,]$minus_starts <- starters
        genome_df[genome_df$Position_num == i,]$minus_ends <- enders
      }
    }
  }
  genome_df$plus_starts_percent <- 100*genome_df$plus_starts/genome_df$plus_depth
  genome_df$plus_ends_percent <- 100*genome_df$plus_ends/genome_df$plus_depth
  genome_df$minus_starts_percent <- 100*genome_df$minus_starts/genome_df$minus_depth
  genome_df$minus_ends_percent <- 100*genome_df$minus_ends/genome_df$minus_depth
  
  genome_df$total_depth <- genome_df$plus_depth+genome_df$minus_depth
  genome_df$total_starts <- genome_df$plus_starts+genome_df$minus_starts
  genome_df$total_ends <- genome_df$plus_ends+genome_df$minus_ends
  genome_df$starts_percent <- 100*genome_df$total_starts/genome_df$total_depth
  genome_df$ends_percent <- 100*genome_df$total_ends/genome_df$total_depth
  
  write_csv(genome_df, genome)
  if (graph_images){
    draft <- ggbarplot(data = genome_df, 
                       x = "Position_num",
                       y = "starts_percent")+
      ylim(0,max_percent)+
      labs(x = "Genome position",
           y = "Start frequency")+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    print(draft)
    if (save_images){
      ggsave("start_percents.png", dpi = 300)
    }
    draft <- ggbarplot(data = genome_df, 
                       x = "Position_num",
                       y = "ends_percent")+
      ylim(0,max_percent)+
      labs(x = "Genome position",
           y = "End frequency")+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    print(draft)
    if (save_images){
      ggsave("end_percents.png", dpi = 300)
    }
  }
}
