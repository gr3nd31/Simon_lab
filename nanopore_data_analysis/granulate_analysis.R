#Contains the following functions:
#-------------
# granugraph()
# Analysis of the genome_df.csv file
#-------------
#stock_o_late()
# Analysis of the reads_df.csv file
#-------------
# pull_ids()
# Pulls read ids from a given df and saves to 'reads.txt'
#------------
#de_missed()
# Analyses the mismatch data from the genome_df.csv file.
#-----------
# fault_line
# Compares the read_df and genome_df to calculate the starting
# and ending frequency for reads
#-----------

library(ggpubr)
library(tidyverse)

#-------------------------------
# granugraph()
#   Analyzes and summarizes error data and generates plots
#-------------------------------
granugraph <- function(file_name="genome_df.csv",
                       window_length = 2,
                       graph_them = T,
                       save_them=TRUE,
                       prefix=""){
  datum <- read_csv(file_name)
  if (!"plume_deleted" %in% names(datum) &
      !"plume_mismatched" %in% names(datum) &
      !"plume_inserted" %in% names(datum)){
    datum$rel_deleted <- datum$number_deleted/datum$number_hit
    datum$rel_mismatched <- datum$number_mismatched/datum$number_hit
    datum$rel_inserted <- datum$number_inserted/datum$number_hit
    datum$plume_deleted <- 0
    datum$plume_mismatched <- 0
    datum$plume_inserted <- 0
    datum$Freq <- datum$number_hit/max(datum$number_hit)
    for (i in 1:nrow(datum)){
      min_num <- i-window_length
      max_num <- i+window_length
      if (min_num<1) {
        min_num<-1
      }
      if (max_num > nrow(datum)){
        max_num <- nrow(datum)
      }
      datum[i,"plume_deleted"] <- mean(datum[min_num:max_num,]$rel_deleted)
      datum[i,"plume_mismatched"] <- mean(datum[min_num:max_num,]$rel_mismatched)
      datum[i,"plume_inserted"] <- mean(datum[min_num:max_num,]$rel_inserted)
    }
    write_csv(datum, file_name)
  } else {
    print("Granugraph data already detected. Skipping analysis")
  }
  
  if (graph_them){
    draft <- ggplot(data = datum, 
                    aes(x = Position_num, y=Freq))+
      geom_bar(stat = "identity", color = "black", fill = "black")+
      labs(x = "Alignment Position", y = "Coverage")+
      theme_bw()+
      theme(legend.position = "top")
    print(draft)
    if(save_them){
      ggsave(paste0(prefix, "coverage.png"), dpi=500)
    }
    
    draft <- ggplot(data = datum, 
                    aes(x = Position_num, y=number_hit,
                        fill = rel_mismatched,
                        color = rel_mismatched))+
      geom_bar(stat = "identity")+
      scale_fill_viridis_c(option="inferno")+
      scale_color_viridis_c(option="inferno")+
      labs(x = "Alignment Position", y = "Coverage")+
      theme_bw()+
      theme(legend.position = "top")
    print(draft)
    if(save_them){
      ggsave(paste0(prefix, "rel_mismatched.png"), dpi=500)
    }
    draft <- ggplot(data = datum, 
                    aes(x = Position_num, y=number_hit,
                        fill = rel_deleted,
                        color = rel_deleted))+
      geom_bar(stat = "identity")+
      scale_fill_viridis_c(option="inferno")+
      scale_color_viridis_c(option="inferno")+
      theme_bw()+
      labs(x = "Alignment Position", y = "Coverage")+
      theme(legend.position = "top")
    print(draft)
    if(save_them){
      ggsave(paste0(prefix, "rel_deleted.png"), dpi=500)
    }
    draft <- ggplot(data = datum, 
                    aes(x = Position_num, y=number_hit,
                        fill = rel_inserted,
                        color = rel_inserted))+
      geom_bar(stat = "identity")+
      scale_fill_viridis_c(option="inferno")+
      scale_color_viridis_c(option="inferno")+
      theme_bw()+
      labs(x = "Alignment Position", y = "Coverage")+
      theme(legend.position = "top")
    print(draft)
    if(save_them){
      ggsave(paste0(prefix, "rel_inserted.png"), dpi=500)
    }
    
    draft <- ggplot(data = datum, 
                    aes(x = Position_num, y=number_hit,
                        fill = plume_mismatched,
                        color = plume_mismatched))+
      geom_bar(stat = "identity")+
      scale_fill_viridis_c(option="inferno")+
      scale_color_viridis_c(option="inferno")+
      theme_bw()+
      labs(x = "Alignment Position", y = "Coverage")+
      theme(legend.position = "top")
    print(draft)
    if(save_them){
      ggsave(paste0(prefix, "plume_mismatched.png"), dpi=500)
    }
    draft <- ggplot(data = datum, 
                    aes(x = Position_num, y=number_hit,
                        fill = plume_deleted,
                        color = plume_deleted))+
      geom_bar(stat = "identity")+
      scale_fill_viridis_c(option="inferno")+
      scale_color_viridis_c(option="inferno")+
      theme_bw()+
      labs(x = "Alignment Position", y = "Coverage")+
      theme(legend.position = "top")
    print(draft)
    if(save_them){
      ggsave(paste0(prefix, "plume_deleted.png"), dpi=500)
    }
    draft <- ggplot(data = datum, 
                    aes(x = Position_num, y=number_hit,
                        fill = plume_inserted,
                        color = plume_inserted))+
      geom_bar(stat = "identity")+
      scale_fill_viridis_c(option="inferno")+
      scale_color_viridis_c(option="inferno")+
      theme_bw()+
      labs(x = "Alignment Position", y = "Coverage")+
      theme(legend.position = "top")
    print(draft)
    if(save_them){
      ggsave(paste0(prefix, "plume_inserted.png"), dpi=500)
    }
  }
}

#-------------------------------
# stack_o_late()
#  Analyzes read data to annotate multi-alignments and generate
#  read alignment plots
#-------------------------------
stack_o_late <- function(file_name = "reads_df.csv",
                         subset_num = 0,
                         split_hsps=F,
                         split_strand = F,
                         set_line = 0,
                         graph_it = T,
                         save_it = T,
                         prefix = "",
                         alph = 0.2){
  # Read the file
  datum <- read_csv(file_name)
  # Skips aligned_length generation if already present
  if (!"aligned_length" %in% names(datum)){
    # Sums the length of all alignments found on the read
    datum$aligned_length <- 0
    datum$hsps_count <- 1
    for (i in unique(datum$read_id)){
      datum[datum$read_id ==i,]$hsps_count <- nrow(datum[datum$read_id==i,])
      datum[datum$read_id ==i,]$aligned_length <- sum(datum[datum$read_id==i,]$length)
    }
    datum$deletion_ave <- datum$deletion_num/datum$length
    datum$mismatch_ave <- datum$mismatch_num/datum$length
    datum$insert_ave <- datum$insert_num/datum$length
    #Orders by negative length
    datum <- datum[order(-datum$aligned_length),]
    datum$a_id <- ordered(-datum$aligned_length)
    # Saves the file
    write_csv(datum, file_name)
  } else {
    # Orders by negative length in case of graphing
    print("Aggregate aligment length detected. Skipping analysis.")
    datum <- datum[order(-datum$aligned_length),]
    datum$a_id <- ordered(-datum$aligned_length)
  }
  # Graphs is told to
  if(graph_it){
    #Subset the dataframe for very large files
    if (subset_num > 0){
      print(paste0("Subsetting graph to ", subset_num, " reads."))
      # Caps sampling at max of reads
      if (subset_num  > nrow(datum)){
        print("Subset given exceeds number of reads. Defaulting to max")
        subset_num <- nrow(datum)
      }
      # Randomly samples reads
      datum <- datum[sample(nrow(datum), subset_num),]
    }
    # Graphs the data as a stacked-alignment plot
    draft <- ggplot(datum)+
      geom_segment(aes(x=h_start, xend=h_end, 
                       y=a_id, yend=a_id), alpha=alph)+
      ylab("Alignment length")+
      xlab("Genome")+
      theme_bw()+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    # Adds facet_wraps
    if (split_hsps | split_strand){
      if (split_hsps & !split_strand){
        draft <- draft+facet_wrap(~hsps_count)
      } else if (!split_hsps & split_strand){
        draft <- draft+facet_wrap(~h_strand)
      } else {
        draft <- draft+facet_wrap(~h_strand*hsps_count)
      }
    }
    
    if (set_line !=0){
      draft <- draft+ geom_vline(xintercept = set_line, color = "red")
    }
    print(draft)
    if (save_it){
      ggsave(paste0(prefix, "stacks.png"), width = 8, height = 4, dpi = 500)
    }
  }
}

#-------------------------------
# pull_ids
#  Pulls read_ids from a dataframe
#-------------------------------
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

#-------------------------------
# de_missed()
#   Analyzes genome_df mismatch data to summarize mismatches
#-------------------------------
de_missed <- function(df = datum,
                      threshold = 0.2,
                      rel_percent_missed = "rel_mismatched",
                      list_of_misses="Other_nt", 
                      tag = ""){
  interim <- df[unlist(df[rel_percent_missed]) > threshold & !is.na(unlist(df[rel_percent_missed])),]
  for (i in 1:nrow(interim)){
    billies <- interim[i,list_of_misses]
    billies <- toupper(billies)
    new_df <- tibble("Position" = interim[i,]$Position_num,
                     "Identity" = interim[i,]$Position_nt,
                     "Count" = nchar(billies),
                     "A"=100*str_count(billies, "A")/nchar(billies),
                     "T"=100*str_count(billies, "T")/nchar(billies),
                     "G"=100*str_count(billies, "G")/nchar(billies),
                     "C"=100*str_count(billies, "C")/nchar(billies))
    if (!exists("gc_data")){
      gc_data <- new_df
    } else{
      gc_data <- rbind(gc_data, new_df)
    }
  }
  print(paste0("Found ", nrow(interim), " bases with the given threshold"))
  print(paste0("% changed to A: ", round(mean(gc_data$A),1)))
  print(paste0("% changed to T: ", round(mean(gc_data$T), 1)))
  print(paste0("% changed to G: ", round(mean(gc_data$G), 1)))
  print(paste0("% changed to C: ", round(mean(gc_data$C), 1)))
  #print(head(gc_data))
  #print(unique(gc_data$Identity))
  for (i in c("A", "T", "G", "C")){
    #print(i)
    interim <- gc_data[gc_data$Identity == i,]
    riff <- tibble("Identity" = i,
                   "Freq in sequence" = 100*nrow(df[df$Position_nt == i,])/nrow(df),
                   "Freq changed" = 100*sum(interim$Count)/sum(gc_data$Count),
                   "% to A" = mean(interim$A),
                   "% to T" = mean(interim$T),
                   "% to G" = mean(interim$G),
                   "% to C" = mean(interim$C),)
    if (!exists("nt_data")){
      nt_data <- riff
    } else {
      nt_data <- rbind(nt_data, riff)
    }
  }
  print(nt_data)
  nt_data <- pivot_longer(data = nt_data, cols = c(4:7), names_to = "Conversion", values_to = "Percent")
  draft <- ggbarplot(data = nt_data, x = "Identity", y = "Percent", fill = "Conversion",
                     position = position_dodge(0.8))+
    theme(legend.position = "right")
  print(draft)
  ggsave(paste0(tag, "de_missed.png"), dpi = 300, height = 2, width = 4)
  return(gc_data)
}

#-------------------------------
# fault_line
# Compares the read_df and genome_df to calculate the starting
# and ending frequency for reads
#-------------------------------

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

