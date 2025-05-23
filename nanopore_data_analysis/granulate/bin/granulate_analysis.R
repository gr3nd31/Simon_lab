#Contains the following functions:
#-------------
# graph_coverage_map()
# Analysis of the genome_df.csv file
#-------------
# graph_reads_map()
# Analysis of the reads_df.csv file
#-------------
# foldback_finder()
#  Analyzes reads_df.csv dataframe to find potential foldback seqeunces
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

suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(cowplot))

#-------------------------------
# graph_coverage_map()
#   Analyzes and summarizes error data and generates plots
#-------------------------------
graph_coverage_map <- function(file_name="genome_df.csv",
                       window_length = 2,
                       graph_them = T,
                       save_them=TRUE,
                       log_transform = F,
                       graph_skew = FALSE,
                       prefix=""){
  datum <- read_csv(file_name, show_col_types = F)
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
    print("Coverage map data already detected. Skipping analysis")
  }
  
  if(log_transform){
    datum$number_hit <- log(datum$number_hit, 10)
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
    if (graph_skew){
      datum$Position_num_2 <- as.character(datum$Position_num)
      draft <- ggplot(data = datum, 
                      aes(x = Position_num, y=Freq,
                          fill = Position_num,
                          color = Position_num))+
        geom_bar(stat = "identity")+
        scale_fill_viridis_b(option="viridis")+
        scale_color_viridis_b(option="viridis")+
        theme_bw()+
        labs(x = "Alignment Position", y = "Frequency")+
        theme(legend.position = "top",
              axis.text.x = element_blank(), 
              axis.ticks.x = element_blank())
      print(draft)
      if(save_them){
        ggsave(paste0(prefix, "coverage_binned.png"), dpi=500)
      }
      
      draft <- ggplot(data = datum, 
                      aes(x = fct_reorder(Position_num_2, -Freq), y=Freq,
                          fill = Position_num,
                          color = Position_num))+
        geom_bar(stat = "identity")+
        scale_fill_viridis_b(option="viridis")+
        scale_color_viridis_b(option="viridis")+
        theme_bw()+
        labs(x = "Alignment Position", y = "Frequency")+
        theme(legend.position = "top",
              axis.text.x = element_blank(), 
              axis.ticks.x = element_blank())
      print(draft)
      if(save_them){
        ggsave(paste0(prefix, "coverage_sorted.png"), dpi=500)
      }
    }
  }
}

#-------------------------------
# graph_reads_map()
#  Analyzes read data to annotate multi-alignments and generate
#  read alignment plots
#-------------------------------
graph_reads_map <- function(file_name = "reads_df.csv",
                         subset_num = 1000,
                         save_subset = F,
                         split_hsps=F,
                         split_strand = F,
                         color_strand = F,
                         set_line = 0,
                         log_transform = F,
                         graph_it = T,
                         save_it = T,
                         prefix = "",
                         alph = 0.2,
                         width = 1,
                         pre_sub = F,
                         make_widget = F,
                         skip_skips = F,
                         genome_size = 0,
                         min_set = 0.1,
                         max_set = 0.85){
  reanalyze <- T
  run_it <- T
  
  # Read the file
  if (typeof(file_name) == "list"){
    datum <- file_name
    reanalyze <- F
  } else if (typeof(file_name) == "character"){
    datum <- read_csv(file_name, show_col_types = F)
  } else{
    print("Unable to read input dataframe. Aborting.")
    run_it <- F
  }
  
  if (subset_num > length(unique(datum$read_id))){
    print(paste0("Requested subset number (", subset_num, ") is greater than existing reads (", length(unique(datum$read_id)), "). Defaulting to all reads."))
    subset_num <- length(unique(datum$read_id))
  }
  
  if (genome_size == 0){
    genome_size <- max(datum$length)
  }
  
  if (run_it & reanalyze){
    lcounter <- 0
    # Skips aligned_length generation if already present
    if (!"aligned_length" %in% names(datum)){
      if (pre_sub & subset_num > 0){
        print(paste0("Pre-subsetting to ", subset_num, " reads."))
        # Randomly samples reads
        subset_list <- sample(datum$read_id, subset_num)
        datum <- datum[datum$read_id %in% subset_list,]
      }
      # Sums the length of all alignments found on the read
      datum$aligned_length <- 0
      datum$hsps_count <- 1
      datum$north_mer <- "Fragment"
      datum$seq_mer <- "Fragment"
      datum$sense_mer <- "NONE"
      
      for (i in unique(datum$read_id)){
        lcounter <- lcounter+1
        if(lcounter%%5000 == 0){
          print(paste0("Analyzing read: ", lcounter, " (", round(lcounter/length(unique(datum$read_id)),2)*100, "%)"))
        }
        # Counts the number of HSPSs on the read
        datum[datum$read_id ==i,]$hsps_count <- nrow(datum[datum$read_id==i,])
        # Adds the total aligned length on the read
        datum[datum$read_id ==i,]$aligned_length <- sum(datum[datum$read_id==i,]$length)
        
        # Defines the read as plus only, minus only, or mixed
        interim <- datum[datum$read_id == i,]
        interim$length <- abs(interim$h_end - interim$h_start)
        if (length(unique(interim$h_strand))==1){
          if (unique(interim$h_strand) == "Plus"){
            datum[datum$read_id ==i,]$sense_mer <- "Plus"
          } else if (unique(interim$h_strand) == "Minus") {
            datum[datum$read_id ==i,]$sense_mer <- "Minus"
          }
        } else if (length(unique(interim$h_strand)) == 2){
          datum[datum$read_id ==i,]$sense_mer <- "Mixed"
        }
        
        # gets the number times the genome size fits in the aligned length of the read
        n_check <- sum(interim$length)/genome_size
        # If the read is large enough to not be a fragment, its mer-classified
        if (n_check >= max_set){
          # Gets the remainder of the length of the read
          m_check <- (sum(interim$length)%%genome_size)/genome_size
          # if the remainder is long enough to count as a full length, its added as a full length
          if (m_check >= max_set) {
            datum[datum$read_id == i,]$north_mer <- paste0(as.character(sum(interim$length)%/%genome_size+1), "-mer")
          } else if (m_check <= min_set){ #if its less than a certain amount, its ignored
            datum[datum$read_id == i,]$north_mer <- paste0(as.character(sum(interim$length)%/%genome_size), "-mer")
          } else { # if its between, its a +-mer
            datum[datum$read_id == i,]$north_mer <- paste0(as.character(sum(interim$length)%/%genome_size), "+-mer")
          }
        }
        
        # counts the number of HSPSs that are longer than the minimal length to be considered 'full'
        s_check <- as.numeric(nrow(interim[interim$length >= genome_size*max_set,]))
        # if its longer than 1, its mer'd
        if (s_check >= 1){
          # counts the number of HSPSs that are longer than a fragment but shorter than full length
          m_check <- nrow(interim[interim$length > genome_size*min_set & interim$length < genome_size*max_set,])
          # intermediate fragments are detected, read is a +-mer
          if (m_check > 0){
            datum[datum$read_id == i,]$seq_mer <- paste0(as.character(s_check), "+-mer")
          } else {
            datum[datum$read_id == i,]$seq_mer <- paste0(as.character(s_check), "-mer")
          }
        }
      }
      datum$deletion_ave <- datum$deletion_num/datum$length
      datum$mismatch_ave <- datum$mismatch_num/datum$length
      datum$insert_ave <- datum$insert_num/datum$length
      #Orders by negative length
      datum <- datum[order(-datum$aligned_length),]
      datum$a_id <- ordered(-datum$aligned_length)
      # Saves the file
      if (pre_sub){
        write_csv(datum, paste0("pre_subbed_to_", subset_num, "_",file_name))
      } else {
        write_csv(datum, file_name)
      }
    } else {
      # Orders by negative length in case of graphing
      print("Aggregate aligment length detected. Skipping analysis.")
      subset_list <- sample(datum$read_id, subset_num)
      datum <- datum[datum$read_id %in% subset_list,]
      
      datum <- datum[order(-datum$aligned_length),]
      datum$a_id <- ordered(-datum$aligned_length)
    }
  } else {
    # Orders by negative length in case of graphing
    print("Aggregate aligment length detected. Skipping analysis.")
    subset_list <- sample(datum$read_id, subset_num)
    datum <- datum[datum$read_id %in% subset_list,]
    datum <- datum[order(-datum$aligned_length),]
    datum$a_id <- ordered(-datum$aligned_length)
  }
  
  if (skip_skips){
    datum <- datum[datum$skipped == "FALSE",]
  }
  if (log_transform){
    datum$a_id <- -log(datum$aligned_length, 10)
  }
  # Graphs data
  if(graph_it){
    # Graphs the data as a stacked-alignment plot
    if (color_strand){
      draft <- ggplot(datum)+
        geom_segment(aes(x=h_start, xend=h_end, 
                         y=a_id, yend=a_id,
                         color = h_strand), alpha=alph, linewidth = width)+
        ylab("Alignment length")+
        xlab("Genome")+
        theme_bw()+
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid = element_blank())
    } else {
      draft <- ggplot(datum)+
        geom_segment(aes(x=h_start, xend=h_end, 
                         y=a_id, yend=a_id), alpha=alph, linewidth = width)+
        ylab("Alignment length")+
        xlab("Genome")+
        theme_bw()+
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid = element_blank())
    }
    
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
    if (make_widget){
      htmlwidgets::saveWidget(ggplotly(draft), paste0(prefix, "stacks.html"), selfcontained = T)
    }
    if(save_subset){
      write_csv(datum, "subsetReads_df.csv")
      subsetReads <<- datum
      print("Saved subset as variable 'subsetReads'")
    }
    
    interim <- as.data.frame(table(datum$north_mer))
    interim$Freq <- interim$Freq/sum(interim$Freq)
    n_graph <- ggplot(data = interim, aes(x=Var1, y = Freq))+
      geom_bar(stat = "identity")+
      labs(x = "Read Type", y = "Frequency", title = "Northern definitions")+
      theme_cowplot(12)
    sinterim <- as.data.frame(table(datum$seq_mer))
    sinterim$Freq <- sinterim$Freq/sum(sinterim$Freq)
    s_graph <- ggplot(data = sinterim, aes(x=Var1, y = Freq))+
      geom_bar(stat = "identity")+
      labs(x = "Read Type", y = "Frequency", title = "Sequence definitions")+
      theme_cowplot(12)
    print(plot_grid(n_graph, s_graph, ncol = 1))
    ggsave(paste0(prefix,"mer_graph.png"), height = 8, width = 7, dpi = 300)
  }
}

#-------------------------------
# graph_read_hsps
#  graphs the position of hspss within each read
#-------------------------------

graph_hsps <- function(file_name = "reads_df.csv",
                                as_percent = F,
                                color_by = "number",
                                save_it = T,
                                save_as = "read_hsps.png",
                                subset_data = T,
                                save_subset = F,
                                subset_number = 20,
                                skips = T,
                                plotly_it = F){
  # Checks to see if the file name is a dataframe or a file name
  if (typeof(file_name)=="character"){
    datum <- read_csv(file_name, show_col_types = F)
  } else if (typeof(file_name)=="list"){
    datum <- file_name
  }
  
  runIt <- TRUE
  
  # Subsets the data
  if (subset_data & typeof(subset_number) == "double"){
    print("Sampling data...")
    read_ids <- unique(datum$read_id)
    if (length(read_ids) > subset_number){
      read_ids <- sample(read_ids, subset_number)
      datum <- datum[datum$read_id %in% read_ids,]
    } else {
      print(paste0("Only ", length(read_ids), " reads found but ", subset_number, " requested. Returning total reads."))
    }
  }
  if (as_percent){
    print("Adjusting reads as percent.")
  }
  if (!skips){
    datum <- datum[datum$skipped == "False",]
  }
  
  if (!"aligned_length" %in% names(datum)){
    print("Reads have not been mapped yet. Please run graph_reads_map() first.")
    runIt <- FALSE
  }
  
  if (runIt){
    # Orders the dataframe and applies levels to the hit strand sense
    datum <- datum[order(-datum$aligned_length),]
    datum$h_strand <- ordered(datum$h_strand, levels = c("Plus", "Minus"))
    datum$hsps <- as.character(datum$hsps)
    # Creates a second dataframe to record the total read information and assigns the same y-axis number for HSPS within those reads
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
    # Changes the end number if the data is to be presented as a percent
    if (as_percent){
      interim$end <- 100
    }
    
    # Graphs data with various coloring schemes
    draft <- ggplot()+
      geom_segment(data = interim, aes(x = start, y = read_id, xend = end, yend = read_id),
                   color = "darkgray", linewidth = 5, alpha = 1, lineend = "round")
    if (color_by == "sense"){
      draft <- draft+
        geom_segment(data = datum, aes(x = q_start, y = read_id, xend = q_end, yend = read_id, color = h_strand),
                     linewidth = 4, alpha = 1, lineend = "round")+
        scale_color_manual(values = c("orange", "lightblue"))+
        labs(color = "Alignment Sense")
    } else if (color_by == "start"){
      draft <- draft+
        geom_segment(data = datum, aes(x = q_start, y = read_id, xend = q_end, yend = read_id, color = h_start),
                     linewidth = 4, alpha = 1, lineend = "round")+
        scale_color_viridis_b(option = "inferno")+
        labs(color = "Alignment Start")
    } else if (color_by == "end"){
      draft <- draft+
        geom_segment(data = datum, aes(x = q_start, y = read_id, xend = q_end, yend = read_id, color = h_end),
                     linewidth = 4, alpha = 1, lineend = "round")+
        scale_color_viridis_b(option = "inferno")+
        labs(color = "Alignment End")
    } else if (color_by == "mismatch"){
      draft <- draft+
        geom_segment(data = datum, aes(x = q_start, y = read_id, xend = q_end, yend = read_id, color = mismatch_ave*100),
                     linewidth = 4, alpha = 1, lineend = "round")+
        scale_color_viridis_b(option = "inferno")+
        labs(color = "Mismatch Freq")
    } else if (color_by == "deletion"){
      draft <- draft+
        geom_segment(data = datum, aes(x = q_start, y = read_id, xend = q_end, yend = read_id, color = deletion_ave*100),
                     linewidth = 4, alpha = 1, lineend = "round")+
        scale_color_viridis_b(option = "inferno")+
        labs(color = "Deletion Freq")
    } else if (color_by == "insertion"){
      draft <- draft+
        geom_segment(data = datum, aes(x = q_start, y = read_id, xend = q_end, yend = read_id, color = insert_ave*100),
                     linewidth = 4, alpha = 1, lineend = "round")+
        scale_color_viridis_b(option = "inferno")+
        labs(color = "Insertion Freq")
    } else if (color_by == "number"){
      draft <- draft+
        geom_segment(data = datum, aes(x = q_start, y = read_id, xend = q_end, yend = read_id, color = hsps),
                     linewidth = 4, alpha = 1, lineend = "round")+
        labs(color = "HSPS number")+
        theme(legend.position = element_blank())
    } else if (color_by == "reference"){
      draft <- draft+
        geom_segment(data = datum, aes(x = q_start, y = read_id, xend = q_end, yend = read_id, color = reference),
                     linewidth = 4, alpha = 1, lineend = "round")+
        labs(color = "Reference")
    }else {
      print("Coloring parameter not recognized. Please try with 'number', 'sense', 'start', 'end', 'mismatch', 'deletion', 'insertion', or 'reference'. Defaulting to 'number'.")
      draft <- draft+
        geom_segment(data = datum, aes(x = q_start, y = read_id, xend = q_end, yend = read_id, color = hsps),
                     linewidth = 4, alpha = 1, lineend = "round")+
        labs(color = "HSPS Number")+
        theme(legend.position = element_blank())
    }
    draft <- draft+
      theme_bw()+
      theme(axis.text.y = element_blank(),
            panel.grid.major.y = element_blank(), 
            axis.ticks.y = element_blank(),
            legend.position = "top",
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 14))
    if (as_percent){
      draft <- draft+
        labs(x = "Relative Read Position (%)",
             y = "")
    } else {
      draft <- draft+
        labs(x = "Read Position (nt)",
             y = "")
    }
    
    # Prints graph and saves
    print(draft)
    #x <<- draft
    if (save_it & nchar(save_as) > 0){
      ggsave(save_as)
    }
    if(save_subset){
      write_csv(datum, "subsetReads_HSPS_df.csv")
      subsetHSPS <<- datum
      print("Saved subset as variable 'subsetHSPS'")
    }
  }
}

#-------------------------------
# graph_northern
#  graphs reads as they would appear in a northern blot
#  using probes designated in the parameters
#-------------------------------

graph_northern <- function(fileName = "reads_df.csv",
                           probe_start = 30,
                           probe_end = 130,
                           strands = "Plus",
                           color_by = "none",
                           max_size = 5000,
                           min_size = 50,
                           subset_number = F,
                           tag = "",
                           transform_length = F,
                           transform_number = 2,
                           ylims = FALSE,
                           line_wd = 1,
                           alph = 0.01, 
                           save_it = T){
  runIt <- TRUE
  if(typeof(fileName) == "character"){
    datum <- read_csv(fileName, show_col_types = F)
  } else if (typeof(fileName) == "list"){
    datum <- fileName
  } else {
    print("Unable to open or read data. Aborting.")
    runIt <- FALSE
  }
  if (runIt){
    if (strands == "Plus" & probe_start > 0 & probe_end > 0){
      interim <- datum[datum$h_start <= probe_start &
                         datum$h_end >= probe_end & 
                         datum$h_strand == strands &
                         datum$read_length <= max_size &
                         datum$read_length >= min_size,]
    } else if (strands == "Minus" & probe_start > 0 & probe_end > 0){
      interim <- datum[datum$h_start >= probe_end &
                         datum$h_end <= probe_start & 
                         datum$h_strand == strands &
                         datum$read_length <= max_size &
                         datum$read_length >= min_size,]
      if(nchar(tag) > 0){
        tag <- paste0(tag, "_Minus")
      } else{
        tag <- paste0(tag, "Minus")
      }
    } else if (probe_start == 0 & probe_end == 0){
      interim <- datum[datum$read_length <= max_size &
                         datum$read_length >= min_size,]
    } else {
      print("Unable to reckonize strand request, returning 'Plus'")
      interim <- datum[datum$h_start <= probe_start &
                         datum$h_end >= probe_end & 
                         datum$h_strand == strands &
                         datum$read_length <= max_size &
                         datum$read_length >= min_size,]
    }
    
    if (transform_length){
      interim$read_length <- log(interim$read_length, transform_number)
    }
    if (typeof(subset_number) == "double"){
      if (length(unique(interim$read_id)) > subset_number){
        interim <- interim[sample(nrow(interim), subset_number),]
      } else {
        print(paste0("Only ", length(unique(interim$read_id)), " reads matching the parameters. Only graphing that number"))
      }
    }
    draft <- ggplot(data = interim)
    if (color_by == "none"){
      draft <- draft+
        geom_segment(aes(x=0.5, xend = 1.5, y=read_length, yend = read_length), alpha = alph, linewidth = line_wd)
    } else if(color_by == "start"){
      draft <- draft+
        geom_segment(aes(x=0.5, xend = 1.5, y=read_length, yend = read_length, color = h_start), alpha = alph, linewidth = line_wd)+
        scale_color_viridis_c(option = "turbo")
    } else if(color_by == "end"){
      draft <- draft+
        geom_segment(aes(x=0.5, xend = 1.5, y=read_length, yend = read_length, color = h_end), alpha = alph, linewidth = line_wd)+
        scale_color_viridis_c(option = "turbo")
    } else if(color_by=="hsps"){
      draft <- draft+
        geom_segment(aes(x=0.5, xend = 1.5, y=read_length, yend = read_length, color = hsps_count), alpha = alph, linewidth = line_wd)+
        scale_color_viridis_b(option = "turbo")
    } else if(color_by == "sense"){
      draft <- draft+
        geom_segment(aes(x=0.5, xend = 1.5, y=read_length, yend = read_length, color = h_strand), alpha = alph, linewidth = line_wd)+
        scale_color_manual(values = c("darkred", 'darkblue'))
    } else {
      print("Unable to color by given argument and defaulting to 'none'. Valid argumnets are 'start', 'end', 'sene', or 'hsps'.")
      draft <- draft+
        geom_segment(aes(x=0.5, xend = 1.5, y=read_length, yend = read_length), alpha = alph, linewidth = line_wd)
    }
    draft <- draft+
      theme_bw()+
      labs(y = "Read Length")+
      theme(panel.grid = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 10),
            axis.text.y = element_text(size = 8, angle = 90, vjust = 1, hjust = 1),
            legend.position = "top")
    if (typeof(ylims) == "double"){
      draft <- draft+ylim(ylims[1], ylims[2])
    }
    print(draft)
    if (nchar(tag) > 0){
      tag <- paste0("_", tag)
    }
    if (save_it){
      ggsave(paste0("Northern_probe_",probe_start, "-", probe_end ,tag,".png"), dpi = 300, width = 1, height = 5)
    }
  }
}

#-------------------------------
# foldback_finder
#  Analyzes reads_df.csv dataframe to find potential foldback seqeunces
#-------------------------------

foldback_finder <- function(df,
                            checkForHalves = TRUE,
                            minRatio = 1.7,
                            maxRatio = 2.1,
                            minSize = 500,
                            graph_reads = TRUE,
                            outList = "PossibleFBs.txt",
                            outCSV = "PossibleFBs.csv"){
  runIt <- TRUE
  if (typeof(df) == "character"){
    datum <- read_csv(df)
  } else if (typeof(df) == "list"){
    datum <- df
  } else {
    print("Unable to read given reads list Abortin'")
    runIt <- FALSE
  }
  if (runIt){
    datum <- datum[datum$aligned_length > minSize,]
    theList <- ""
    for (i in unique(datum$read_id)){
      if (length(unique(datum[datum$read_id == i,]$h_strand)) > 1){
        theList <- paste0(theList, i,"\n")
        if (!exists("FBs")){
          FBs <- datum[datum$read_id == i,]
        } else {
          FBs <- rbind(FBs, datum[datum$read_id == i,])
        }
      } else if (checkForHalves){
        readLength <- unique(datum[datum$read_id ==i,]$read_length)[1]
        alignedLength <- unique(datum[datum$read_id ==i,]$aligned_length)[1]
        if (readLength/alignedLength >= minRatio &
            readLength/alignedLength <= maxRatio){
          theList <- paste0(theList, i,"\n")
          if (!exists("FBs")){
            FBs <- datum[datum$read_id == i,]
          } else {
            FBs <- rbind(FBs, datum[datum$read_id == i,])
          }
        }
      }
    }
    
    if (nchar(theList) > 0){
      cat(theList, file = outList)
      write_csv(FBs, outCSV)
    }
    if (graph_reads & exists("graph_read_hsps") & nchar(theList) > 0){
      graph_read_hsps(FBs, color_by = "sense", as_percent = T)
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
  interim <- df[unlist(df[rel_percent_missed]) > threshold & 
                  !is.na(unlist(df[rel_percent_missed])),]
  #print(interim)
  for (i in 1:nrow(interim)){
    billies <- unlist(interim[i,list_of_misses])
    the_pick <- strsplit(billies, "_")[[1]]
    the_g <- as.numeric(strsplit(the_pick[1], ":")[[1]][2])
    the_c <- as.numeric(strsplit(the_pick[2], ":")[[1]][2])
    the_a <- as.numeric(strsplit(the_pick[3], ":")[[1]][2])
    the_t <- as.numeric(strsplit(the_pick[4], ":")[[1]][2])
    the_u <- as.numeric(strsplit(the_pick[5], ":")[[1]][2])
    totes <- the_g+the_c+the_a+the_t+the_u
    new_df <- tibble("Position" = interim[i,]$Position_num,
                     "Identity" = interim[i,]$Position_nt,
                     "Count" = totes,
                     "A"=100*the_a/totes,
                     "T"=100*(the_t+the_u)/totes,
                     "G"=100*the_g/totes,
                     "C"=100*the_c/totes)
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
  nt_data <- nt_data %>% replace(is.na(.), 0)
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
                       prefix = "",
                       graph_images = T,
                       save_images = T,
                       max_percent = 100,
                       min_size = F,
                       max_size = F,
                       quant = 0.9){
  genome_df <- read_csv(genome)
  reads_df <- read_csv(reads)
  
  print(paste0("Total reads loaded: ", nrow(reads_df)))
  if (max_size){
    maxie <- readline(prompt = paste0("What max size should be applied? (genome is ", nrow(genome_df), "): "))
  } else {
    maxie <- nrow(genome_df)
  }
  if (min_size){
    minnie <- readline(prompt = paste0("What min size should be applied? (Minumum alignment is ", min(reads_df$aligned_length), "): "))
  } else {
    minnie <- 1
  }
  reads_df <- reads_df[reads_df$aligned_length <= as.numeric(maxie) &
                         reads_df$aligned_length >= as.numeric(minnie),]
  if (max_size | min_size){
    print(paste0("Reads after size trimming: ", nrow(reads_df)))
  }
  
  genome_df$plus_depth <- 0
  genome_df$minus_depth <- 0
  genome_df$plus_starts <- 0
  genome_df$plus_ends <- 0
  genome_df$minus_starts <- 0
  genome_df$minus_ends <- 0
  
  for (i in c(1:nrow(genome_df))){
    for (j in unique(reads_df$h_strand)){
      if (j == "Plus"){
        crossers <- nrow(reads_df[reads_df$h_start <= i &
                                    reads_df$h_end >= i &
                                    reads_df$h_strand == j,])
        starters <- nrow(reads_df[reads_df$h_start ==i &
                                    reads_df$h_strand == j,])
        enders <- nrow(reads_df[reads_df$h_end ==i &
                                  reads_df$h_strand == j,])
        genome_df[genome_df$Position_num == i,]$plus_depth <- crossers
        genome_df[genome_df$Position_num == i,]$plus_starts <- starters
        genome_df[genome_df$Position_num == i,]$plus_ends <- enders
      } else {
        crossers <- nrow(reads_df[reads_df$h_start >= i &
                                    reads_df$h_end <= i &
                                    reads_df$h_strand == j,])
        starters <- nrow(reads_df[reads_df$h_end == i &
                                    reads_df$h_strand == j,])
        enders <- nrow(reads_df[reads_df$h_start ==i &
                                  reads_df$h_strand == j,])
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
  genome_df$total_starts <- genome_df$plus_ends+genome_df$minus_starts
  genome_df$total_ends <- genome_df$plus_starts+genome_df$minus_ends
  genome_df$starts_percent <- 100*genome_df$total_starts/genome_df$total_depth
  genome_df$ends_percent <- 100*genome_df$total_ends/genome_df$total_depth
  genome_df$ends_quantile <- genome_df$ends_percent/quantile(genome_df[!is.na(genome_df$ends_percent),]$ends_percent, probs = seq(0,1,quant))[2]
  
  #print(names(genome_df))
  trixie <- genome_df[,c(1:3,15:31)]
  
  write_csv(trixie, paste0(prefix,"fault_lines.csv"))
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
      ggsave(paste0(prefix, "start_percents.png"), dpi = 300)
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
      ggsave(paste0(prefix, "end_percents.png"), dpi = 300)
    }
    draft <- ggbarplot(data = genome_df, 
                       x = "Position_num",
                       y = "ends_quantile")+
      labs(x = "Genome position",
           y = "Relative end frequency")+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    print(draft)
    if (save_images){
      ggsave("end_quant.png", dpi = 300)
    }
  }
}

fault_align <- function(shaped = "fault_lines.csv",
                        control = "dmso_fault_lines.csv",
                        prefix = "",
                        quant = 0.9,
                        save_images = T){
  shape_data <- read_csv(shaped)
  shape_name <- strsplit(shaped, ".csv")[[1]][1]
  print(shape_name)
  control_data <- read_csv(control)
  control_name <- strsplit(control, ".csv")[[1]][1]
  
  shape_data$relative_ends <- shape_data$ends_percent-control_data$ends_percent
  shape_data$relative_ends_quantile <- shape_data$relative_ends/quantile(shape_data[!is.na(shape_data$relative_ends),]$relative_ends, probs = seq(0,1,quant))[2]
  draft <- ggbarplot(data = shape_data, 
                     x = "Position_num",
                     y = "relative_ends_quantile")+
    ylim(-1, 5)+
    labs(x = "Genome position",
         y = "Relative end frequency")+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  print(draft)
  if (save_images){
    ggsave("rel_end_quant.png", dpi = 300)
  }
  write_csv(shape_data, paste0(shape_name, "_normalized_to_", prefix, ".csv"))
}

