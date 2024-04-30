library(tidyverse)
library(ggpubr)
library(stringi)

primer_spot <- function(file_name,
                        genome_file = "",
                        q_start = 1,
                        q_end = 2690,
                        primer_length = 20,
                        primer_align_threshold = 10,
                        min_amplicon_length = 100,
                        max_amplicon_length = 150,
                        number_of_amplicons = 30,
                        allowed_hits = 0,
                        match_primers = T,
                        map_primers = T){
  print("Pulling alignments...")
  #print(primer_length - primer_align_threshold)
  datum <- read_csv(file_name, show_col_types = F)
  names(datum) <- c("Query", 
                    "Seq_ID", 
                    "Percent_identity", 
                    "Length", 
                    "Gaps", 
                    "Gap_percent", 
                    "q_start", 
                    "q_end", 
                    "sub_start", 
                    "sub_end",
                    "Expect", 
                    "Score")
  
  if (genome_file != ""){
    print("Genome file detected...")
    the_genome <- unlist(read.table(genome_file, header = T, sep = "\n")[1,1])
  }
  datum$sub_Orientation <- "Plus"
  datum[datum$sub_start > datum$sub_end,]$sub_Orientation <- "Minus" 
  datum$q_Orientation <- "Plus"
  datum[datum$q_start > datum$q_end,]$q_Orientation <- "Minus"
  aligns <- tibble("Position" = c(q_start:q_end),
                   "Aligned_depth_plus" = 0,
                   "Plus_list" = "",
                   "Aligned_depth_minus" = 0,
                   "Minus_list" = "",
                   "Clean_plus" = "No",
                   "Clean_minus" = "No")
  for (i in c(q_start:q_end)){
    zim <- nrow(datum[datum$q_start <= i &
                        datum$q_end >= i &
                        datum$sub_Orientation == "Plus",])
    zim_list <- unique(datum[datum$q_start <= i &
                             datum$q_end >= i &
                             datum$sub_Orientation == "Plus",]$Seq_ID)
    miz <- nrow(datum[datum$q_start <= i &
                        datum$q_end >= i &
                        datum$sub_Orientation == "Minus",])
    miz_list <- unique(datum[datum$q_start <= i &
                             datum$q_end >= i &
                             datum$sub_Orientation == "Minus",]$Seq_ID)
    z_list <- ""
    m_list <- ""
    if (length(zim_list) > 0){
      for (j in zim_list){
        z_list <- paste0(z_list, " ",j)
      }
    }
    if (length(miz_list) > 0){
      for (j in miz_list){
        m_list <- paste0(m_list, " ", j)
      }
    }
    
    aligns[aligns$Position == i,]$Aligned_depth_plus <- zim
    aligns[aligns$Position == i,]$Aligned_depth_minus <- miz
    aligns[aligns$Position == i,]$Plus_list <- z_list
    aligns[aligns$Position == i,]$Minus_list <- m_list
    #print(i)
    #print(zim)
  }
  
  for (i in c(q_start:q_end)){
    #print(i)
    zim_list <- strsplit(aligns[aligns$Position == i,]$Plus_list, " ")[[1]]
    #print("1")
    miz_list <- strsplit(aligns[aligns$Position == i,]$Minus_list, " ")[[1]]
    #print("2")
    for (j in c(i:i+primer_align_threshold)){
      if (j <= q_end){
        zim_list_plus <- strsplit(aligns[aligns$Position == j,]$Plus_list, " ")[[1]]
        zim_list <- intersect(zim_list, zim_list_plus)
        miz_list_plus <- strsplit(aligns[aligns$Position == j,]$Minus_list, " ")[[1]]
        miz_list <- intersect(miz_list, miz_list_plus)
      }
    }
    #print("3")
    if (length(zim_list) <= allowed_hits){
      #print("4")
      #print(i)
      aligns[aligns$Position == i,]$Clean_plus <- "Yes"
      #print(i)
      if (i+primer_align_threshold <= q_end &
          i-(primer_length - primer_align_threshold) >= q_start){
        #print("5")
        interim <- tibble("Primer_start" = i-(primer_length - primer_align_threshold),
                          "Primer_end" = i+primer_align_threshold,
                          "Sense" = "Forward")
        if (exists("the_genome")){
          interim$Sequence <- substr(the_genome, i-(primer_length - primer_align_threshold), i+primer_align_threshold)
        }
        if (!exists("primer_log")){
          primer_log <- interim
        } else {
          primer_log <- rbind(primer_log, interim)
        }
      }
    }
    #print("done plus")
    if (length(miz_list) <= allowed_hits){
      aligns[aligns$Position == i,]$Clean_minus <- "Yes"
      if (i+(primer_length - primer_align_threshold) <= q_end &
          i-primer_align_threshold >= q_start){
        interim <- tibble("Primer_start" = i+(primer_length - primer_align_threshold),
                          "Primer_end" = i-primer_align_threshold,
                          "Sense" = "Reverse")
        if (exists("the_genome")){
          r_check <- substr(the_genome, 
                            i-primer_align_threshold, 
                            i+(primer_length - primer_align_threshold))
          r_check <- str_replace_all(r_check, "G", "c")
          r_check <- str_replace_all(r_check, "C", "g")
          r_check <- str_replace_all(r_check, "A", "t")
          r_check <- str_replace_all(r_check, "T", "a")
          r_check <- toupper(r_check)
          interim$Sequence <- stri_reverse(r_check)
        }
        if (!exists("primer_log")){
          primer_log <- interim
        } else {
          primer_log <- rbind(primer_log, interim)
        }
      }
    }
    #print("Done Minus")
  }
  write_csv(aligns, "align_report.csv")
  if (!exists("primer_log")){
    print("No primers found that avoid alignments.")
  } else if(match_primers) {
    print("Primers identified...")
    
    # Primer pairing and amplicon check
    primer_log$ID <- c(1:nrow(primer_log))
    #print(unique(primer_log$Sense))
    pluses <- primer_log[primer_log$Sense == "Forward",]
    #print(nrow(pluses))
    minuses <- primer_log[primer_log$Sense == "Reverse",]
    #print(nrow(minuses))
    primer_log$Matched <- ""
    if (nrow(pluses) > 0 & nrow(minuses) > 0){
      amp_count <- 0
      for (i in c(1:nrow(pluses))){
        if (amp_count <= number_of_amplicons){
          amp_count <- amp_count + 1
          p_start <- pluses[i,]$Primer_start
          p_id <- unique(pluses[i,]$ID)
          for (j in c(1:nrow(minuses))){
            m_id <- unique(minuses[j,]$ID)
            if (minuses[j,]$Primer_start - p_start >= min_amplicon_length &
                minuses[j,]$Primer_start - p_start <= max_amplicon_length &
                minuses[j,]$Primer_start - p_start > 0){
              primer_log[primer_log$ID == p_id,]$Matched <- paste0(primer_log[primer_log$ID == p_id,]$Matched,
                                                                   m_id, " ")
              primer_log[primer_log$ID == m_id,]$Matched <- paste0(primer_log[primer_log$ID == m_id,]$Matched,
                                                                   p_id, " ")
              interim <- tibble("Foward" = p_id,
                                "F_start" = p_start,
                                "F_seq" =  pluses[i,]$Sequence,
                                "F_length" = nchar(pluses[i,]$Sequence),
                                "F_GC" = round(100*(str_count(pluses[i,]$Sequence, "G")+str_count(pluses[i,]$Sequence, "C"))/nchar(pluses[i,]$Sequence),2),
                                "Reverse" = m_id,
                                "R_start" = minuses[j,]$Primer_start,
                                "R_seq" = minuses[j,]$Sequence,
                                "R_length" = nchar(minuses[j,]$Sequence),
                                "R_GC" = round(100*(str_count(minuses[j,]$Sequence, "G")+str_count(minuses[j,]$Sequence, "C"))/nchar(minuses[j,]$Sequence),2),
                                "Amplicon" = minuses[j,]$Primer_start - p_start,
                                "GC_diff" = abs(round(100*(str_count(pluses[i,]$Sequence, "G")+str_count(pluses[i,]$Sequence, "C"))/nchar(pluses[i,]$Sequence),2)-round(100*(str_count(minuses[j,]$Sequence, "G")+str_count(minuses[j,]$Sequence, "C"))/nchar(minuses[j,]$Sequence),2)))
              if (!exists("the_runs")){
                the_runs <- interim
              } else {
                the_runs <- rbind(the_runs, interim)
              }
            }
          }
        }
      }
    }
  } else {
    print("Skipping primer match...")
  }
  if(match_primers){
    print("Amplicons mapped...")
  }
  
  draft <- ggplot(data = aligns, aes(x = Position, y = Aligned_depth_plus))+
    geom_bar(stat = "identity", aes(fill = Clean_plus, color = Clean_plus))+
    scale_color_manual(values = c("black", "red"))+
    scale_fill_manual(values = c("black", "red"))+
    theme_bw()
  if(exists("primer_log") & map_primers){
    turn <- 0
    if (nrow(primer_log[primer_log$Sense == "Forward",]) > 0){
      chip <- primer_log[primer_log$Sense == "Forward",]
      for (i in c(1:nrow(chip))){
        if (chip[i,]$Matched != ""){
          draft <- draft + geom_segment(x = chip[i,]$Primer_start, 
                                        xend = chip[i,]$Primer_end,
                                        y = max(aligns$Aligned_depth_plus)+1+turn,
                                        yend = max(aligns$Aligned_depth_plus)+1+turn, color = "green")
        } else if(map_primers) {
          draft <- draft + geom_segment(x = chip[i,]$Primer_start, 
                                        xend = chip[i,]$Primer_end,
                                        y = max(aligns$Aligned_depth_plus)+1+turn,
                                        yend = max(aligns$Aligned_depth_plus)+1+turn, color = "red")
        }
        turn <- turn+1
      }
    }
    draft <- draft+ylim(0, max(aligns$Aligned_depth_plus)+2+nrow(primer_log[primer_log$Sense == "Forward",]))
  } else {
    draft <- draft+ylim(0, max(aligns$Aligned_depth_plus)+10)
  }
  print(draft)
  ggsave("forward_primers.png", dpi = 150, width = 3000, height = 750, units = "px")
  draft <- ggplot(data = aligns, aes(x = Position, y = Aligned_depth_minus))+
    geom_bar(stat = "identity", aes(fill = Clean_minus, color = Clean_minus))+
    scale_color_manual(values = c("black", "blue"))+
    scale_fill_manual(values = c("black", "blue"))+
    theme_bw()
  if(exists("primer_log") & map_primers){
    turn <- 0
    if (nrow(primer_log[primer_log$Sense == "Reverse",]) > 0){
      chip <- primer_log[primer_log$Sense == "Reverse",]
      for (i in c(1:nrow(chip))){
        if (chip[i,]$Matched != ""){
          draft <- draft + geom_segment(x = chip[i,]$Primer_start, 
                                        xend = chip[i,]$Primer_end,
                                        y = max(aligns$Aligned_depth_minus)+1+turn,
                                        yend = max(aligns$Aligned_depth_minus)+1+turn, color = "green")
        } else if(map_primers) {
          draft <- draft + geom_segment(x = chip[i,]$Primer_start, 
                                        xend = chip[i,]$Primer_end,
                                        y = max(aligns$Aligned_depth_minus)+1+turn,
                                        yend = max(aligns$Aligned_depth_minus)+1+turn, color = "blue")
        }
        
        turn <- turn+1
      }
    }
    draft <- draft+ylim(0, max(aligns$Aligned_depth_minus)+2+nrow(primer_log[primer_log$Sense == "Reverse",]))
  } else {
    draft <- draft+ylim(0, max(aligns$Aligned_depth_minus)+5)
  }
  print(draft)
  ggsave("reverse_primers.png", dpi = 150, width = 3000, height = 750, units = "px")
  
  if (exists("primer_log")){
    #print(primer_log)
    write_csv(primer_log, "primer_report.csv")
    if (exists("the_runs")){
      write_csv(the_runs, "amplicon_report.csv")
    } else {
      print("No primer pairs detected for current parameters.")
    }
  }
  print("Analysis complete.")
}