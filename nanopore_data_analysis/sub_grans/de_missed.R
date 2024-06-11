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
