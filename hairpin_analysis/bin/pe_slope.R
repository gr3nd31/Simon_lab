suppressPackageStartupMessages(library(tidyverse))

pe_slope <- function(df = datum,
                     base_count = 30,
                     make_graphs = F){
  df$bp_percent <- 2*df$bp/df$Length
  df$pe_slope <- 0
  df$basal_pe_slope <- 0
  df$basal_pe_intercept <- 0
  df$apical_pe_slope <- 0
  df$apical_pe_intercept <- 0
  df$pe_intercept <- 0
  df$compressed_percent <- 0
  df$GC_compressed_percent <- 0
  df$AU_compressed_percent <- 0
  for (i in c(1:nrow(df))){
    j <- df[i,]$PE
    j <- str_remove(string = j, pattern = "\\(")
    j <- str_remove(string = j, pattern = "\\)")
    j <- str_split(j, " ")[[1]]
    j <- j[2:length(j)]
    j <- as.numeric(j)
    interim <- tibble("Position" = c(1:length(j)),
                      "PE" = j)
    midpoint <- as.integer(length(j)/2)
    interim$relative_Position <- interim$Position
    interim[interim$Position > midpoint,]$relative_Position <- 1+abs(interim[interim$Position > midpoint,]$Position-max(interim$Position))
    pe_slope <- lm(PE~relative_Position, data = interim)
    df[i,]$pe_slope <- pe_slope$coefficients[2]
    df[i,]$pe_intercept <- pe_slope$coefficients[1]
    
    pe_slope <- lm(PE~relative_Position, data = interim[interim$relative_Position <= base_count,])
    df[i,]$basal_pe_slope <- pe_slope$coefficients[2]
    df[i,]$basal_pe_intercept <- pe_slope$coefficients[1]
    
    pe_slope <- lm(PE~relative_Position, data = interim[interim$relative_Position >= (max(interim$relative_Position)-base_count),])
    df[i,]$apical_pe_slope <- pe_slope$coefficients[2]
    df[i,]$apical_pe_intercept <- pe_slope$coefficients[1]
    
    if (make_graphs){
      draft <- ggplot(data = interim, aes(x = relative_Position, y = PE))+
        geom_point()+
        geom_smooth(method = "lm")+
        theme_bw()
      print(draft)
      ggsave(paste0(str_replace_all(df[i,]$Name, ">", ""),".png"), dpi = 300)
    }
    
    the_string <- df[i,]$Sequence
    the_string <- str_replace_all(the_string, "C+", "C")
    the_string <- str_replace_all(the_string, "G+", "G")
    df[i,]$GC_compressed_percent <- nchar(the_string)/df[i,]$Length
    the_string <- df[i,]$Sequence
    the_string <- str_replace_all(the_string, "A+", "A")
    the_string <- str_replace_all(the_string, "U+", "U")
    df[i,]$AU_compressed_percent <- nchar(the_string)/df[i,]$Length
    the_string <- str_replace_all(the_string, "C+", "C")
    the_string <- str_replace_all(the_string, "G+", "G")
    df[i,]$compressed_percent <- nchar(the_string)/df[i,]$Length
  }
  return(df)
}

find_a_hairpin <- function(df = datum,
                           shifter = 0.1,
                           target_ape = 0.2,
                           target_dGslope = -0.41, 
                           target_PEslope = 0.1,
                           target_GC = 0.5,
                           bulge = "None",
                           bulge_size = 0.1,
                           bulge_position = 0.5,
                           complementarity = 100,
                           target_structure = "",
                           cutoff = 3, 
                           out_file = "targets.csv"){
  print("Starting search")
  maxer <- 1
  miner <- 1
  if (!"pe_slope" %in% names(df)){
    cat("Running PE slope calculations...\n")
    df <- pe_slope(df)
  }
  if (target_structure != ""){
    df <- df[df$Structure == target_structure,]
  }
  
  if (nrow(df) < 1){
    cat("No sequences with given structure. Aborting search.\n")
  } else {
    cat("Running initial search...\n")
    interim <- df[df$APE == target_ape &
                    df$pe_slope == target_PEslope &
                    df$GC == target_GC & 
                    df$Bulge == bulge & 
                    df$BulgeSize == bulge_size &
                    df$Complementarity == complementarity,]
    if (nrow(interim) < 1){
      cat("Initial run did not find a sequence. Relaxing search parameters...")
      enough <- FALSE
    } else {
      if_enough <- readline(prompt = paste0("Found ", nrow(interim), " sequences with desired targets. Is this enough (Y/n)?"))
      if (if_enough != "n"){
        enough <- TRUE
      } else{
        enough <- FALSE
      }
    }
    # Iterates through target parameter relaxation until hits are found
    while (nrow(interim) < 1 | enough == FALSE){
      maxer <- maxer+shifter
      miner <- miner-shifter
      interim <- df[df$APE < target_ape*maxer & df$APE > target_ape*miner &
                     df$pe_slope < target_PEslope*maxer & df$pe_slope > target_PEslope*miner &
                     df$GC < target_GC*maxer & df$GC > target_GC*miner &
                     df$Bulge == bulge & 
                     df$BulgeSize < bulge_size*maxer & df$BulgeSize > bulge_size*miner &
                     df$Complementarity < complementarity*maxer & df$Complementarity > complementarity*miner,]
      # Confirms there are enough target sequences
      if (nrow(interim) > 0){
        if_enough <- readline(prompt = paste0("Found ", nrow(interim), " sequences with parameters within ", (maxer-1)*100,"% of targets. Is this enough (Y/n)?"))
        if (if_enough != "n"){
          enough <- TRUE
        }
      } else if (maxer >= cutoff){
        cat("Relaxation of parameters has surpassed cutoff threshold. Aborting search.\n")
        break
      }
    }
    if (enough | nrow(interim) > 0){
      cat(paste0("Found ", nrow(interim), " sequences within ", (maxer-1)*100, "% of the target parameters.\n"))
      write_csv(interim, out_file)
    } else {
      cat("Unable to find matching sequences. Run more iterations and try again.\n")
      print(maxer)
    }
  }
}
#find_a_hairpin(datum)
