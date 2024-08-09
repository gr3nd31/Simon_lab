suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))

alignify <- function(file = "sequences_aligned_matrix.csv",
                     threshold = 70,
                     nameer = "none", 
                     return_tops = F) {
  datum <- read_csv(file)
  counterNum <- 0
  occurs <- length(unique(datum$Name))
  
  for (i in unique(datum$Position)){
    most_common <- names(sort(table(datum[datum$Position==i,]$Base), decreasing = T)[1])[1]
    common_count <- sort(table(datum[datum$Position==i,]$Base), decreasing = T)[1]
    if ("CY1" %in% unique(datum[datum$Position ==i & datum$Base != "-",]$Name)){
      counterNum <- counterNum+1
    }
    interim <- tibble("Position" = i,
                      "Percent" = round(100*common_count/occurs, 2),
                      "Consensus" = most_common,
                      "Counter" = counterNum)
    if (!exists("the_data")){
      the_data <- interim
    } else {
      the_data <- rbind(the_data, interim)
    }
  }
  the_data$Threshold <- "No"
  the_data[the_data$Percent > threshold & the_data$Consensus != "-",]$Threshold <- "Yes"
  write_csv(the_data, file)
  if (return_tops & nameer != "none"){
    the_list <- ""
    the_data <- the_data[the_data$Threshold == "Yes",]
    for (i in unique(the_data$Counter)){
      the_list <- paste0(the_list, i,"+")
    }
    cat(the_list)
  }
}
alignify()

# threshold <- 99
# datum <- read_csv("sequences_aligned_matrix.csv")
# the_list <- ""
# datum$Threshold <- "No"
# datum[datum$Percent > threshold & datum$Consensus != "-",]$Threshold <- "Yes"
# datum <- datum[datum$Threshold == "Yes",]
# for (i in unique(datum$Counter)){
#   the_list <- paste0(the_list, i,"+")
# }
# the_list <- substr(the_list, 1, nchar(the_list)-1)
# cat(the_list)
