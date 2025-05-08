suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))

brainbox <- function(input_file,
                     data_file,
                     viridis_palette = "inferno",
                     threshold=0,
                     output_file="hexed_data.csv"){
  datum <- read_csv(data_file)
  if (ncol(datum) > 2){
    print(names(datum))
    num_column <- as.integer(readline(prompt = paste0("Data file contains ", ncol(datum), " columns. Which column number represents the base number: ")))
    dat_column <- as.integer(readline(prompt = paste0("Data file contains ", ncol(datum), " columns. Which column number represents the data to color by: ")))
    datum <- datum[,c(num_column, dat_column)]
  }
  gg <- ggplot(data = datum, aes(x = unlist(datum[,2]), y = 1, color = unlist(datum[,2]), fill = unlist(datum[,2])))+
    geom_bar(stat = "identity")+
    scale_color_viridis_c(option = viridis_palette)+
    scale_fill_viridis_c(option = viridis_palette)+
    theme_bw()
  print(gg)
  ggsave("legend.png", dpi = 300, width = 6, height = 6)
  datum$hex <- ggplot_build(gg)[[1]][[1]]$colour
  datum[unlist(datum[,2]) < threshold,]$hex <- "#000000"
  write_csv(datum, output_file)
}

brainbox(input_file = "CY1 Structure.rnacanvas", data_file = "hexed_data_t1.csv", viridis_palette = "plasma", threshold = 0, output_file = "t1.csv")
