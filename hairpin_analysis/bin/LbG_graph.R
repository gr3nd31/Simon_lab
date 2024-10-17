#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (system.file(package = 'ggpubr') == ""){
  install.packages('ggpubr')
}

if (system.file(package = 'tidyverse') == ""){
  install.packages('tidyverse')
}

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))

if (length(args) == 0){
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (!grepl("csv", args[1])){
  stop("Argument file should be a CSV file.n", call.=FALSE)
}

datum <- read_csv(args[1], show_col_types = F)
the_name <- str_replace(args[1], ".csv", "_")

if (nrow(datum) > nrow(datum[str_count(datum$Structure, "\\(") > 3,])){
  cat("Detected sequences with less than 3 basepairs. These are being trimmed")
  datum <- datum[str_count(datum$Structure, "\\(")>= 3,]
}

m = lm(dG~Length, data = datum)[[1]][2]
b = lm(dG~Length, data = datum)[[1]][1]

draft <- ggplot(data = datum, aes(x = Length, y = dG))+
  geom_point(size = 5, alpha = 0.7)+
  theme_bw()+
  geom_smooth(method = "lm", linewidth = 2, alpha = 0.5, linetype = 2)+
  scale_x_continuous(position = "top")+
  labs(x = "Length (nt)",
       y = "Minimum Free Energy",
       title = paste0("y = ", round(m, 2), "x + ", round(b,2)))+
  theme(line = element_line(linewidth = 1),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.line = element_line(linewidth = 1),
        title = element_text(size = 20))
ggsave(paste0(the_name,"Length_by_dG.png"), dpi = 300, plot = draft, width = 6, height = 6)

m = lm(dG~APE, data = datum)[[1]][2]
b = lm(dG~APE, data = datum)[[1]][1]

draft <- ggplot(data = datum, aes(x = Length, y = APE))+
  geom_point(size = 5, alpha = 0.7)+
  theme_bw()+
  geom_smooth(method = "lm", linewidth = 2, alpha = 0.5, linetype = 2)+
  labs(x = "Length (nt)",
       y = "Average Positional Entropy",
       title = paste0("y = ", round(m, 2), "x + ", round(b,2)))+
  theme(line = element_line(linewidth = 1),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.line = element_line(linewidth = 1),
        title = element_text(size = 20))
ggsave(paste0(the_name,"Length_by_APE.png"), dpi = 300, plot = draft, width = 6, height = 6)

m = lm(dG~log(APE), data = datum)[[1]][2]
b = lm(dG~log(APE), data = datum)[[1]][1]
draft <- ggplot(data = datum, aes(x = dG_Length, y = log(APE)))+
  geom_point(size = 5, alpha = 0.7)+
  theme_bw()+
  geom_smooth(method = "lm", linewidth = 2, alpha = 0.5, linetype = 2)+
  labs(x = "MFE/Length",
       y = "Average Positional Entropy (log10)",
       title = paste0("y = ", round(m, 2), "x + ", round(b,2)))+
  theme(line = element_line(linewidth = 1),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.line = element_line(linewidth = 1),
        title = element_text(size = 20))
ggsave(paste0(the_name,"dGLength_by_APE.png"), dpi = 300, plot = draft, width = 6, height = 6)
