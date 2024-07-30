#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0){
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1){
  args[2] = "rectangular"
}

#print(length(args))
#print(args[2])

suppressMessages(library(ggtree))
suppressMessages(library(ggplot2))
suppressMessages(library(treeio))
tree <- read.tree(args[1])
draft <- ggtree(tree, 
                layout = args[2], 
                size = 1.5, 
                linetype = 1,)+
  geom_treescale(linesize = 1.5)+
  geom_tiplab()
print(draft)
ggsave("sequences_tree.svg", width = 14, height = 7)