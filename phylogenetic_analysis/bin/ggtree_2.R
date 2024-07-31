suppressMessages(library(ggtree))
suppressMessages(library(ggplot2))
suppressMessages(library(treeio))
suppressMessages(library(tidyverse))

the_db <- read_csv("nucleotide_list.csv")
tree <- read.tree("RAxML_bipartitions.sequences_tree")
the_db <- the_db[the_db$Name %in% tree$tip.label,]
colorz <- c("lightgreen", "darkgreen", "purple", "gold", "blue")

draft <- ggtree(tree, size = 1.5, linetype = 1)
clades <- c()
if (length(the_db$Tag > 1)){
  for (i in c(1:length(unique(the_db$Tag)))){
    trick <- unique(the_db$Tag)[i]
    print(trick)
    print(unlist(the_db[the_db$Tag == trick,]$Name))
    clade_call <- MRCA(tree, unlist(the_db[the_db$Tag == trick,]$Name))
    clades <- append(clades, clade_call)
  }
}
tree2 <- groupClade(tree, clades[c(1:5)])
draft <- ggtree(tree, size = 1.5, linetype = 1)+
  geom_tiplab(align = T, hjust = 0.2)+
  geom_treescale(linesize = 1.5)+
  #geom_text(aes(label=node, hjust=-0.3))+
  theme_tree()

print(draft)
ggsave("RdRp_tree.svg")
