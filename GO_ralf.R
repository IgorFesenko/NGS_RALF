library(GOplot)
library(readr)

rm(list=ls())


# загружаем данные
david <- read_csv('david_de_ralf.csv')
head(david)

genelist <- read_csv('genelist_de_ralf.csv')
head(genelist)

## Generate the plotting object
circ <- circle_dat(david, genelist)
head(circ)

#Визуализируем
options(repr.plot.width = 4, repr.plot.height = 7)

GOBubble(circ, labels = 4, table.legend = T,bg.col = T,)


#Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.8)
GOBubble(reduced_circ, labels = 2.5, table.legend = T,)
