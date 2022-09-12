library(GOplot)
library(readr)

rm(list=ls())


# загружаем данные
david <- read_csv('david_de_wt.csv')
head(david)

genelist <- read_csv('genelist_de_wt.csv')
head(genelist)

## Generate the plotting object
circ <- circle_dat(david, genelist)
head(circ)

#Визуализируем
GOBubble(circ, labels = 2.7, table.legend = T)

GOBar(subset(circ, category == 'BP'))



#Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.9)
GOBubble(reduced_circ, labels = 2.5)

GOBar(subset(reduced_circ, category == 'BP'))


