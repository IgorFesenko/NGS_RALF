library(GOplot)
library(readr)

rm(list=ls())


# загружаем данные
david <- read_csv('david_de_wt_common.csv')
head(david)

genelist <- read_csv('genelist_de_wt_common.csv')
head(genelist)

## Generate the plotting object
circ <- circle_dat(david, genelist)
head(circ,10)

#Визуализируем
GOBubble(circ, labels = 3, table.legend = T, colour = c('orange', 'gold','darkred'))


#Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.8)
GOBubble(reduced_circ, labels = 3)
