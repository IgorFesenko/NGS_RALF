if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

rm(list=ls())

library(edgeR)  # загрузка библиотеки

# Загружаем данные по каунтам после featureCounts
fc <- read.table('RALF_ngs', sep='\t', header=TRUE)
colnames(fc)
dim(fc)

# имена колонок содержать ".bam", уберем:
colnames(fc) <-  sub('_ME.bam','',colnames(fc))
colnames(fc) <-  sub('..bam.','',colnames(fc))
colnames(fc)

# оставляем нужные столбцы
filt_fc <-  fc[,c("Geneid", "VB_1_S1","VB_2_S2","VB_3_S3","VB_4_S4","VB_5_S5",
                  "VB_6_S6","VB_7_S7","VB_8_S8","VB_9_S9","VB_10_S10","VB_11_S11","VB_12_S12")]
row.names(filt_fc) <- filt_fc$Geneid
filt_fc$Geneid <- NULL

# Сохраняем промежуточные данные
saveRDS(filt_fc,'counts.Rdata') # save counts
write.csv(filt_fc, file = "count_values.csv", row.names = TRUE) #выводим таблицу

## Предварительный анализ данных для MDS PCA
# Нормализуем данные
edger = DGEList(filt_fc) # создаем объект edgeR хранящий каунты
# нормируем методом RLE
edger = calcNormFactors(edger,method='RLE')
edger$samples
edger$dispersion # смотрим на дисперсию
# get RLE-norm cpms
cpm = cpm(edger) # посчитаем cpm с учетом RLE нормировки
write.csv(cpm, file = "cpm_values.csv", row.names = TRUE) #выводим таблицу

#убираем второй образец и подсчитываем все данные заново
filt_fc_no2 <-  filt_fc[,c("VB_1_S1","VB_3_S3","VB_4_S4","VB_5_S5", "VB_6_S6",
                           "VB_7_S7","VB_8_S8","VB_9_S9","VB_10_S10",
                           "VB_11_S11","VB_12_S12")]

# filtration based on coverage
# distribution according to reads coverage
hist(log10(1+apply(filt_fc_no2,1,mean)),20)
abline(v=log10(1+10),col='red') # threshold in 10 reads

# genes without expression
table(apply(filt_fc_no2,1,mean) == 0)

# genes having coverage above 10 reads
table(apply(filt_fc_no2,1,mean) >= 10) 

# filtering genes  - above 10 reads per sample
fc10 = filt_fc_no2[apply(filt_fc_no2,1,mean) >= 10,]
table(apply(fc10,1,mean) >= 10)

# distribution according to reads coverage
hist(log10(1+apply(fc10,1,mean)),20)

# save data
saveRDS(fc10,'filtered_counts.Rdata')

# создадим табличку с информацией об образцах. 
meta <-  data.frame(infection=c("mock","fusarium","fusarium",
                                "mock","mock","fusarium","fusarium",
                                "mock","mock","fusarium","fusarium"),
                    genotype=c('WT','WT','RALF3',
                               'WT','RALF3','WT','RALF3',
                               'WT','RALF3','WT','RALF3'),
                   color=c("green","orange","red",
                            "green","blue","orange","red",
                            "green","blue","orange","red"))

rownames(meta) = colnames(filt_fc_no2)
head(meta)
# save data to file
saveRDS(meta,'metadata.Rdata')
write.csv(meta, file = "metadata.csv", row.names = TRUE) #выводим таблицу

# Нормализуем данные
edger = DGEList(fc10) # создаем объект edgeR хранящий каунты
# нормируем методом RLE
edger = calcNormFactors(edger,method='RLE')
edger$samples
edger$dispersion # смотрим на дисперсию

cpm = cpm(edger) # посчитаем cpm с учетом RLE нормировки
write.csv(cpm, file = "cpm_values_above10.csv", row.names = TRUE) #выводим таблицу


#СОЗДАЕМ МАТРИЦУ ДИЗАЙНА
#мы будем искать отличия между зараженными растениями и контролем (infection) 
#для двух генотипов

design = model.matrix(~infection + genotype + infection:genotype,data = meta)
design

# biological coefficient of variation (BCV)
edger <-  estimateDisp(edger,design)
plotBCV(edger)
edger$common.dispersion


