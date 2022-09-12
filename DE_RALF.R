library(edgeR)  # загрузка библиотеки

# загружаем данные которые сохранили в прошлый раз
counts <-  readRDS('counts.Rdata')
head(counts)

################# DE ANALYSIS ##############

#выбираем данные для ральфов
fc_ralf = counts[, c("VB_4_S4","VB_6_S6", "VB_8_S8","VB_10_S10", "VB_12_S12")]
head(fc_ralf)

saveRDS(fc_ralf,'fc_ralf.Rdata')

# genes without expression
table(apply(fc_ralf,1,mean) == 0)

# genes having coverage above 10 reads
table(apply(fc_ralf,1,mean) >= 5)

#задаем метаданные для образцов
meta <-  data.frame(infection=c("fusarium","mock", "fusarium", "mock", "fusarium"))
rownames(meta) = colnames(fc_ralf)
table(meta)
meta

## Загружаем данные в edgeR и нормализуем
edger <-  DGEList(counts = fc_ralf, group = c(1,0,1,0,1)) # создаем объект edgeR хранящий каунты

#фильтруем по покрытию
keep <- filterByExpr(edger)
table(keep)
edger <- edger[keep, , keep.lib.sizes=FALSE]

# нормируем методом RLE
edger = calcNormFactors(edger,method='TMM')
edger$samples

cpm = cpm(edger) # посчитаем cpm с учетом TMM нормировки
write.csv(cpm, file = "cpm_values_ralf.csv", row.names = TRUE) #выводим таблицу

#создаем матрицу дизайна
#мы будем искать отличия между зараженными растениями и контролем (infection),
design = model.matrix(~ infection, data = meta)
design

# biological coefficient of variation (BCV)
# рисуем зависимость биологической вариабельности от средней экспрессии, 
#немного падает с ростом экпсрессии
edger <-  estimateDisp(edger,design)
plotBCV(edger)
edger$common.dispersion

plotMDS(edger)

# fit GLM model
glm <-  glmFit(edger,design)
#plotQLDisp(glm)

# calculating DE genes during infection
glf <- glmLRT(glm)

topTags(glf)#наиболее меняющиеся гены

#summary statistics
summary(decideTests(glf))

#Plot log-fold change against log-counts per million, with DE genes highlighted:
plotMD(glf)
abline(h=c(-1, 1), col="blue")

# calculate FDR
glf$table$FDR <- p.adjust(glf$table$PValue, method="BH")
#распределение p-value
hist(glf$table$PValue)

#CPM top DE genes
top <- rownames(topTags(glf, n=20))
cpm[top,]

#Выводим данные
write.csv(glf$table, file = "de_test_ralf.csv", row.names = TRUE) #выводим таблицу

#Всего меняющихся генов
sum(glf$table$FDR < 0.05)

