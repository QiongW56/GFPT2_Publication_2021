setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/Reactome/Metabolic Proteins_H vs M")
library(readxl)
dat <- read_excel('ReactomeOutputTable.xlsx', sheet = 2)
dat <- data.frame(dat)
dat <- dat[-1:-2, ] # delete "Metabolism of proteins" and "Metabolism of RNA"
dat <- dat[, c(2, 3, 7)] # "Pathway.name", "X.Entities.found", "Entities.FDR"
colnames(dat)[3] <- c('FDR')

# install.packages('treemap')

library(treemap)
library(grDevices)

setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018/QIONG_Functions")
load('treemap.plot.Rdata')

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/Reactome/Metabolic Proteins_H vs M")
treemap.plot(dat, 
             title = 'Metabolic Dysregulation in D492HER2 vs. D492M', 
             color =  "Blues",
             label.size = 12,
             lgd = 16,
             width = 8)
