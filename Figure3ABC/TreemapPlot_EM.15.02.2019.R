setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/Reactome/Metabolic Proteins_E vs M")
library(readxl)
dat <- read_excel('ReactomeResult_FDR0.05_Proteins different in EMT_LFQ&SILAC.xlsx', sheet = 2)
dat <- data.frame(dat)
dat <- dat[-1:-2, ] # delete "Metabolism of RNA" and "Metabolism of proteins"
dat <- dat[, c(2, 3, 7)]
colnames(dat)[3] <- c('FDR')

# install.packages('treemap')

library(treemap)
library(grDevices)

setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018/QIONG_Functions")
load('treemap.plot.Rdata')

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/Reactome/Metabolic Proteins_E vs M")
treemap.plot(dat, 
             title = 'Metabolic Dysregulation in D492M vs. D492', 
             color =  "YlOrBr",
             label.size = 12,
             lgd = 16,
             width = 8)
