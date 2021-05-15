setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/IPA/ImportFile")
library(readxl)
em <- read_excel("Proteomics_LFQ&SILAC_EM.xlsx", sheet = 1)
em <- data.frame(em) # 883

he <- read_excel("Proteomics_LFQ&SILAC_HE.xlsx", sheet = 1)
he <- data.frame(he) # 672

hm <- read_excel("Proteomics_LFQ&SILAC_HM.xlsx", sheet = 1)
hm <- data.frame(hm) # 565

table(duplicated(em$Gene.Name)) # 22 duplicates
table(duplicated(he$Gene.Name)) # 15 duplicates
table(duplicated(hm$Gene.Name)) # 13 duplicates

em <- em[-which(duplicated(em$Gene.Name)), ] # 861
he <- he[-which(duplicated(he$Gene.Name..GN.)), ] # 657
hm <- hm[-which(duplicated(hm$Gene.Name..GN.)), ] # 552

#----------------------------------------------------------------------------#
# Import SILAC data
setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018")
silac <- read_excel("Proteomics_SILAC_summary_score_28.12.2018.xlsx", sheet = 1)
silac <- data.frame(silac)
silac <- silac[-which(duplicated(silac$Gene.Name)), ] # 7082
silac <- silac[-which(is.na(silac$Gene.Name)), ] # 7081

#----------------------------------------------------------------------------#
for (i in 1:nrow(silac)){
  if (silac$Gene.Name[i] %in% em$Gene.Name){
    silac$"EM"[i] <- "YES"
  } else {
    silac$"EM"[i] <- "NO"
  }
  
  if (silac$Gene.Name[i] %in% he$Gene.Name){
    silac$"HE"[i] <- "YES"
  } else {
    silac$"HE"[i] <- "NO"
  }
  
  if (silac$Gene.Name[i] %in% hm$Gene.Name){
    silac$"HM"[i] <- "YES"
  } else {
    silac$"HM"[i] <- "NO"
  }
}

# check
table(silac$EM) # 859/861
table(silac$HE) # 656/657
table(silac$HM) # 551/552

# Find out which are not 
silac_em <- silac[which(silac$EM == "YES"), ] # 859
silac_he <- silac[which(silac$HE == "YES"), ] # 656
silac_hm <- silac[which(silac$HM == "YES"), ] # 551

`%!in%` <- Negate(`%in%`)
em$Gene.Name[em$Gene.Name %!in% silac_em$Gene.Name] # "GLUL" (A8YXX4), "ACSL1" (A8K9T3)
he$Gene.Name[he$Gene.Name %!in% silac_he$Gene.Name] # "GLUL" (A8YXX4)
hm$Gene.Name[hm$Gene.Name %!in% silac_hm$Gene.Name] # "ACSL1" (A8K9T3)

#-----------------------------------------------------------------#
# export the data from Perseus
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/Perseus/EM&HE&HM_06.02.2019_LFQ&SILAC_significant")
write.table(silac, "PerseusInPutTable_SILAC_EM&HE&HM_07.02.2019.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

#----------------------------------------------------------------------------#
# add one categorical column
# D492M vs. D492
for (i in 1:nrow(em)){
  if (em$Mean.D492M.D492._LFQ_SILAC[i] >= 1 | em$Mean.D492M.D492._LFQ_SILAC[i] <= -1){
    em$"Grouping"[i] <- "Group 1"
  } else if ((em$Mean.D492M.D492._LFQ_SILAC[i] >= 0.5 & em$Mean.D492M.D492._LFQ_SILAC[i] < 1) | (em$Mean.D492M.D492._LFQ_SILAC[i] <= -0.5 & em$Mean.D492M.D492._LFQ_SILAC[i] > -1)){
    em$"Grouping"[i] <- "Group 2"
  } else {
    em$"Grouping"[i] <- "Group 3"
  }
}

for (i in 1:nrow(he)){
  if (he$Mean.D492HER2.D492._LFQ_SILAC[i] >= 1 | he$Mean.D492HER2.D492._LFQ_SILAC[i] <= -1){
    he$"Grouping"[i] <- "Group 1"
  } else if ((he$Mean.D492HER2.D492._LFQ_SILAC[i] >= 0.5 & he$Mean.D492HER2.D492._LFQ_SILAC[i] < 1) | (he$Mean.D492HER2.D492._LFQ_SILAC[i] <= -0.5 & he$Mean.D492HER2.D492._LFQ_SILAC[i] > -1)){
    he$"Grouping"[i] <- "Group 2"
  } else {
    he$"Grouping"[i] <- "Group 3"
  }
}

for (i in 1:nrow(hm)){
  if (hm$Mean.D492HER2.D492M._LFQ_SILAC[i] >= 1 | hm$Mean.D492HER2.D492M._LFQ_SILAC[i] <= -1){
    hm$"Grouping"[i] <- "Group 1"
  } else if ((hm$Mean.D492HER2.D492M._LFQ_SILAC[i] >= 0.5 & hm$Mean.D492HER2.D492M._LFQ_SILAC[i] < 1) | (hm$Mean.D492HER2.D492M._LFQ_SILAC[i] <= -0.5 & hm$Mean.D492HER2.D492M._LFQ_SILAC[i] > -1)){
    hm$"Grouping"[i] <- "Group 2"
  } else {
    hm$"Grouping"[i] <- "Group 3"
  }
}

#-------------------------------------------------------------------#
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/Perseus/EM&HE&HM_06.02.2019_LFQ&SILAC_significant")
write.table(em, "PerseusInPutTable_EM_07.02.2019.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(he, "PerseusInPutTable_HE_07.02.2019.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(hm, "PerseusInPutTable_HM_07.02.2019.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
