#------------------------------------------------------------------------------------#

# the same data used for survival curve for metabolic targets

#------------------------------------------------------------------------------------#
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/SurvivalCurveGeneration_17.09.2019")
# subtype, survival ("sub_vec", "age", "patient", "status", "survival", "cluster", "ER")
dd_subtype <- readRDS(file = "Subtyped_TCGA_data.rds") # 1066 * 7

# all genes for all cases (20247 * 1066)
# NON-z-scored TCGA data (but normalized)
vsd_tcga <- readRDS(file = 'TCGA_data_whole_vst_transformed_170519.rds') #DEseq package


# GFPT2
vsd_tcga_GFPT2 <- vsd_tcga[which(rownames(vsd_tcga) == "GFPT2"), ]

vsd_tcga_GFPT2 <- as.data.frame(t(vsd_tcga_GFPT2), stringsAsFactors = FALSE)

vsd_tcga_GFPT2$patient <- row.names(vsd_tcga_GFPT2)

library(reshape2)
dat <- merge(dd_subtype, vsd_tcga_GFPT2)

# Z-score 
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

dat[, ncol(dat)] <- cal_z_score(dat[, ncol(dat)])

bl <- mean(dat$GFPT2[which(dat$sub_vec == "Basal-like")])
her2 <- mean(dat$GFPT2[which(dat$sub_vec == "HER2-enriched")])
LA <- mean(dat$GFPT2[which(dat$sub_vec == "Luminal A")])
LB <- mean(dat$GFPT2[which(dat$sub_vec == "Luminal B")])
nl <- mean(dat$GFPT2[which(dat$sub_vec == "Normal-like")])

table(dat$sub_vec)

dat_mean <- data.frame(names(table(dat$sub_vec)), c( bl, her2, LA, LB, nl))
colnames(dat_mean) <- c("Tumor_Type", "GFPT2_Avg.Expression")


#------------------------------------------------------------------------------------#

# Metaboric: C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\
# 5 ProteomicPaper\Figures&Tables in the paper\Figures\Figure 10\RawData\brca_metabric

#------------------------------------------------------------------------------------#
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Figures/Figure 10/RawData/brca_metabric")
datMetExp <- read.delim("data_expression_median.txt")

# GFPT2
datMetGF <- datMetExp[which(datMetExp$Hugo_Symbol == "GFPT2"), ]
datMetGF <- as.data.frame(t(datMetGF), stringsAsFactors = FALSE)
datMetGF$patient <- row.names(datMetGF)
datMetGF <- datMetGF[-1:-2, ]

datMetGF[, 1] <- as.numeric(datMetGF[, 1])

# patient tumor types
datMetType <- read.delim("data_clinical_patient.txt")

for (i in 1:ncol(datMetType)){
  datMetType[i] <- as.character.factor(datMetType[, i])
}

datMetType <- datMetType[-1:-4, ]

datMetType[, 1] <- sub("-", ".", datMetType[, 1])

colnames(datMetType)[1] <- "patient"

library(reshape2)
dat <- merge(datMetType, datMetGF)
colnames(dat)[ncol(dat)] <- "GFPT2"

# Z-score 
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

dat[, ncol(dat)] <- cal_z_score(dat[, ncol(dat)])

colnames(dat)[15] <- "Subtype"

table(dat$Subtype)

bl <- median(dat$GFPT2[which(dat$Subtype == "Basal")], na.rm = TRUE)
clauLow <- median(dat$GFPT2[which(dat$Subtype == "claudin-low")], na.rm = TRUE)
her2 <- median(dat$GFPT2[which(dat$Subtype == "Her2")], na.rm = TRUE)
LA <- median(dat$GFPT2[which(dat$Subtype == "LumA")], na.rm = TRUE)
LB <- median(dat$GFPT2[which(dat$Subtype == "LumB")], na.rm = TRUE)
NC <- median(dat$GFPT2[which(dat$Subtype == "NC")], na.rm = TRUE)
nl <- median(dat$GFPT2[which(dat$Subtype == "Normal")], na.rm = TRUE)

dat_median <- data.frame(names(table(dat$Subtype)), c( bl, clauLow, her2, LA, LB, NC, nl))
colnames(dat_median) <- c("Tumor_Type", "GFPT2_Avg.Expression")

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Figures/Figure 10/Output")
library("xlsx")
write.xlsx(dat_median, "Metabric.xlsx", sheetName = "GFPT2_BC_Type_median", row.names = FALSE, append = TRUE)

#------------------------------------------------------------------------------------#

# Pan-dataset: C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\
# 5 ProteomicPaper\Figures&Tables in the paper\Figures\Figure 10\RawData\brca_tcga_pan_can_atlas_2018

#------------------------------------------------------------------------------------#
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Figures/Figure 10/RawData/brca_tcga_pan_can_atlas_2018")
datPanExp <- read.delim("data_RNA_Seq_v2_expression_median.txt")

# GFPT2
datPanGF <- datPanExp[which(datPanExp$Hugo_Symbol == "GFPT2"), ]
datPanGF <- as.data.frame(t(datPanGF), stringsAsFactors = FALSE)
datPanGF$patient <- row.names(datPanGF)
datPanGF <- datPanGF[-1:-2, ]
datPanGF[, 2] <- sub(".01", "", datPanGF[, 2])

datPanGF[, 1] <- as.numeric(datPanGF[, 1])

# patient tumor types
datPanType <- read.delim("data_clinical_patient.txt")

for (i in 1:ncol(datPanType)){
  datPanType[i] <- as.character.factor(datPanType[, i])
}

datPanType <- datPanType[-1:-4, ]

datPanType[, 1] <- sub("-", ".", datPanType[, 1])
datPanType[, 1] <- sub("-", ".", datPanType[, 1])

colnames(datPanType)[1] <- "patient"

library(reshape2)
dat <- merge(datPanType, datPanGF)
colnames(dat)[ncol(dat)] <- "GFPT2"

# Z-score 
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

dat[, ncol(dat)] <- cal_z_score(dat[, ncol(dat)])

table(dat$Subtype)

bl <- mean(dat$GFPT2[which(dat$Subtype == "BRCA_Basal")], na.rm = TRUE)
her2 <- mean(dat$GFPT2[which(dat$Subtype == "BRCA_Her2")], na.rm = TRUE)
LA <- mean(dat$GFPT2[which(dat$Subtype == "BRCA_LumA")], na.rm = TRUE)
LB <- mean(dat$GFPT2[which(dat$Subtype == "BRCA_LumB")], na.rm = TRUE)
nl <- mean(dat$GFPT2[which(dat$Subtype == "BRCA_Normal")], na.rm = TRUE)

dat_mean <- data.frame(names(table(dat$Subtype)), c( bl, her2, LA, LB, nl))
colnames(dat_mean) <- c("Tumor_Type", "GFPT2_Avg.Expression")
