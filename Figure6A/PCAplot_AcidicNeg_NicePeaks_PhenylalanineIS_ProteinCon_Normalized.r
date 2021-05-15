# Copy FJfunctions folder to R working directory.
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentRepeat_16.07.2019")
source("./FJfunctions/loadAll.r")
loadAll()

# Load Data for Acidic metabolomics:
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentRepeat_16.07.2019/QIONG")
tdataA <- loadData("AcidicNeg_Knockdown_Metabolomics_3Scr_18.07.2019.xlsx","TargetLynx",1)
# only 1 character column

# normalization for Acidic metabolomics
for (i in 10:ncol(tdataA)){
  tdataA[, i] <- tdataA[, i]/tdataA[, 2]/tdataA[, 3]*1000000
}

# delete the duplicated scramble samples
n <- seq(1,length(grep("Scr",tdataA$Sample.ID)),2)
tdataA <- tdataA[-c(grep("Scr",tdataA$Sample.ID)[n]), ]

# upload Basic metabolomics data
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentBasic_03.12.2019")
tdataB <- loadData("BasicNeg_Knockdown_Metabolomics_NOV2019.xlsx", "TargetLynx", 1)

# normalization for Basic metabolomics
for (i in 10:ncol(tdataB)){
  tdataB[, i] <- tdataB[, i]/tdataB[, 2]/tdataB[, 3]*1000000
}

# valid peaks for Basic metabolomics
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentBasic_03.12.2019")
library(readxl)
valMetB <- read_excel("TargetedMetabolites-IS Database_Basic_NOV2019.xlsx")
valMetB <- data.frame(valMetB)

ValMetBName <- valMetB$Compound[which(valMetB$Basic_Neg == "OK" & is.na(valMetB$Acidic_Neg) == TRUE)]

# `%!in%` <- Negate(`%in%`)
tdataValB <- tdataB[, which(colnames(tdataB) %in% ValMetBName)]

# combine Acidic and Basic datasets (already normalized)
tdata <- cbind(tdataA, tdataValB)

# delete the Protein con. and IS columns  
tdata <- tdata[, -2:-9]

# add "cell" column for PCA
tdata$Cell <- factor(c(rep("D492", 21), rep("D492M", 21), rep("D492HER2", 21)), 
                       levels = c("D492", "D492M", "D492HER2"))

# add "class" column for PCA
class_1 <- c(rep("wt1", 3), rep("wt2", 3), rep("scr", 3), rep("siGFPT2", 3), rep("siGALE", 3), rep("siPGM2L1", 3), rep("siUGDH", 3))
class <- rep(class_1, 3)
tdata$Class <- factor(class, levels = c("wt1", "wt2", "scr", "siGFPT2", "siGALE", "siPGM2L1", "siUGDH"))

# rearrange the column order
tdata <- tdata[, c(1, ncol(tdata), (ncol(tdata)-1), 2:(ncol(tdata)-2))]

# valid peaks for Acidic&Basic metabolomics
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentRepeat_16.07.2019/QIONG")
library(readxl)
valMet <- read_excel("TargetedMetabolites-IS Database_17.07.2019.xlsx")
valMet <- data.frame(valMet)

ValMetAName <- valMet$Metabolites[which(valMet$Acidic_Neg == "OK")]

ValMetName <- c(ValMetBName, ValMetAName)

# `%!in%` <- Negate(`%in%`)
tdataVal <- tdata[, c(1:3, which(colnames(tdata) %in% ValMetName))]

# delete Scramble the first measurement and wt2
# tdataVal <- tdataVal[-grep("wt2",tdataVal$Sample.ID), ]

# For PCA:
p <- metPCA(tdataVal,
            "Cell",
            minvalue = T, 
            log = T, 
            colors = c("steelblue", "red", "orange1"),
            pointsize = 5,
            CI = TRUE)

library(ggplot2)
library(ggrepel)
p1 <- p +
  # geom_text_repel(aes(label = tdata[, 3]), fontface = "bold") +
  labs(title = "") +
  theme_classic() +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 25, face = "bold"),
        legend.position = "right",
        title = element_text(""),
        axis.title = element_text(size = 30, face = "bold"),
        axis.text = element_text(size = 30, face = "bold", color = "black"),
        axis.text.x.bottom = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y.left = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.line = element_line(size = 2))
p1

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/Combined_Acidic_Basic")
file <- paste0(format(Sys.time(), "%F %H-%M-%S"), " ", "PCAplot.tiff")
ggsave(filename = file, units = "in", dpi = 300, width = 9, height = 6)
