# this file is to find the EMT metabolic makers from paper "Dihydropyrimidine accumulation is required for the epithelial-mesenchymal transition" in our datasets 
# Upload the markers from the paper
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/EMTMetaboliteMarker_heatmap")
MarkerMet <- data.frame(read_excel("Markers.xlsx", col_names = FALSE))
colnames(MarkerMet) <- "GeneName"

# Upload LFQ data
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/1 DundeeProteomicDataSet_14112017/ProteomicDataAnalysisV3")
library(readxl)
datLFQ <- read_excel("Proteomics_LFQ_summary_16.04.2018.xlsx")
datLFQ <- data.frame(datLFQ)
datLFQ <- datLFQ[, c(5, 8:16)]
datLFQ[, 2:10] <- log2(datLFQ[, 2:10])

# find the markers in LFQ
datLFQ <- datLFQ[which(datLFQ$Gene.Name..GN. %in% MarkerMet$GeneName), ]
names(datLFQ) <- c("GeneName", 
                   paste0("D492", "_", 1:3), 
                   paste0("D492M", "_", 1:3), 
                   paste0("D492HER2", "_", 1:3))

is.na(datLFQ[, 2:10]) <- sapply(datLFQ[, 2:10], is.infinite)

for (i in 1:nrow(datLFQ)){
  for (j in 2:ncol(datLFQ)){
   if (is.na(datLFQ[i, j]) == TRUE){
     datLFQ[i, j] <- min(datLFQ[, 2:10], na.rm = T)
   }
  } # is.na(dat) <- sapply(dat, is.infinite)
}

# SILAC
funSILAC <- function(sheetNO){
  # Upload SILAC data
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/2 SILIC_DundeePhosphoproteomicsDataset/Results/Proteomics/Analysis Results - OR _QIONG_Organized")
  datSILAC <- data.frame(read_excel("1 Proteomics_SILAC_summary_11.10.2018.xlsx", sheet = sheetNO))
  datSILAC <- datSILAC[, c(3, grep("LOG2", names(datSILAC)))]
  
  for (i in 2:4){
    datSILAC[, i] <- as.numeric(datSILAC[, i])
  }
  
  # Find the markers in SILAC
  datSILAC <- datSILAC[which(datSILAC$Gene.Name %in% MarkerMet$GeneName), ]
  for (i in 1:nrow(datSILAC)){
    if (rowSums(!is.na(datSILAC[i, 2:4])) >= 2){
      datSILAC$Group[i] <- "valid"
    } else {
      datSILAC$Group[i] <- "invalid"
    }
  }
  
  datSILAC <- datSILAC[-which(datSILAC$Group == "invalid"), ]
  
  datSILAC <- datSILAC[, -ncol(datSILAC)]
  
  for (i in 1:nrow(datSILAC)){
    datSILAC$Mean[i] <- mean(as.numeric(datSILAC[i, 2:4]), na.rm = TRUE)
  }
  
  datSILAC <- datSILAC[, c(1, ncol(datSILAC))]
  
  if (sheetNO == 2){
    colnames(datSILAC) <- c("GeneName", "D492") 
  } else if (sheetNO == 3){
    colnames(datSILAC) <- c("GeneName", "D492HER2")
  }
  
  return(datSILAC)
  
}

silacEM <- funSILAC(2)
silacHM <- funSILAC(3)

library(reshape2)
datSILAC <- merge(silacEM, silacHM)
datSILAC$D492M <- rep(0, nrow(datSILAC))
datSILAC <- datSILAC[, c(1:2, 4, 3)]

# both in LFQ and in SILAC
datLFQ <- datLFQ[which(datLFQ$GeneName %in% datSILAC$GeneName), ]
datSILAC <- datSILAC[which(datSILAC$GeneName %in% datLFQ$GeneName), ]

datLFQ <- datLFQ[order(datLFQ$GeneName), ]

#------------------------------------------------------------------------------------------#

# Heatmap

#------------------------------------------------------------------------------------------#
# z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

datLFQ[,-1] <- t(apply(datLFQ[,-1], 1, cal_z_score))
datSILAC[,-1] <- t(apply(datSILAC[,-1], 1, cal_z_score))

# plotting
rownames(datLFQ) <- datLFQ[, 1]
datLFQ <- datLFQ[, -1]

rownames(datSILAC) <- datSILAC[, 1]
datSILAC <- datSILAC[, -1]

#------------------------------------------------------------------------------------------#

# source the heatmap R files

#------------------------------------------------------------------------------------------#
# Run heatmap R file in the same directory
datLFQ <- as.matrix(datLFQ)
datSILAC <- as.matrix(datSILAC)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/EMTMetaboliteMarker_heatmap/Rfiles")
source("heatmap-plotting.R")

#------------------------------------------------------------------------------------------#

# LFQ_Heatmap (not used)

#------------------------------------------------------------------------------------------#
max(datLFQ) # 1.87289
min(datLFQ) #  -1.7566
col_fun_1 = circlize::colorRamp2(c(-1.88, 0, 1.77), c("steelblue", "white", "red"))

datLFQ <- as.matrix(datLFQ)

library(ComplexHeatmap)
library(grDevices)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/EMTMetaboliteMarker_heatmap/Figures")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "HeatmapPlot.tiff")
tiff(filename = file, units="in", width = 7, height = 10, res = 300)
ht <- Heatmap(datLFQ,
              col_fun_1,
              cluster_columns = FALSE, # keep the column order as in the matrix
              cluster_rows = FALSE,
              heatmap_legend_param = list(title = "Expression",
                                          title_gp = gpar(fontsize = 20),
                                          labels_gp = gpar(fontsize = 20),
                                          at = c(-1.88, 0, 1.77),
                                          labels = c("low", "zero", "high"),
                                          legend_height = unit(4, "cm")),
              column_names_side = "top",
              row_names_side = "left",
              show_row_names = TRUE,
              show_row_dend = FALSE,
              show_column_dend = TRUE,
              column_names_gp = gpar(cex = 2),
              row_names_gp = gpar(cex = 2),
              column_dend_height = unit(3, "cm"),
              km = 1,
              # combined_name_fun = NULL, 
              # row_title = c("C1", "C2", "C3", "C4", "C5", "C6"),
              row_title_gp = gpar(fontsize = 20),
              column_title = "LFQ",
              column_title_gp = gpar(fontsize = 40),
              column_order = sort(colnames(datLFQ), decreasing = FALSE))
print(ht)

# export the heatmap plot into .tiff
dev.off()

#------------------------------------------------------------------------------------------#

# SILAC_Heatmap (not used)

#------------------------------------------------------------------------------------------#
max(datSILAC) # 1.14932
min(datSILAC) # -1.151205
col_fun_2 = circlize::colorRamp2(c(-1.15, 0, 1.16), c("steelblue", "white", "red"))

datSILAC <- as.matrix(datSILAC)

library(ComplexHeatmap)
library(grDevices)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/EMTMetaboliteMarker_heatmap/Figures")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "HeatmapPlot.tiff")
tiff(filename = file, units="in", width = 5, height = 10, res = 300)
ht <- Heatmap(datSILAC,
              col_fun_2,
              cluster_columns = F, # keep the column order as in the matrix
              cluster_rows = FALSE,
              heatmap_legend_param = list(title = "Expression",
                                          title_gp = gpar(fontsize = 25),
                                          labels_gp = gpar(fontsize = 25),
                                          at = c(-1.15, 0, 1.16),
                                          labels = c("low", "zero", "high"),
                                          legend_height = unit(4, "cm")),
              column_names_side = "top",
              row_names_side = "left",
              show_row_names = TRUE,
              show_row_dend = FALSE,
              show_column_dend = TRUE,
              column_names_gp = gpar(cex = 2),
              row_names_gp = gpar(cex = 2),
              column_dend_height = unit(3, "cm"),
              km = 1,
              # combined_name_fun = NULL, 
              #  row_title = c("C1", "C2", "C3", "C4", "C5", "C6"),
              row_title_gp = gpar(fontsize = 25),
              column_title = "SILAC",
              column_title_gp = gpar(fontsize = 40),
              column_order = sort(colnames(datSILAC), decreasing = FALSE))
print(ht)

# export the heatmap plot into .tiff
dev.off()

#------------------------------------------------------------------------------------------#

# Export (not used)

#------------------------------------------------------------------------------------------#
# export
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/EMTMetaboliteMarker_heatmap/OutputData")
library(xlsx)
write.xlsx(datLFQ, "Metabolic_EMT_Markers_LFQ.xlsx", sheetName = "LFQ", row.names = FALSE)
write.xlsx(datSILAC, "Metabolic_EMT_Markers_SILAC.xlsx", sheetName = "SILAC", row.names = FALSE)

setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018/QIONG_Functions")
load("heatmap.plot.Rdata")

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/EMTMetaboliteMarker_heatmap/Figures")
heatmap(datLFQ, datSILAC, km = 1, title1 = "LFQ", title2 = "SILAC")
