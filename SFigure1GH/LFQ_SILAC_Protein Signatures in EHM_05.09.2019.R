# load the LFQ data from Perseus (after the heatmap with clusters)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/ProteinSignatureInCells/LFQ_ProteinSignatureInCells")
library(readxl)
lfq_e <- read_excel("Annotations for E & M & HER2_LFQ.xlsx", sheet = 1)
lfq_e <- data.frame(lfq_e[, -17:-18])
lfq_e <- lfq_e[which(!is.na(lfq_e$Protein.ID)), ] # 820

lfq_m <- read_excel("Annotations for E & M & HER2_LFQ.xlsx", sheet = 2)
lfq_m <- data.frame(lfq_m[, -17:-18])
lfq_m <- lfq_m[which(!is.na(lfq_m$Protein.ID)), ] # 482

lfq_h <- read_excel("Annotations for E & M & HER2_LFQ.xlsx", sheet = 3)
lfq_h <- data.frame(lfq_h[, -17:-18])
lfq_h <- lfq_h[which(!is.na(lfq_h$Protein.ID)), ] # 537

# load the SILAC data from Perseus (after the heatmap with clusters)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/ProteinSignatureInCells/SILAC_ProteinSignatureInCells")
silac_e <- read_excel("Annotations for E & M & HER2_SILAC.xlsx", sheet = 1)
silac_e <- data.frame(silac_e[, -16:-17]) # 1792

silac_m <- read_excel("Annotations for E & M & HER2_SILAC.xlsx", sheet = 2)
silac_m <- data.frame(silac_m[, -16:-17]) # 701

silac_h <- read_excel("Annotations for E & M & HER2_SILAC.xlsx", sheet = 3)
silac_h <- data.frame(silac_h[, -16:-17]) # 344

silac_h_low <- read_excel("Annotations for SILAC_HER2 LOW.xlsx", sheet = 1)
silac_h_low <- data.frame(silac_h_low) # 108
silac_h_low <- silac_h_low[, c(14:16, 1:12)]

# LFQ data, Genes in both LFQ and SILAC
lfq_silac_e <- lfq_e[which(lfq_e$Gene.Name..GN. %in% silac_e$Gene.Name),] # 332

lfq_silac_m <- lfq_m[which(lfq_m$Gene.Name..GN. %in% silac_m$Gene.Name..GN.),] 
lfq_silac_m <- lfq_silac_m[-which(is.na(lfq_silac_m$Gene.Name..GN.)),] # 103

lfq_silac_h <- lfq_h[which(lfq_h$Gene.Name..GN. %in% silac_h$Gene.Name..GN._R1),]
lfq_silac_h <- lfq_silac_h[-which(is.na(lfq_silac_h$Gene.Name..GN.)),] # 74

lfq_silac_h_low <- lfq_h[which(lfq_h$Gene.Name..GN. %in% silac_h_low$Gene.Name..GN._R1),]
lfq_silac_h_low <- lfq_silac_h_low[-which(is.na(lfq_silac_h_low$Gene.Name..GN.)),] # 21

# SILAC data, Genes in both LFQ and SILAC
silac_lfq_e <- silac_e[which(silac_e$Gene.Name %in% lfq_e$Gene.Name..GN.),] # 334

silac_lfq_m <- silac_m[which(silac_m$Gene.Name..GN. %in% lfq_m$Gene.Name..GN.),] 
silac_lfq_m <- silac_lfq_m[-which(is.na(silac_lfq_m$Gene.Name..GN.)),] # 105

silac_lfq_h <- silac_h[which(silac_h$Gene.Name..GN._R1 %in% lfq_h$Gene.Name..GN.),] 
silac_lfq_h <- silac_lfq_h[-which(is.na(silac_lfq_h$Gene.Name..GN._R1)),] # 74

silac_lfq_h_low <- silac_h_low[which(silac_h_low$Gene.Name..GN._R1 %in% lfq_h$Gene.Name..GN.),]
silac_lfq_h_low <- silac_lfq_h_low[-which(is.na(silac_lfq_h_low$Gene.Name..GN._R1)),] # 21

table(duplicated(silac_lfq_e$Gene.Name)) # 2/334 duplicates
table(duplicated(silac_lfq_m$Gene.Name..GN.)) # 2/105
table(duplicated(silac_lfq_h$Gene.Name..GN._R1)) # 0/74
table(duplicated(silac_lfq_h_low$Gene.Name..GN._R1)) # 0/21

# Delete the duplicated genes after checking they were having similar expression changes 
silac_lfq_e[which(duplicated(silac_lfq_e$Gene.Name)),]
silac_lfq_e[which(silac_lfq_e$Gene.Name == "SEC22B"),]
silac_lfq_e[which(silac_lfq_e$Gene.Name == "FLII"),]
# Delete row 560 and 1384
silac_lfq_e <- silac_lfq_e[-which(rownames(silac_lfq_e) == "560"),]
silac_lfq_e <- silac_lfq_e[-which(rownames(silac_lfq_e) == "1384"),] # 332

silac_lfq_m[which(duplicated(silac_lfq_m$Gene.Name..GN.)),]
silac_lfq_m[which(silac_lfq_m$Gene.Name..GN. == "TMPO"),]
silac_lfq_m[which(silac_lfq_m$Gene.Name..GN. == "FLNA"),]
# Delete row 53 and 425
silac_lfq_m <- silac_lfq_m[-which(rownames(silac_lfq_m) == "53"),]
silac_lfq_m <- silac_lfq_m[-which(rownames(silac_lfq_m) == "425"),] # 103

# Calculate the mean z-score for LFQ data
for (i in 1:length(lfq_silac_e$Protein.ID)) {
  lfq_silac_e$"avg. D492"[i] <- mean(lfq_silac_e$D492_3[i], 
                                     lfq_silac_e$D492_2[i], 
                                     lfq_silac_e$D492_1[i])
  lfq_silac_e$"avg. D492M"[i] <- mean(lfq_silac_e$D492M_2[i], 
                                     lfq_silac_e$D492M_3[i], 
                                     lfq_silac_e$D492M_1[i])
  lfq_silac_e$"avg. D492HER2"[i] <- mean(lfq_silac_e$D492_Her2_3[i], 
                                     lfq_silac_e$D492_Her2_2[i], 
                                     lfq_silac_e$D492_Her2_1[i])
}

for (i in 1:length(lfq_silac_m$Protein.ID)) {
  lfq_silac_m$"avg. D492"[i] <- mean(lfq_silac_m$D492_3[i], 
                                     lfq_silac_m$D492_2[i], 
                                     lfq_silac_m$D492_1[i])
  lfq_silac_m$"avg. D492M"[i] <- mean(lfq_silac_m$D492M_2[i], 
                                      lfq_silac_m$D492M_3[i], 
                                      lfq_silac_m$D492M_1[i])
  lfq_silac_m$"avg. D492HER2"[i] <- mean(lfq_silac_m$D492_Her2_3[i], 
                                         lfq_silac_m$D492_Her2_2[i], 
                                         lfq_silac_m$D492_Her2_1[i])
}

for (i in 1:length(lfq_silac_h$Protein.ID)) {
  lfq_silac_h$"avg. D492"[i] <- mean(lfq_silac_h$D492_3[i], 
                                     lfq_silac_h$D492_2[i], 
                                     lfq_silac_h$D492_1[i])
  lfq_silac_h$"avg. D492M"[i] <- mean(lfq_silac_h$D492M_2[i], 
                                      lfq_silac_h$D492M_3[i], 
                                      lfq_silac_h$D492M_1[i])
  lfq_silac_h$"avg. D492HER2"[i] <- mean(lfq_silac_h$D492_Her2_3[i], 
                                         lfq_silac_h$D492_Her2_2[i], 
                                         lfq_silac_h$D492_Her2_1[i])
}

for (i in 1:length(lfq_silac_h_low$Protein.ID)) {
  lfq_silac_h_low$"avg. D492"[i] <- mean(lfq_silac_h_low$D492_3[i], 
                                     lfq_silac_h_low$D492_2[i], 
                                     lfq_silac_h_low$D492_1[i])
  lfq_silac_h_low$"avg. D492M"[i] <- mean(lfq_silac_h_low$D492M_2[i], 
                                      lfq_silac_h_low$D492M_3[i], 
                                      lfq_silac_h_low$D492M_1[i])
  lfq_silac_h_low$"avg. D492HER2"[i] <- mean(lfq_silac_h_low$D492_Her2_3[i], 
                                         lfq_silac_h_low$D492_Her2_2[i], 
                                         lfq_silac_h_low$D492_Her2_1[i])
}

# delete the redundant columns
lfq_silac_e <- lfq_silac_e[, c(1:3, 17, 16)] # 332
lfq_silac_m <- lfq_silac_m[, c(1:3, 18, 16)] # 103
lfq_silac_h <- lfq_silac_h[, c(1:3, 19, 16)] # 74
lfq_silac_h_low <- lfq_silac_h_low[, c(1:3, 19, 16)] # 21
lfq_silac_h <- rbind(lfq_silac_h, lfq_silac_h_low) # 95

silac_lfq_e <- silac_lfq_e[, c(1:3, 6, 14)] # 332
silac_lfq_m <- silac_lfq_m[, c(1:3, 4)] # 103
silac_lfq_m$"Coefficient.of.variation_D492M" <- rep(0, nrow(silac_lfq_m))
silac_lfq_h <- silac_lfq_h[, c(1:3, 5, 15)] # 74
silac_lfq_h_low <- silac_lfq_h_low[, c(1:3, 5, 15)] # 21
silac_lfq_h <- rbind(silac_lfq_h, silac_lfq_h_low) # 95

# no "NA" in all gene names
table(is.na(lfq_silac_e$Gene.Name..GN.))
table(is.na(silac_lfq_e$Gene.Name))

table(is.na(lfq_silac_m$Gene.Name..GN.))
table(is.na(silac_lfq_m$Gene.Name..GN.))

table(is.na(lfq_silac_h$Gene.Name..GN.))
table(is.na(silac_lfq_h$Gene.Name..GN._R1))

# rename the columsn for merging
colnames(lfq_silac_e) <- c("LFQ_Protein ID", "LFQ_Protein Name", "Gene Name", "LFQ_D492", "ANOVA.FDR_D492")
colnames(silac_lfq_e) <- c("Protein ID", "Protein Name", "Gene Name", "SILAC_D492", "CV_D492")
colnames(lfq_silac_m) <- c("LFQ_Protein ID", "LFQ_Protein Name", "Gene Name", "LFQ_D492M", "ANOVA.FDR_D492M")
colnames(silac_lfq_m) <- c("Protein ID", "Protein Name", "Gene Name", "SILAC_D492M", "CV_D492M")
colnames(lfq_silac_h) <- c("LFQ_Protein ID", "LFQ_Protein Name", "Gene Name", "LFQ_D492HER2", "ANOVA.FDR_D492HER2")
colnames(silac_lfq_h) <- c("Protein ID", "Protein Name", "Gene Name", "SILAC_D492HER2", "CV_D492HER2")
colnames(lfq_silac_h_low) <- c("LFQ_Protein ID", "LFQ_Protein Name", "Gene Name", "LFQ_D492HER2", "ANOVA.FDR_D492HER2")
colnames(silac_lfq_h_low) <- c("Protein ID", "Protein Name", "Gene Name", "SILAC_D492HER2", "CV_D492HER2")

# merge based on the same GENE column
mer_e <- merge(lfq_silac_e, silac_lfq_e)
mer_m <- merge(lfq_silac_m, silac_lfq_m)
mer_h <- merge(lfq_silac_h, silac_lfq_h)

# genes need to change the same (both up- or down-) in LFQ and SILAC
mer_e_final <- mer_e[which((mer_e$LFQ_D492 >= 0 & mer_e$SILAC_D492 >= 0) |
                       (mer_e$LFQ_D492 < 0 & mer_e$SILAC_D492 < 0) ),] # 312

mer_m_final <- mer_m[which((mer_m$LFQ_D492M >= 0 & mer_m$SILAC_D492M >= 0) |
                             (mer_m$LFQ_D492M < 0 & mer_m$SILAC_D492M < 0) ),] # 97

mer_h_final <- mer_h[which((mer_h$LFQ_D492HER2 >= 0 & mer_h$SILAC_D492HER2 >= 0) |
                             (mer_h$LFQ_D492HER2 < 0 & mer_h$SILAC_D492HER2 < 0) ),] # 84

mer_e_final <- mer_e_final[, c(6:7, 1, 4:5, 8:9)]
mer_m_final <- mer_m_final[, c(6:7, 1, 4:5, 8:9)]
mer_h_final <- mer_h_final[, c(6:7, 1, 4:5, 8:9)]

# Export
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/ProteinSignatureInCells")
library(xlsx)
write.xlsx(mer_e_final, "SupplementaryData3_SignatureProteinsInCells.xlsx", sheetName = "D492", row.names = FALSE)
write.xlsx(mer_m_final, "SupplementaryData3_SignatureProteinsInCells.xlsx", sheetName = "D492M", row.names = FALSE, append = TRUE)
write.xlsx(mer_h_final, "SupplementaryData3_SignatureProteinsInCells.xlsx", sheetName = "D492HER2", row.names = FALSE, append = TRUE)

# Load the results from before (did manully in excel)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/ProteinSignatureInCells")
library(readxl)
mer_e_final_before <- read_excel("LFQ_SILAC_ProteinSignaturesInCells.xlsx", sheet = 1)
mer_e_final_before <- data.frame(mer_e_final_before) # 317

mer_m_final_before <- read_excel("LFQ_SILAC_ProteinSignaturesInCells.xlsx", sheet = 2)
mer_m_final_before <- data.frame(mer_m_final_before) # 97

mer_h_final_before <- read_excel("LFQ_SILAC_ProteinSignaturesInCells.xlsx", sheet = 3)
mer_h_final_before <- data.frame(mer_h_final_before) # 87

`%!in%` <- Negate(`%in%`)
mer_e_final_before[which(mer_e_final_before$Gene.Name..GN. %!in% mer_e_final$`Gene Name`),]
# AP3B1, PPP1R12A, SAP18, SNRPE, WDR6

mer_m_final_before[which(mer_m_final_before$Gene.Name..GN. %!in% mer_m_final$`Gene Name`),]
# 0

mer_m_final[which(mer_m_final$`Gene Name` %!in% mer_m_final_before$Gene.Name..GN...3),]
# CNOT2

mer_h_final_before[which(mer_h_final_before$Gene.Name..GN. %!in% mer_h_final$`Gene Name`),]
# TBL3, TUBG1, YLPM1, TMED7-TICAM2

#----------------------------------------------------------------------------------#

                             # Heatmap (LFQ and SILAC)

#----------------------------------------------------------------------------------#
# Heatmap for LFQ
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/ProteinSignatureInCells/LFQ_ProteinSignatureInCells")
Hlfq <- read.delim("PerseusOutPutTable_LFQ_ForHeatmapPlottingInR_12.02.2019.txt")
Hlfq <- data.frame(Hlfq)
Hlfq <- Hlfq[-1:-2, ]
Hlfq <- Hlfq[, c(14, 1:9)]
colnames(Hlfq)[8:10] <- paste0("D492HER2", "_", 1:3)

for (i in 2:10){
  Hlfq[, i] <- as.numeric(levels(Hlfq[, i]))[Hlfq[, i]]
}

table(duplicated(Hlfq$Protein.ID))

rownames(Hlfq) <- Hlfq[, 1]
Hlfq <- Hlfq[, -1]

library(ComplexHeatmap)
library(grDevices)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/ProteinSignatureInCells/LFQ_ProteinSignatureInCells")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "HeatmapPlot.tiff")
tiff(filename = file, units="in", width = 7, height = 10, res = 300)

max(Hlfq) # 2.408921
min(Hlfq) # -2.530834
col_fun = circlize::colorRamp2(c(-2.53, 0, 2.41), c("steelblue", "white", "red"))

Hlfq <- as.matrix(Hlfq)

ht <- Heatmap(Hlfq,
              col_fun,
              cluster_columns = FALSE, # keep the column order as in the matrix
              heatmap_legend_param = list(title = "Expression",
                                          title_gp = gpar(fontsize = 25),
                                          labels_gp = gpar(fontsize = 25),
                                          at = c(-2.53, 0, 2.41),
                                          labels = c("low", "zero", "high"),
                                          legend_height = unit(4, "cm")),
              column_names_side = "top",
              show_row_names = FALSE,
              show_row_dend = FALSE,
              column_names_gp = gpar(cex = 2),
              row_names_gp = gpar(cex = 2),
              column_dend_height = unit(3, "cm"),
              km = 6,
             # combined_name_fun = NULL, 
              row_title = c("C1", "C2", "C3", "C4", "C5", "C6"),
              row_title_gp = gpar(fontsize = 25),
              column_title = "LFQ",
              column_title_gp = gpar(fontsize = 25))
print(ht)
# export the heatmap plot into .tiff
dev.off()

#-------------------------------------------------------------------------------------#
# Heatmap for SILAC
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/ProteinSignatureInCells/SILAC_ProteinSignatureInCells")
Hsilac <- read.delim("PerseusOutPutTable_SILAC_ForHeatmapPlottingInR_12.02.2019.txt")
Hsilac <- data.frame(Hsilac)
Hsilac <- Hsilac[-1, ]
Hsilac <- Hsilac[, c(12, 1, 3, 2)]
colnames(Hsilac)[2:4] <- c("D492", "D492M", "D492HER2")

for (i in 2:4){
  Hsilac[, i] <- as.numeric(levels(Hsilac[, i]))[Hsilac[, i]]
}

table(duplicated(Hsilac$Protein.ID_total))

rownames(Hsilac) <- Hsilac[, 1]
Hsilac <- Hsilac[, -1]

library(ComplexHeatmap)
library(grDevices)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/ProteinSignatureInCells/SILAC_ProteinSignatureInCells")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "HeatmapPlot.tiff")
tiff(filename = file, units="in", width = 4, height = 10, res = 300)

max(Hsilac) # 1.1547
min(Hsilac) # -1.154701
col_fun = circlize::colorRamp2(c(-1.155, 0, 1.155), c("steelblue", "white", "red"))

Hsilac <- as.matrix(Hsilac)

ht <- Heatmap(Hsilac,
              col_fun,
              cluster_columns = F, # keep the column order as in the matrix
              heatmap_legend_param = list(title = "Expression",
                                          title_gp = gpar(fontsize = 25),
                                          labels_gp = gpar(fontsize = 25),
                                          at = c(-1.155, 0, 1.155),
                                          labels = c("low", "zero", "high"),
                                          legend_height = unit(4, "cm")),
              column_names_side = "top",
              #row_title_side = "right",
              show_row_names = FALSE,
              show_row_dend = FALSE,
              column_names_gp = gpar(cex = 2),
              row_names_gp = gpar(cex = 2),
              column_dend_height = unit(4, "cm"),
              km = 6,
              # combined_name_fun = NULL, 
              row_title = c("C1", "C2", "C3", "C4", "C5", "C6"),
              row_title_gp = gpar(fontsize = 25),
              column_title = "SILAC",
              column_title_gp = gpar(fontsize = 25),
              column_order = sort(colnames(Hsilac), decreasing = FALSE))
print(ht)
# export the heatmap plot into .tiff
dev.off()

#---------------------------------------------------------------------------#



#---------------------------------------------------------------------------#
# Load the results from before (did manully in excel)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/ProteinSignatureInCells")
library(readxl)
mer_e_final_before <- read_excel("LFQ_SILAC_ProteinSignaturesInCells.xlsx", sheet = 1)
mer_e_final_before <- data.frame(mer_e_final_before) # 317

mer_m_final_before <- read_excel("LFQ_SILAC_ProteinSignaturesInCells.xlsx", sheet = 2)
mer_m_final_before <- data.frame(mer_m_final_before) # 97

mer_h_final_before <- read_excel("LFQ_SILAC_ProteinSignaturesInCells.xlsx", sheet = 3)
mer_h_final_before <- data.frame(mer_h_final_before) # 87

# LFQ
setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018")
lfq <- read_excel("Dundee_LFQ_valid_Log2_imputated_p.value_t.stat_q.value_mean_22.11.2018.xlsx")
lfq <- data.frame(lfq)

lfqE <- lfq[which(lfq$Gene.Name..GN. %in% mer_e_final_before$Gene.Name..GN.), ] # 318
lfqM <- lfq[which(lfq$Gene.Name..GN. %in% mer_m_final_before$Gene.Name..GN...3), ] # 98
lfqH <- lfq[which(lfq$Gene.Name..GN. %in% mer_h_final_before$Gene.Name..GN.), ] # 87

# Delete the duplicated ones
lfqE <- lfqE[-168, ] # 317
lfqM <- lfqM[-98, ] # 97

# delete the useless columns
lfqE <- lfqE[, c(3, 18, 15, 24, 21)]
lfqM <- lfqM[, c(3, 18, 15, 30, 27)]
lfqH <- lfqH[, c(3, 24, 21, 30, 27)]

# SILAC
silacEM <- read_excel("Dundee_proteo_SILAC_log2_imput_t-test_median_valid_22.11.2018.xlsx", sheet = 1)
silacEM <- data.frame(silacEM)
silacHE <- read_excel("Dundee_proteo_SILAC_log2_imput_t-test_median_valid_22.11.2018.xlsx", sheet = 2)
silacHE <- data.frame(silacHE)
silacHM <- read_excel("Dundee_proteo_SILAC_log2_imput_t-test_median_valid_22.11.2018.xlsx", sheet = 3)
silacHM <- data.frame(silacHM)

# D492
silacE_EM <- silacEM[which(silacEM$Gene.Name %in% mer_e_final_before$Gene.Name..GN.), ]
silacE_EM <- silacE_EM[, c(1:3, 7, 8)]
colnames(silacE_EM)[5] <- "EM_p.value"
silacE_HE <- silacHE[which(silacHE$Gene.Name %in% mer_e_final_before$Gene.Name..GN.), ]
silacE_HE <- silacE_HE[, c(1:3, 7, 8)]
colnames(silacE_HE)[5] <- "HE_p.value"

silacE <- merge(silacE_EM, silacE_HE)
silacE <- silacE[-which(duplicated(silacE$Gene.Name)), ] # 317

# D492M
silacM_EM <- silacEM[which(silacEM$Gene.Name %in% mer_m_final_before$Gene.Name..GN...3), ]
silacM_EM <- silacM_EM[, c(1:3, 7, 8)]
colnames(silacM_EM)[5] <- "EM_p.value"
silacM_HM <- silacHM[which(silacHM$Gene.Name %in% mer_m_final_before$Gene.Name..GN...3), ]
silacM_HM <- silacM_HM[, c(1:3, 7, 8)]
colnames(silacM_HM)[5] <- "HM_p.value"

silacM <- merge(silacM_EM, silacM_HM)
silacM <- silacM[-which(duplicated(silacM$Gene.Name)), ] # 97

# D492HER2
silacH_HE <- silacHE[which(silacHE$Gene.Name %in% mer_h_final_before$Gene.Name..GN.), ]
silacH_HE <- silacH_HE[, c(1:3, 7, 8)]
colnames(silacH_HE)[5] <- "HE_p.value"
silacH_HM <- silacHM[which(silacHM$Gene.Name %in% mer_h_final_before$Gene.Name..GN.), ]
silacH_HM <- silacH_HM[, c(1:3, 7, 8)]
colnames(silacH_HM)[5] <- "HM_p.value"

silacH <- merge(silacH_HE, silacH_HM)
silacH <- silacH[-which(duplicated(silacH$Gene.Name)), ] # 87

# change column name for LFQ
colnames(lfqE)[1] <- colnames(silacE)[3]
colnames(lfqM)[1] <- colnames(silacM)[3]
colnames(lfqH)[1] <- colnames(silacH)[3]

# merge SILAC and LFQ
lfq_silac_e <- merge(lfqE, silacE)
lfq_silac_m <- merge(lfqM, silacM)
lfq_silac_h <- merge(lfqH, silacH)

# rename 
lfq_silac_e <- lfq_silac_e[, c(6:7, 1, 2:3, 8:9, 4:5, 10:11)]
lfq_silac_m <- lfq_silac_m[, c(6:7, 1, 2:3, 8:9, 4:5, 10:11)]
lfq_silac_h <- lfq_silac_h[, c(6:7, 1, 2:3, 8:9, 4:5, 10:11)]

lfq_silac_e$"Mean_Log2.M.E." <- -lfq_silac_e$Mean_Log2.E.M. # D492M/D492
lfq_silac_e$"Median_Log2.M.E." <- -lfq_silac_e$Median_Log2.E.M. # D492M/D492
lfq_silac_e <- lfq_silac_e[, c(-4, -6)]
lfq_silac_e <- lfq_silac_e[, c(1:3, 10, 4, 11, 5:9)]

lfq_silac_m$"Mean_Log2.M.E." <- -lfq_silac_m$Mean_Log2.E.M. # D492M/D492
lfq_silac_m$"Median_Log2.M.E." <- -lfq_silac_m$Median_Log2.E.M. # D492M/D492
lfq_silac_m <- lfq_silac_m[, c(-4, -6)]
lfq_silac_m <- lfq_silac_m[, c(1:3, 10, 4, 11, 5:9)]

colnames(lfq_silac_e) <- c("Protein ID", 
                           "Protein Names", 
                           "Gene Names", 
                           "LFQ_Log2(D492M/D492)",
                           "LFQ_FDR_D492M/D492",
                           "SILAC_Log2(D492M/D492)",
                           "SILAC_pvalue_D492M/D492",
                           "LFQ_Log2(D492HER2/D492)",
                           "LFQ_FDR_D492HER2/D492",
                           "SILAC_Log2(D492HER2/D492)",
                           "SILAC_pvalue_D492HER2/D492")

colnames(lfq_silac_m) <- c("Protein ID", 
                           "Protein Names", 
                           "Gene Names", 
                           "LFQ_Log2(D492M/D492)",
                           "LFQ_FDR_D492M/D492",
                           "SILAC_Log2(D492M/D492)",
                           "SILAC_pvalue_D492M/D492",
                           "LFQ_Log2(D492HER2/D492M)",
                           "LFQ_FDR_D492HER2/D492M",
                           "SILAC_Log2(D492HER2/D492M)",
                           "SILAC_pvalue_D492HER2/D492M")

colnames(lfq_silac_h) <- c("Protein ID", 
                           "Protein Names", 
                           "Gene Names", 
                           "LFQ_Log2(D492HER2/D492)",
                           "LFQ_FDR_D492HER2/D492",
                           "SILAC_Log2(D492HER2/D492)",
                           "SILAC_pvalue_D492HER2/D492",
                           "LFQ_Log2(D492HER2/D492M)",
                           "LFQ_FDR_D492HER2/D492M",
                           "SILAC_Log2(D492HER2/D492M)",
                           "SILAC_pvalue_D492HER2/D492M")

# export
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/ProteinSignatureInCells")
library(xlsx)
write.xlsx(lfq_silac_e, "SupplementaryData3_SignatureProteinsInCells.xlsx",
           row.names = FALSE, sheetName = "D492")
write.xlsx(lfq_silac_m, "SupplementaryData3_SignatureProteinsInCells.xlsx",
           row.names = FALSE, sheetName = "D492M", append = TRUE)
write.xlsx(lfq_silac_h, "SupplementaryData3_SignatureProteinsInCells.xlsx",
           row.names = FALSE, sheetName = "D492HER2", append = TRUE)
