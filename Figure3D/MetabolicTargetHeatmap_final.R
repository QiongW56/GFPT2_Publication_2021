# metabolic targets
setwd("F:/DundeeAnalysisResults/All datasets_organized")
met <- read.delim("MetabolicTargetsFromProteomicDataAnalysis.txt")
met <- data.frame(met)
met <- unique(met)
met$GENE <- as.character.factor(met$GENE)

# LFQ
setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018")
library(readxl)
lfq <- read_excel("Dundee_LFQ_valid_Log2_imputated_p.value_t.stat_q.value_mean_22.11.2018.xlsx")
lfq <- data.frame(lfq)
lfq <- lfq[, c(3, 10:12, 4:6, 7:9)]
lfq_met <- lfq[which(lfq$Gene.Name..GN. %in% met$GENE), ]
rownames(lfq_met) <- lfq_met$Gene.Name..GN.
lfq_met <- lfq_met[order(lfq_met$Gene.Name..GN., decreasing = TRUE), ]

# SILAC
setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018")
silacEM <- read_excel("Dundee_proteo_SILAC_log2_imput_t-test_median_valid_22.11.2018.xlsx", sheet = 1)
silacEM <- data.frame(silacEM)
silacEM$"Ratio.L.M.TotProt.R1.LOG2" <- -silacEM$Ratio.M.L.TotProt.R1..LOG2
silacEM$"Ratio.L.M.TotProt.R2.LOG2" <- -silacEM$Ratio.M.L.TotProt.R2.LOG2
silacEM$"Ratio.L.M.TotProt.R3.LOG2" <- -silacEM$Ratio.M.L.TotProt.R3.LOG2
for (i in 1:nrow(silacEM)){
  silacEM$"D492M"[i] <- median(as.numeric(silacEM[i, 10:12]))
}
silacEM <- silacEM[, c(3, 13)]

silacHE <- read_excel("Dundee_proteo_SILAC_log2_imput_t-test_median_valid_22.11.2018.xlsx", sheet = 2)
silacHE <- data.frame(silacHE)
for (i in 1:nrow(silacHE)){
  silacHE$"D492HER2"[i] <- median(as.numeric(silacHE[i, 4:6]))
}
silacHE <- silacHE[, c(3, 10)]

colnames(silacEM)[1] <- c("Gene.Name")
colnames(silacHE)[1] <- c("Gene.Name")

silacEM_met <- silacEM[which(silacEM$Gene.Name %in% met$GENE), ]
silacHE_met <- silacHE[which(silacHE$Gene.Name %in% met$GENE), ]

# merge SILAC
silac_met <- merge(silacEM_met, silacHE_met)
silac_met$"D492" <- rep(0, nrow(silac_met))
rownames(silac_met) <- silac_met$Gene.Name
silac_met <- silac_met[order(silac_met$Gene.Name, decreasing = TRUE), ]

# z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

lfq_met[,-1] <- t(apply(lfq_met[,-1], 1, cal_z_score))
lfq_met_matrix <- as.matrix(lfq_met[,-1])

silac_met[,-1] <- t(apply(silac_met[,-1], 1, cal_z_score))
silac_met_matrix <- as.matrix(silac_met[,-1])

# Heatmap
library(ComplexHeatmap)
library(grDevices)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/MetabolicTargetsHeatmap")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "HeatmapPlot.tiff")
tiff(filename = file, units="in", width = 12, height = 12, res = 300)

col_fun_1 = circlize::colorRamp2(c(-1.83, 0, 1.70), c("steelblue", "white", "red"))
ht_1 <- Heatmap(lfq_met_matrix, 
                col = col_fun_1,
                heatmap_legend_param = list(title = "Expression",
                                            title_gp = gpar(fontsize = 25),
                                            labels_gp = gpar(fontsize = 25),
                                            at = c(-1.83, 0, 1.70),
                                            labels = c("low", "zero", "high"),
                                            legend_height = unit(4, "cm")),
                column_names_side = "top", 
                row_names_side = "left",
                row_dend_side = "left",
                column_names_gp = gpar(cex = 2),
                row_names_gp = gpar(cex = 1.6),
                row_dend_width = unit(6,"cm"),
                column_dend_height = unit(3, "cm"),
                km = 6,
                combined_name_fun = NULL,
                row_title = paste0("C", 1:6),
                row_title_gp = gpar(fontsize = 25),
                column_title = "LFQ",
                column_title_gp = gpar(fontsize = 25))

col_fun_2 = circlize::colorRamp2(c(-1.16, 0, 1.16), c("steelblue", "white", "red"))
ht_2 <- Heatmap(silac_met_matrix, 
                col = col_fun_2,
                heatmap_legend_param = list(title = "Expression",
                                            title_gp = gpar(fontsize = 25),
                                            labels_gp = gpar(fontsize = 25),
                                            at = c(-1.16, 0, 1.16),
                                            labels = c("low", "zero", "high"),
                                            legend_height = unit(4, "cm")),
                column_names_side = "top", 
                row_names_side = "left",
                row_dend_side = "left",
                column_names_gp = gpar(cex = 2),
                row_names_gp = gpar(cex = 2),
                row_dend_width = unit(6,"cm"),
                column_dend_height = unit(3, "cm"),
                km = 6,
                combined_name_fun = NULL,
                row_title_gp = gpar(fontsize = 25),
                column_title = "SILAC",
                column_title_gp = gpar(fontsize = 25),
                show_heatmap_legend = F,
                show_row_names = F)

ht_list <- ht_1 + ht_2
draw(ht_list)

dev.off() # export the heatmap plot into .tiff

#----------------------------------------------------------------------# 
# check the pvalues for SILAC, all SILAC p values are < 0.05
# metabolic targets
setwd("F:/DundeeAnalysisResults/All datasets_organized")
met <- read.delim("MetabolicTargetsFromProteomicDataAnalysis.txt")
met <- data.frame(met)
met <- unique(met)
met$GENE <- as.character.factor(met$GENE)

setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018")
library(readxl)
lfq <- read_excel("Dundee_LFQ_valid_Log2_imputated_p.value_t.stat_q.value_mean_22.11.2018.xlsx")
lfq <- data.frame(lfq)
lfq_met <- lfq[which(lfq$Gene.Name..GN. %in% met$GENE), ]

# SILAC
setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018")
silacEM <- read_excel("Dundee_proteo_SILAC_log2_imput_t-test_median_valid_22.11.2018.xlsx", sheet = 1)
silacEM <- data.frame(silacEM)
silacEM_met <- silacEM[which(silacEM$Gene.Name %in% met$GENE), ]

silacHE <- read_excel("Dundee_proteo_SILAC_log2_imput_t-test_median_valid_22.11.2018.xlsx", sheet = 2)
silacHE <- data.frame(silacHE)
silacHE_met <- silacHE[which(silacHE$Gene.Name %in% met$GENE), ]

silacHM <- read_excel("Dundee_proteo_SILAC_log2_imput_t-test_median_valid_22.11.2018.xlsx", sheet = 3)
silacHM <- data.frame(silacHM)
silacHM_met <- silacHM[which(silacHM$Gene.Name %in% met$GENE), ]
