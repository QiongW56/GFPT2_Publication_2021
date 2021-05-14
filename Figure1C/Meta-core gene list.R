# Load EMT gene list from literature
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Figures/Figure 1/1A")
metaEMT <- read.delim("Meta-core gene list.txt", col.names = "GENE")
metaEMT <- data.frame(metaEMT) # 129

# LFQ with EM and HE
setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018")
library(readxl)
lfq <- read_excel("Dundee_LFQ_valid_Log2_imputated_p.value_t.stat_q.value_mean_22.11.2018.xlsx")
lfq <- data.frame(lfq)

# SILAC with EM
silacEM <- read_excel("Dundee_proteo_SILAC_log2_imput_t-test_median_valid_22.11.2018.xlsx", sheet = 1)
silacEM <- data.frame(silacEM)

# SILAC with HE
silacHE <- read_excel("Dundee_proteo_SILAC_log2_imput_t-test_median_valid_22.11.2018.xlsx", sheet = 2)
silacHE <- data.frame(silacHE)

# genes in the literature also in LFQ and SILAC
meta_lfq <- lfq[which(lfq$Gene.Name..GN. %in% metaEMT$GENE), ] # 32
meta_silacEM <- silacEM[which(silacEM$Gene.Name %in% metaEMT$GENE), ] # 37
meta_silacHE <- silacHE[which(silacHE$Gene.Name %in% metaEMT$GENE), ] # 37

# get useful columns, gene name, pvalue, FDR and FCs for LFQ (EM and HE)
meta_lfq <- meta_lfq[, c(3, 14:15, 18, 20:21, 24)] # 32

# get useful columns, gene name, pvalue and FCs for SILAC
meta_silacEM <- meta_silacEM[, c(3, 8, 7)] # 37
colnames(meta_silacEM)[2] <- "SILAC_EM_p.value"
meta_silacHE <- meta_silacHE[, c(3, 8, 7)] # 37
colnames(meta_silacHE)[2] <- "SILAC_HE_p.value"

# merge two datasets into one
meta_silac <- merge(meta_silacEM, meta_silacHE) # 41

# in both SILAC and LFQ
meta_lfq_s <- meta_lfq[which(meta_lfq$Gene.Name..GN. %in% meta_silac$Gene.Name), ]
meta_silac_l <- meta_silac[which(meta_silac$Gene.Name %in% meta_lfq$Gene.Name..GN.), ]

# Delete the duplicated genes
meta_lfq_s <- meta_lfq_s[-which(duplicated(meta_lfq_s$Gene.Name..GN.)), ] # 26/27

meta_silac_l <- meta_silac_l[-c(29, 31, 32), ] # there were 4 TPMs and deleted 3 of them
meta_silac_l <- meta_silac_l[-which(duplicated(meta_silac_l$Gene.Name)), ] # 26/29

# change columns to merge SILAC and LFQ
colnames(meta_lfq_s)[1] <- "Gene.Name"
colnames(meta_lfq_s)[2] <- "LFQ_EM_p.value"
colnames(meta_lfq_s)[5] <- "LFQ_HE_p.value"
colnames(meta_silac_l)[1] <- "Gene.Name"

s_l <- merge(meta_lfq_s, meta_silac_l) # 26

# deleted "ALDH1A3", "EPCAM", "MAP1B", "MAP7", the LFQ FDR of all these genes are more than 0.05
# s_l  <- s_l[-grep("ALDH1A3", s_l$Gene.Name), ]
# s_l  <- s_l[-grep("EPCAM", s_l$Gene.Name), ]
# s_l  <- s_l[-grep("MAP1B", s_l$Gene.Name), ]
# s_l  <- s_l[-grep("MAP7", s_l$Gene.Name), ]

# note, PRKCA and FGF2, the LFQ FDR are also more than 0.05, but they are still in the final heatmap plot

s_l <- s_l[, c(1, 4, 9, 7, 11, 2:3, 5:6, 8, 10)]


#---------------------------------------------------------------------------------------#

                                   # generate heatmap

#---------------------------------------------------------------------------------------#
# generate Heatmap
gene <- data.frame(s_l$Gene.Name)

# load data
setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018")
library(readxl)
lfq <- read_excel("Dundee_LFQ_valid_Log2_imputated_p.value_t.stat_q.value_mean_22.11.2018.xlsx")
lfq <- data.frame(lfq)

# SILAC with EM
silacEM <- read_excel("Dundee_proteo_SILAC_log2_imput_t-test_median_valid_22.11.2018.xlsx", sheet = 1)
silacEM <- data.frame(silacEM)

# SILAC with HE
silacHE <- read_excel("Dundee_proteo_SILAC_log2_imput_t-test_median_valid_22.11.2018.xlsx", sheet = 2)
silacHE <- data.frame(silacHE)

# both in LFQ and "gene"
lfq <- lfq[which(lfq$Gene.Name..GN. %in% gene$s_l.Gene.Name), ]
lfq <- lfq[-which(duplicated(lfq$Gene.Name..GN.)), ]
lfq <- lfq[, c(3, 10:12, 4:6, 7:9)] # D492, D492M, D492HER2, already Log2 values
lfq <- lfq[order(lfq$Gene.Name..GN.), ] # alphatical order

# both in SILAC and "gene"
silacEM <- silacEM[which(silacEM$Gene.Name %in% gene$s_l.Gene.Name), ]
silacHE <- silacHE[which(silacHE$Gene.Name %in% gene$s_l.Gene.Name), ]

silacEM <- silacEM[-which(duplicated(silacEM$Gene.Name)), ]
silacHE <- silacHE[-which(duplicated(silacHE$Gene.Name)), ]

silacEM <- silacEM[, c(3, 7)]
silacEM$Median_Log2.E.M. <- -silacEM$Median_Log2.E.M. # transform to D492M/D492
silacHE <- silacHE[, c(3, 7)]

silac <- merge(silacEM, silacHE)
silac$"D492" <- rep(0, nrow(silac))
colnames(silac)[2] <- "D492M"
colnames(silac)[3] <- "D492HER2"
silac <- silac[, c(1, 4, 2, 3)]

# heatmap plotting
dat1 <- lfq
dat2 <- silac
TargetList <- gene

rownames(dat1) <- dat1[, 1]
rownames(dat2) <- dat2[, 1]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

dat1[,-1] <- t(apply(dat1[,-1], 1, cal_z_score))
dat1_matrix <- as.matrix(dat1[,-1])

dat2[,-1] <- t(apply(dat2[,-1], 1, cal_z_score))
dat2_matrix <- as.matrix(dat2[,-1])

library(ComplexHeatmap)
library(grDevices)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Figures/Figure 1/1A")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "HeatmapPlot.tiff")
tiff(filename = file, units = "in", width = 12, height = 12, res = 300)

col_fun_1 = circlize::colorRamp2(c(-2.2, 0, 1.75), c("steelblue", "white", "red"))
ht_1 <- Heatmap(dat1_matrix, 
                col = col_fun_1,
                heatmap_legend_param = list(title = "Expression",
                                            title_gp = gpar(fontsize = 25),
                                            labels_gp = gpar(fontsize = 25),
                                            at = c(-2.2, 0, 1.75),
                                            labels = c("low", "zero", "high"),
                                            legend_height = unit(4, "cm")),
                column_names_side = "top", 
                row_names_side = "left",
                row_dend_side = "left",
                column_names_gp = gpar(cex = 2),
                row_names_gp = gpar(cex = 1.6),
                row_dend_width = unit(6,"cm"),
                column_dend_height = unit(3, "cm"),
                # row_dend_gp = gpar(col = "red"),
                km = 4,
                combined_name_fun = NULL,
                row_title = paste0("C", 1:4),
                row_title_gp = gpar(fontsize = 25),
                column_title = "LFQ",
                column_title_gp = gpar(fontsize = 25))

col_fun_2 = circlize::colorRamp2(c(-1.2, 0, 1.2), c("steelblue", "white", "red"))
ht_2 <- Heatmap(dat2_matrix, 
                col = col_fun_2,
                heatmap_legend_param = list(title = "Expression",
                                            title_gp = gpar(fontsize = 25),
                                            labels_gp = gpar(fontsize = 25),
                                            at = c(-1.2, 0, 1.2),
                                            labels = c("low", "zero", "high"),
                                            legend_height = unit(4, "cm")),
                column_names_side = "top", 
                row_names_side = "left",
                row_dend_side = "left",
                column_names_gp = gpar(cex = 2),
                row_names_gp = gpar(cex = 2),
                row_dend_width = unit(6,"cm"),
                column_dend_height = unit(3, "cm"),
                # row_dend_gp = gpar(col = "red"),
                km = 5,
                combined_name_fun = NULL,
                row_title_gp = gpar(fontsize = 25),
                column_title = "SILAC",
                column_title_gp = gpar(fontsize = 25),
                show_heatmap_legend = F,
                show_row_names = F)

ht_list <- ht_1 + ht_2
draw(ht_list)
# export the heatmap plot into .tiff
dev.off()
# K-mean clustering, Manhattan, Ward.D, hclust
