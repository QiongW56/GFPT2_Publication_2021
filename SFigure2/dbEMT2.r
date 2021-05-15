# referred to the "Meta-core gene list.R" file from C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5 ProteomicPaper\Figures&Tables in the paper\Figures\Figure 1\1A
# Load EMT gene list from literature
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/dbEMT")
metaEMT <- read.delim("dbemt2.txt")
metaEMT <- data.frame(metaEMT) # 1184
metaEMT <- data.frame(as.character(metaEMT$GeneSymbol))
colnames(metaEMT) <- "GENE"

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
meta_lfq <- lfq[which(lfq$Gene.Name..GN. %in% metaEMT$GENE), ] # 239
meta_silacEM <- silacEM[which(silacEM$Gene.Name %in% metaEMT$GENE), ] # 306
meta_silacHE <- silacHE[which(silacHE$Gene.Name %in% metaEMT$GENE), ] # 306

# get useful columns, gene name, pvalue, FDR and FCs for LFQ (EM and HE)
meta_lfq <- meta_lfq[, c(3, 14:15, 18, 20:21, 24)] # 239, no duplicates

# get useful columns, gene name, pvalue and FCs for SILAC
meta_silacEM <- meta_silacEM[, c(3, 8, 7)] # 306, 4 duplicated genes
colnames(meta_silacEM)[2] <- "SILAC_EM_p.value"
meta_silacHE <- meta_silacHE[, c(3, 8, 7)] # 306, 4 duplicated genes
colnames(meta_silacHE)[2] <- "SILAC_HE_p.value"

# merge two datasets into one
meta_silac <- merge(meta_silacEM, meta_silacHE) # 314, 12 duplicated genes (4 genes in total)

funDup <- function(dat){
  
  GeneDup <- dat[, 1][which(duplicated(dat[, 1]))]
  
  # A function to find the duplicated genes in a dataset and find the best duplicate (FDR < 0.05) for this gene and delete other duplicates
  funDupGene <- function(dat, dupGene){
    
    dat <- dat[which(dat[, 1] == dupGene), ] # find the duplicated dataset
    # dat <- dat[which(dat[, grep("p", names(dat))[1]] < 0.05), ] # first delete the rows with FDR (EM) >= 0.05
    # dat <- dat[which(dat[, grep("p", names(dat))[2]] < 0.05), ] # then delete the rows with FDR (HE) >= 0.05
    dat <- dat[which(abs(dat[, grep("Log2", names(dat))[1]]) == max(abs(dat[, grep("Log2", names(dat))[1]]))), ]
    dat <- dat[which(abs(dat[, grep("Log2", names(dat))[2]]) == max(abs(dat[, grep("Log2", names(dat))[2]]))), ]
    # second, leave the max abs(fold changes) for both EM and HE (log2)
    
    return(dat)
    
  }
  
  # delete all genes with duplicated genes
  dat_No_dup <- dat[-which(dat[, 1] %in% GeneDup), ]
  
  for (i in GeneDup){
    df_temp <- funDupGene(dat, i)
    dat_No_dup <- rbind(dat_No_dup, df_temp)
  }
  
  # since every cycle, the same gene probably will be added once,so there are duplicated rows to delete
  dat_No_dup <- dat_No_dup[-which(duplicated(dat_No_dup[, 1])), ]
  
  return(dat_No_dup)
  
}

meta_silac <- funDup(meta_silac)

# in both SILAC and LFQ
meta_lfq_s <- meta_lfq[which(meta_lfq$Gene.Name..GN. %in% meta_silac$Gene.Name), ] # 173
meta_silac_l <- meta_silac[which(meta_silac$Gene.Name %in% meta_lfq$Gene.Name..GN.), ] # 169

# Delete the duplicated genes
meta_lfq_s <- meta_lfq_s[-which(duplicated(meta_lfq_s$Gene.Name..GN.)), ] # 169/173

# meta_silac_l <- meta_silac_l[-c(29, 31, 32), ] # there were 4 TPMs and deleted 3 of them
# meta_silac_l <- meta_silac_l[-which(duplicated(meta_silac_l$Gene.Name)), ] # 169/181

# change columns to merge SILAC and LFQ
colnames(meta_lfq_s)[1] <- "Gene.Name"
colnames(meta_lfq_s)[2] <- "LFQ_EM_p.value"
colnames(meta_lfq_s)[5] <- "LFQ_HE_p.value"
colnames(meta_silac_l)[1] <- "Gene.Name"

s_l <- merge(meta_lfq_s, meta_silac_l) # 169

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

silac <- silac[order(silac$Gene.Name), ] # alphatical order

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

dat2[, -1] <- t(apply(dat2[,-1], 1, cal_z_score))
dat2_matrix <- as.matrix(dat2[,-1])

#  max(dat1_matrix)
# 2.171548
# min(dat1_matrix)
# -2.643683

# max(dat2_matrix)
# 1.1547
# min(dat2_matrix)
# -1.15469

library(ComplexHeatmap)
library(grDevices)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/dbEMT/Figures")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "HeatmapPlot.tiff")
tiff(filename = file, units = "in", width = 12, height = 24, res = 300)

library(seriation)
library(dendextend)
col_fun_1 = circlize::colorRamp2(c(-2.2, 0, 2.7), c("steelblue", "white", "red"))
o1 = seriate(dist(t(dat1_matrix)), method = "TSP") # package: seriation
ht_1 <- Heatmap(dat1_matrix, 
                col = col_fun_1,
                cluster_columns = F,
                heatmap_legend_param = list(title = NULL,
                                            title_gp = gpar(fontsize = 25),
                                            labels_gp = gpar(fontsize = 25),
                                            at = c(-2.2, 0, 2.7),
                                            labels = c("low", "zero", "high"),
                                            legend_height = unit(4, "cm")),
                column_names_side = "top", 
                row_names_side = "left",
                row_dend_side = "left",
                column_names_gp = gpar(cex = 2),
                row_names_gp = gpar(cex = 0.5),
                row_dend_width = unit(6,"cm"),
                column_dend_height = unit(3, "cm"),
                # row_dend_gp = gpar(col = "red"),
                km = 4,
                # combined_name_fun = NULL,
                row_title = paste0("C", 1:4),
                row_title_gp = gpar(fontsize = 25),
                column_title = "LFQ",
                column_title_gp = gpar(fontsize = 25))

col_fun_2 = circlize::colorRamp2(c(-1.2, 0, 1.2), c("steelblue", "white", "red"))
o2 = seriate(dist(t(dat2_matrix)), method = "TSP")
ht_2 <- Heatmap(dat2_matrix, 
                col = col_fun_2,
                cluster_columns = F,
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
                km = 4,
                # combined_name_fun = NULL,
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
