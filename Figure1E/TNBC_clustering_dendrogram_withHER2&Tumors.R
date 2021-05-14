# Deleted the scripts on how the datasets used below were generated
# if you want to know, check the other R file "TNBC_clustering_dendrogram.R"

#------------------------------------------------------------------------------------#

                                   # Dendrogram plotting (LFQ)

#------------------------------------------------------------------------------------#
# upload the LFQ data (LFQ and Literature)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/1 DundeeProteomicDataSet_14112017/iBAQ_TNBC")
library(readxl)
lfq <- read_excel("7 2541 Proteins in iBAQ and literature and with all detected_z score.xlsx", sheet = 1)
lfq <- data.frame(lfq)

# test if "7 2541 Proteins in iBAQ and literature and with all detected_z score.csv" is 
# the same as "7 2541 Proteins in iBAQ and literature and with all detected_z score.xlsx", sheet = 1
# lfq1 <- read.csv("7 2541 Proteins in iBAQ and literature and with all detected_z score.csv")
# lfq1 <- data.frame(lfq)
# sum(lfq[, 4:59]) == sum(lfq1[, 4:59])

# rename the HER2 names
colnames(lfq)[54:56] <- paste0("D492HER1", "_", 1:3)

lfq <- lfq[, c(-1:-2, -57:-59, -grep("Tumor", names(lfq)))] # delete DEE

# function for dendrogram plotting
funDenLFQ <- function(file){
  library(ggplot2)
  library(ggdendro)
  library(dendextend)
  
  file.cor <- cor(file[, -1], method = "pearson") # a dataset with samples on both x-axis and y-axis
  d.cor <- as.dist(1 - file.cor) # a list of distance between samples in a triangular shape
  file.hclust <- hclust(d.cor, method = "ward.D2") # a list
  
  dend <- as.dendrogram(file.hclust) # not a plot, but a list, already have the hierarchical structure
  dend_data <- dendro_data(dend, type = "rectangle") # a list
  file.cut <- cutree(file.hclust, k = 4) # to color the labels into 4 groups
  file.cut.df <- data.frame(label = names(file.cut), 
                            cluster = factor(file.cut))
  dend_data[["labels"]] <- merge(dend_data[["labels"]], file.cut.df, by="label")
  
  # add extra rectangles
  # dat <- data.frame(x1 = c(7, 18), x2 = c(8, 19), y1 = c(1, 1), y2 = c(2, 2))
  
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Figures/Figure 1/1C")
  file <- paste(format(Sys.time(), "%F %H-%M-%S"), "LFQ_dendrogram.plot.tiff", sep = " ")
  tiff(filename = file, units = "in", res = 300, width = 8, height = 5)
  
  dendrogram <- ggplot(dend_data$segments) +
    geom_segment(aes(x = x, 
                     y = y, 
                     xend = xend, 
                     yend = yend), 
                     size = 1) + # the boldness of the dendrogram
    geom_text(data = dend_data$labels, 
              aes(x, y, label = label, hjust = 0, color = cluster),
              angle = 270, 
              size = 3.5,
              fontface = "bold") +
    scale_colour_manual(values = c("steelblue", "black", "orange1", "red")) +
    ylim(-1, 3) +
    xlim(1, 57) +
    theme_void() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5)) +
    # geom_segment(aes(x = 14.5, y = -0.7, xend = 17.5, yend = -0.7), size = 1) + # x = starting point, xend = ending point
    # geom_segment(aes(x = 35.5, y = -0.52, xend = 38.2, yend = -0.52), size = 1) +
    # geom_segment(aes(x = 35.5, y = -0.52, xend = 38.2, yend = -0.52), size = 1) +
    ggtitle("LFQ") +
    theme(plot.title = element_text(size = 24, face = "bold")) +
    annotate("text", x = 6, y = -1.0, label = "Luminal", size = 4, fontface =2) +
    annotate("text", x = 17, y = -1.0, label = "Basal-like 2", size = 4, fontface =2) +
    annotate("text", x = 29, y = -1.0, label = "Claudin-low", size = 4, fontface =2) +
    # annotate("text", x = 36.5, y = -1.0, label = "Tumor", size = 4, fontface =2) +
    annotate("text", x = 42.5, y = -1.0, label = "Basal-like 1", size = 4, fontface =2)
  
  print(dendrogram)
  
  dev.off()
  
  print(dendrogram)
}

funDenLFQ(lfq)

#----------------------------------------------------------------------------------#

                                  # SILAC plotting

#----------------------------------------------------------------------------------#
# upload the SILAC data (SILAC and literature)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/1 DundeeProteomicDataSet_14112017/iBAQ_TNBC/iBAQ_SILAC_2")
library(readxl)
silac <- read_excel("SILAC_Literature_TNBC_5385.xlsx", sheet = 1)
silac <- data.frame(silac)
colnames(silac)[54:56] <- paste0("D492HER2", "_", 1:3)
silac <- silac[, c(-1:-2, -57:-59, -grep("Tumor", names(silac)))]

# function for dendrogram plotting
funDenSILAC <- function(file){
  library(ggplot2)
  library(ggdendro)
  library(dendextend)
  
  file.cor <- cor(file[, -1], method = "pearson") # a dataset with samples on both x-axis and y-axis
  d.cor <- as.dist(1 - file.cor) # a list of distance between samples in a triangular shape
  file.hclust <- hclust(d.cor, method = "ward.D2") # a list
  
  dend <- as.dendrogram(file.hclust) # not a plot, but a list, already have the hierarchical structure
  dend_data <- dendro_data(dend, type = "rectangle") # a list
  file.cut <- cutree(file.hclust, k = 4) # to color the labels into 4 groups
  file.cut.df <- data.frame(label = names(file.cut), 
                            cluster = factor(file.cut))
  dend_data[["labels"]] <- merge(dend_data[["labels"]], file.cut.df, by="label")
  
  # add extra rectangles
  # dat <- data.frame(x1 = c(7, 18), x2 = c(8, 19), y1 = c(1, 1), y2 = c(2, 2))
  
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Figures/Figure 1/1C")
  file <- paste(format(Sys.time(), "%F %H-%M-%S"), "SILAC_dendrogram.plot.tiff", sep = " ")
  tiff(filename = file, units = "in", res = 300, width = 8, height = 5)
  
  dendrogram <- ggplot(dend_data$segments) +
    geom_segment(aes(x = x, 
                     y = y, 
                     xend = xend, 
                     yend = yend), 
                     size = 1) + # the boldness of the dendrogram
    geom_text(data = dend_data$labels, 
              aes(x, y, label = label, hjust = 0, color = cluster),
              angle = 270, 
              size = 3.5,
              fontface = "bold") +
    scale_colour_manual(values = c("red", "black", "steelblue", "orange1")) +
    ylim(-1, 3) +
    xlim(1, 57) +
    theme_void() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5)) +
    # geom_segment(aes(x = 10.5, y = -0.55, xend = 13.2, yend = -0.55), size = 1) +
    # geom_segment(aes(x = 33.5, y = -0.45, xend = 36.2, yend = -0.45), size = 1) +
    ggtitle("SILAC") +
    theme(plot.title = element_text(size = 24, face = "bold")) +
    annotate("text", x = 5.5, y = -1, label = "Claudin-low", size = 4, fontface =2) +
    # annotate("text", x = 12.5, y = -1, label = "Tumor", size = 4, fontface =2) +
    annotate("text", x = 18.5, y = -1, label = "Basal-like 1", size = 4, fontface =2) +
    annotate("text", x = 30.5, y = -1, label = "Luminal", size = 4, fontface =2) +
    annotate("text", x = 42.5, y = -1, label = "Basal-like 2", size = 4, fontface =2)
  
  print(dendrogram)
  
  dev.off()
  
  print(dendrogram)
}

funDenSILAC(silac)
