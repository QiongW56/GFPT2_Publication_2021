# Volcano plot for SILAC data (proteomics)
# prepare the data
# import data from Perseus, since it has -log(p.value)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/MatchingLFQ&SILAC")
dat <- read.delim("PerseusOutputTable_SILAC_HE.txt") # already a dataframe

for (i in 11:13) {
  dat[, i] <- as.character.factor(dat[ ,i]) # only need to change to character, already are numeric
}

dat$"Median_Log2.E.H." <- -dat$Median_Log2.H.E. # D492/D492HER2
dat <- dat[, c(13, 14, 8, 9, 1:4, 6:7, 10:12)] # 5119

#-------------------------------------------------------------------------------------#

# add labels to the volcano plot

#-------------------------------------------------------------------------------------#
# setwd("F:/DundeeAnalysisResults/All datasets_organized")
# met <- read.delim("MetabolicTargetsFromProteomicDataAnalysis.txt", col.names = FALSE)
# met <- unique(met)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper")
library(readxl)
df_met <- data.frame(read_xlsx("Tables_1.xlsx", sheet = 3))
met <- data.frame(df_met$Gene.Name) # 12

# add some genes with fold changes more than 4 or less than -4
dat_met <- data.frame(dat[which(abs(dat[, 2]) >= 4), 1]) # 11

# rm NA from the gene list
dat_met <- data.frame(dat_met[which(!is.na(dat_met$dat.which.abs.dat...2......4...1.)), ]) # 10

colnames(met) <- "FALSE."
colnames(dat_met) <- "FALSE."

# met <- rbind(met, dat_met)

for (i in 1:nrow(dat)) {
  if (dat$Gene.Name[i] %in% met$FALSE.) {
    dat$"label"[i] <- "YES"
  }  else if (dat$Gene.Name[i] %in% dat_met$FALSE.){
    dat$"label"[i] <- "FC"
  } else {
    dat$"label"[i] <- "NO"
  }
}

table(dat$label) # 12+10/5119

dat <- dat[, c(1:4, 14, 5:13)]

#-------------------------------------------------------------------------------------#

# plotting

#-------------------------------------------------------------------------------------#
volcano.plot <- function(dat, 
                         fdr = T,
                         label = T,
                         Log2Fold = 1, 
                         significance = 0.05, 
                         xintercept = 1, 
                         yintercept = 1.3,
                         col = c("grey","steelblue","orange1","darkgrey"),
                         title = "",
                         labNames = "D492/D492HER2",
                         height = 10,
                         width = 18) {
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(grDevices)
  
  # prepare datasets for plotting
  if (fdr == TRUE) {
    dat <- data.frame(dat)
    for (i in 2:4) {
      dat[, i] <- as.numeric(dat[, i])
    }
    colnames(dat)[1:4] <- c("GeneName",
                            "Log2FoldChanges",
                            "pvalue",
                            "FDR")
    if (label == TRUE) {
      colnames(dat)[5] <- "label"
    } else {
      dat <- dat[, -5]
    }
    
  } else {
    dat <- data.frame(dat[, c(1:3, 5)]) # delete FDR
    for (i in 2:3) {
      dat[, i] <- as.numeric(dat[, i])
    }
    colnames(dat)[1:3] <- c("GeneName",
                            "Log2FoldChanges",
                            "pvalue")
    if (label == TRUE) {
      colnames(dat)[4] <- "label"
    } else {
      dat <- dat[, -4]
    }
  }
  
  # add another "SILAC.Significance" char column
  if (fdr == TRUE) {
    dat <- dat %>%
      mutate(SILAC.significance = factor(case_when((Log2FoldChanges < -Log2Fold) & (FDR < significance) ~ "Higher in D492HER2",
                                          (Log2FoldChanges > Log2Fold) & (FDR < significance) ~ "Higher in D492",
                                          (FDR > significance) ~ "FDR > 0.05",
                                          TRUE ~ "Less than 2-fold")))
  } else {
    dat <- dat %>%
      mutate(SILAC.significance = factor(case_when((Log2FoldChanges < -Log2Fold) & (pvalue > significance) ~ "downregulated",
                                          (Log2FoldChanges > Log2Fold) & (pvalue > significance) ~ "upregulated",
                                          (pvalue > significance) ~ "not significant",
                                          TRUE ~ "significant with less fold changes")))
  }
  
  # function to change the space between keys in legends
  draw_key_polygon3 <- function(data, params, size) {
    lwd <- min(data$size, min(size) / 4)
    
    grid::rectGrob(
      width = grid::unit(0.6, "npc"),
      height = grid::unit(0.6, "npc"),
      gp = grid::gpar(
        col = data$colour,
        fill = alpha(data$fill, data$alpha),
        lty = data$linetype,
        lwd = lwd * .pt,
        linejoin = "mitre"
      ))
  }
  GeomBar$draw_key = draw_key_polygon3
  
  # Plotting
  if (TRUE %in% grepl("label", colnames(dat), ignore.case = F)) {
    p <- ggplot(dat) +
      geom_point(aes(x = Log2FoldChanges, 
                     y = pvalue, 
                     color = SILAC.significance), alpha = 0.5, size = 8) + # alpha for transparency
      scale_color_manual(values = col) +
      scale_y_continuous(breaks = seq(0, 7.5, 1.5)) +
      scale_x_continuous(breaks = seq(-8, 8, 2)) +
      ggtitle(title) +
      xlab(bquote("Log"[2]*"("*.(labNames)*")")) +
      ylab(expression("-Log"[10]*"(p value)")) +
      theme_classic() +
      theme(legend.position = "right", 
            legend.title = element_text(size = rel(3), face = "bold"),
            legend.text = element_text(size = rel(2.8), face = "bold"),
            plot.title = element_text(size = rel(3.6), face = "bold", hjust = 0.5),
            axis.title = element_text(size = rel(4.8), face = "bold"),
            axis.text = element_text(size = rel(4.6), face = "bold", color = "black"),
            axis.line = element_line(size = 4)) +
      geom_hline(yintercept = yintercept, linetype = "dashed", size = 2) +
      geom_vline(xintercept = xintercept, linetype = "dashed", size = 2) +
      geom_vline(xintercept = -xintercept, linetype = "dashed", size = 2) +
      geom_text_repel(aes(x = Log2FoldChanges, y = pvalue, 
                          label = ifelse(label == "YES", GeneName,"")), size = 8) +
      geom_text_repel(aes(x = Log2FoldChanges, y = pvalue, 
                          label = ifelse(label == "FC", GeneName,"")), size = 8, fontface = "bold")
    print(p)
    
    setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/VolcanoPlots/VolcanoPlotForSILAC")
    file <- paste(format(Sys.time(), "%F %H-%M-%S"), "SILAC_HER2_VolcanoPlot.tiff")
    ggsave(filename = file, units="in", width = width, height = height)
    
  } else {
    p <- ggplot(dat) +
      geom_point(aes(x = Log2FoldChanges, 
                     y = pvalue, 
                     color = SILAC.significance), alpha = 0.5, size = 8) + # alpha for transparency
      scale_color_manual(values = col) +
      scale_y_continuous(breaks = seq(0, 7.5, 1.5)) +
      scale_x_continuous(breaks = seq(-8, 8, 2)) +
      ggtitle(title) +
      xlab(bquote("Log"[2]*"("*.(labNames)*")")) +
      ylab(expression("-Log"[10]*"(p value)")) +
      theme_classic() +
      theme(legend.position = "right", 
            legend.title = element_text(size = rel(3), face = "bold"),
            legend.text = element_text(size = rel(2.8), face = "bold"),
            plot.title = element_text(size = rel(3.6), face = "bold", hjust = 0.5),
            axis.title = element_text(size = rel(4.8), face = "bold"),
            axis.text = element_text(size = rel(4.8), face = "bold", color = "black"),
            axis.line = element_line(size = 4)) +
      geom_hline(yintercept = yintercept, linetype = "dashed", size = 2) +
      geom_vline(xintercept = xintercept, linetype = "dashed", size = 2) +
      geom_vline(xintercept = -xintercept, linetype = "dashed", size = 2)
    
    print(p)
    
    setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/VolcanoPlots/VolcanoPlotForSILAC")
    file <- paste(format(Sys.time(), "%F %H-%M-%S"), "SILAC_HER2_VolcanoPlot.tiff")
    ggsave(filename = file, units="in", width = width, height = height)
  }
}

volcano.plot(dat = dat, 
             # title = " SILAC Proteome of HER2-induced Tumorigenesis", 
             yintercept = 1.5)


#-------------------------------------------------------------------------------------#

# statistical analysis for section "Statistical analysis of proteome in non-malignant and HER2-induced malignant differentiation in D492"

#-------------------------------------------------------------------------------------#
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Rfiles_MatchingLFQ&SILAC")
dat <- read.delim("PerseusOutputTable_SILAC_HE.txt")
for (i in 11:13) {
  dat[, i] <- as.character.factor(dat[ ,i])
}
dat$"Median_Log2.E.H." <- -dat$Median_Log2.H.E. # D492/D492HER2
dat <- dat[, c(13, 14, 8, 9, 1:4, 6:7, 10:12)] # 5119

# number of proteins with significance (p value < 0.05)
nrow(dat[which(dat$p.value < 0.05), ]) # 2981

# number of proteins with significance (p value < 0.05)
nrow(dat[which(dat$p.value < 0.05 & dat$Median_Log2.E.H. > 1), ]) # 228 upregulated
nrow(dat[which(dat$p.value < 0.05 & dat$Median_Log2.E.H. < -1), ]) # 221 downregulated
# note, 228 downregulated and 221 upregulated in HER2 overexpression
