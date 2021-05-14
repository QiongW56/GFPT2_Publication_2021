# Volcano plot for LFQ data (proteomics)
# prepare the data
# import data from Perseus, since it has -log(p.value)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/VolcanoPlots/VolcanoPlotForLFQ")
dat <- read.delim("PerseusOutPutTable_D492M vs. D492_HER2.txt")
dat <- data.frame(dat)
dat <- dat[-1:-3, ]

for (i in 1:9){
  dat[, i] <- as.numeric(levels(dat[, i]))[dat[, i]]
}

for (i in 1:nrow(dat)){
  dat$Mean_D492[i] <- mean(as.numeric(dat[i, 1:3]), na.rm = TRUE)
  dat$Mean_D492M[i] <- mean(as.numeric(dat[i, 4:6]), na.rm = TRUE)
  dat$Mean_D492HER2[i] <- mean(as.numeric(dat[i, 7:9]), na.rm = TRUE)
}

dat$Log2_D492HER2.D492M <- dat$Mean_D492HER2 - dat$Mean_D492M

for (i in 16:18) {
  dat[, i] <- as.character.factor(dat[ ,i]) # "protein name", "protein ID" and "gene name"
}

dat <- dat[, c(18, 22, 12, 13)] # 2705 "gene name", "fold changes", "p value" and "FDR"

for (i in 3:4) {
  dat[, i] <- as.numeric(levels(dat[ ,i]))[dat[ ,i]] # NAs introduced by coercion
}

# dat <- dat[-grep("ERBB", dat$Gene.Name), ]

#-------------------------------------------------------------------------------------#

# add labels to the volcano plot

#-------------------------------------------------------------------------------------#
# setwd("F:/DundeeAnalysisResults/All datasets_organized")
# met <- read.delim("MetabolicTargetsFromProteomicDataAnalysis.txt", col.names = FALSE)
# met <- unique(met) # 29
# names(met) <- "GeneName"

# erbb <- data.frame(factor("ERBB2", levels = "ERBB2"))
# names(erbb) <- "GeneName"

# met <- rbind(met, erbb)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper")
library(readxl)
df_met <- data.frame(read_xlsx("Tables_1.xlsx", sheet = 4))
met <- data.frame(df_met$Gene.Name) # 5

# rm PXDN, add it manually later because it overlaps with other gene name
met <- data.frame(met[-which(as.character.factor(met[, 1]) == "PXDN"), ]) # 4

# add some genes with fold changes more than 4 or less than -4
dat_met <- data.frame(dat[which(abs(dat[, 2]) >= 4), 1]) # 4

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

table(dat$label) # 5-1+4/2705

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
                         col = c("grey","orange1","red","darkgrey"),
                         title = "",
                         labNames = "D492HER2/D492M",
                         height = 10,
                         width = 18) {
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(grDevices)
  
  # prepare datasets for plotting
  if (fdr == T) {
    dat <- data.frame(dat)
    for (i in 2:4) {
      dat[, i] <- as.numeric(dat[, i])
    }
    colnames(dat)[1:4] <- c("GeneName",
                            "Log2FoldChanges",
                            "pvalue",
                            "FDR")
    if (label == T) {
      colnames(dat)[5] <- "label"
    } else {
      dat <- dat[, -5]
    }
    
  } else {
    dat <- data.frame(dat[, c(1:3, 5)])
    for (i in 2:3) {
      dat[, i] <- as.numeric(dat[, i])
    }
    colnames(dat)[1:3] <- c("GeneName",
                            "Log2FoldChanges",
                            "pvalue")
    if (label == T) {
      colnames(dat)[4] <- "label"
    } else {
      dat <- dat[, -4]
    }
  }
  
  # add another "LFQ.Significance" char column
  if (fdr == TRUE) {
    dat <- dat %>%
      mutate(LFQ.Significance = factor(case_when((Log2FoldChanges < -Log2Fold) & (FDR < significance) ~ "Higher in D492M",
                                          (Log2FoldChanges > Log2Fold) & (FDR < significance) ~ "Higher in D492HER2",
                                          (FDR > significance) ~ "FDR > 0.05",
                                          TRUE ~ "Less than 2-fold")))
  } else {
    dat <- dat %>%
      mutate(LFQ.Significance = factor(case_when((Log2FoldChanges < -Log2Fold) & (pvalue > significance) ~ "downregulated",
                                          (Log2FoldChanges > Log2Fold) & (pvalue > significance) ~ "upregulated",
                                          (pvalue < significance) ~ "not significant",
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
                     color = LFQ.Significance), alpha = 0.5, size = 8) + # alpha for transparency
      scale_color_manual(values = col) +
      scale_y_continuous(breaks = seq(0, 7.5, 1.5)) +
      scale_x_continuous(breaks = seq(-10, 10, 2)) +
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
    
    setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/VolcanoPlots/VolcanoPlotForLFQ")
    file <- paste(format(Sys.time(), "%F %H-%M-%S"), "LFQ_HM_VolcanoPlot.tiff")
    ggsave(filename = file, units="in", width = width, height = height, dpi = 300)
    
  } else {
    p <- ggplot(dat) +
      geom_point(aes(x = Log2FoldChanges, 
                     y = pvalue, 
                     color = LFQ.Significance), alpha = 0.5, size = 8) + # alpha for transparency
      scale_color_manual(values = col) +
      scale_y_continuous(breaks = seq(0, 7.5, 1.5)) +
      scale_x_continuous(breaks = seq(-10, 10, 2)) +
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
    
    setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/VolcanoPlots/VolcanoPlotForLFQ")
    file <- paste(format(Sys.time(), "%F %H-%M-%S"), "LFQ_HM_VolcanoPlot.tiff")
    ggsave(filename = file, units="in", width = width, height = height, dpi = 300)
  }
}

volcano.plot(dat = dat, 
             # title = "LFQ Proteome of EMT Cell Model", 
             yintercept = 1.5)


#-------------------------------------------------------------------------------------#

# some statistical analysis for section "Statistical analysis of proteome in non-malignant and HER2-induced malignant differentiation in D492"

#-------------------------------------------------------------------------------------#
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/VolcanoPlots/VolcanoPlotForLFQ")
dat <- read.delim("PerseusOutPutTable_D492 vs. D492M.txt")
dat <- data.frame(dat)
dat <- dat[-1:-3, ]

for (i in 21:23) {
  dat[, i] <- as.character.factor(dat[ ,i])
}

dat <- dat[, c(23, 15, 17, 18)] # 2705

for (i in 2:4) {
  dat[, i] <- as.numeric(levels(dat[ ,i]))[dat[ ,i]]
}

# number of proteins with significance (FDR < 0.05)
nrow(dat[which(dat$Student.s.T.test.q.value.D492_D492M < 0.05), ]) # 1205

# number of proteins with significance (FDR < 0.05) and at least 2-fold changes
nrow(dat[which(dat$Student.s.T.test.q.value.D492_D492M < 0.05 & dat$Log2_D492.D492M < -1), ])
# 130 downregulated 2fold
nrow(dat[which(dat$Student.s.T.test.q.value.D492_D492M < 0.05 & dat$Log2_D492.D492M > 1), ])
# 124 upregulated 2fold
# 254 in total
# note, 124 downregulated and 130 upregulated in EMT
