funWT <- function(wt, IS, Met, s, title, acy) {
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentRepeat_16.07.2019/QIONG")
  library(readxl)
  glycanKD <- read_excel("AcidicNeg_Knockdown_Metabolomics_3Scr_18.07.2019.xlsx", sheet = 1)
  glycanKD <- data.frame(glycanKD)
  
  glycanKD <- glycanKD[grep(wt, glycanKD$Sample.ID, ignore.case = TRUE), ]
  
  # calculate mean and SD for plotting, also added columns for indexing
  GlyFun <- function(db, IS, Met) {
    dat <- data.frame(db[, Met]/db[, IS]/db[, 2]*10000)
    colnames(dat) <- colnames(db)[Met]
    dat$"ID" <- factor(c(rep("D492", 3), rep("D492M", 3), rep("D492HER2", 3)), 
                       levels = c("D492", "D492M", "D492HER2"))
    dat <- dat[, c(2, 1)]
    
    # function to calculate mean without taking NA into consideration
    funMean <- function(dat){
      mn <- mean(dat, na.rm = TRUE)
      return(mn)
    }
    
    # function to calculate sd without taking NA into consideration
    funSD <- function(dat){
      mn <- sd(dat, na.rm = TRUE)
      return(mn)
    }
    # mean, sd, sem, number
    m <- tapply(dat[, 2], dat[, 1], funMean)
    stv <- tapply(dat[, 2], dat[, 1], funSD)
    # n <- tapply(dat[, 2], dat[, 1], length)
    # sem <- s/sqrt(n)
    
    dat <- data.frame(rbind(m, stv))
    dat <- data.frame(t(dat))
    
    dat$"Cell.Line" <- factor(rownames(dat), levels = rownames(dat))
    
    dat <- dat[, c(3, 1, 2)]
    
    colnames(dat)[2] <- colnames(db)[Met]
    
    # dat[, 2] <- formatC(signif(dat[, 2],digits = 3)) # 3 significant digits
    
    return(dat)
  }
  
  Plot <- GlyFun(glycanKD, IS, Met)
  
  # function for plotting
  ggplot.dat <- function(datPlot, acy){
    library(ggplot2)
    
    # if a number is < 0.5, it will ceiling to 1, the y axis will be less ticks
    if ((max(datPlot[, 2], na.rm = TRUE) + max(datPlot[, 3], na.rm = TRUE)) > 0.5){
      ymax <- ceiling(max(datPlot[, 2], na.rm = TRUE) + max(datPlot[, 3], na.rm = TRUE))
    } else {
      ymax <- max(datPlot[, 2], na.rm = TRUE) + max(datPlot[, 3], na.rm = TRUE)
    }
    
    p <- ggplot(datPlot, aes(x = Cell.Line, y = datPlot[, 2], fill = Cell.Line)) + # linetype = 1, 2...
      geom_bar(stat = "identity", position = position_dodge(), width = 0.6, 
               color = "black", size = 4) +
      geom_errorbar(aes(ymin = datPlot[, 2] - stv, ymax = datPlot[, 2] + stv), 
                    width = 0.2,
                    size = 4,
                    position = position_dodge(0.9)) + 
      scale_y_continuous(breaks = seq(0, ymax, 
                                      ymax/5),
                         labels = scales::number_format(accuracy = acy)) +
      scale_fill_manual(values = c("steelblue", "red", "orange1")) + 
      theme_classic() +
      labs(title = title, 
           x = NULL, 
           y = "Normalized Expression") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 96), 
            axis.line = element_line(colour = 'black', size = 4),
            legend.position = "none",
            legend.key.size = unit(1, "cm"),
            legend.text = element_text(size = 30),
            legend.title = element_text(size = 30),
            axis.title.y.left = element_text(size = 72, face = "bold",
                                             margin = margin(t = 0, r = 50, b = 0, l = 0)),
            axis.text = element_text(size = 72, face = "bold", color = "black"),
            axis.ticks.y.left = element_line(colour = "black", size = 4))
    
    setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentRepeat_16.07.2019/QIONG/Figures")
    file <- paste(format(Sys.time(), "%F %H-%M-%S"), colnames(glycanKD)[Met], colnames(glycanKD)[IS], wt, "bar.plot.tiff", sep = " ")
    ggsave(filename = file, units = "in", dpi = s, width = 18, height = 14)
    
    print(p)
    
  }
  
  ggplot.dat(Plot, acy = acy)
  
}

#------------------------------------------------------------------------------------#

                           # to print out all the names of metabolites

#------------------------------------------------------------------------------------#
funMet <- function(){
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentRepeat_16.07.2019/QIONG")
  library(readxl)
  glycanKD <- read_excel("AcidicNeg_Knockdown_Metabolomics_3Scr_18.07.2019.xlsx", sheet = 1)
  glycanKD <- data.frame(glycanKD)
  
  colna <- colnames(glycanKD)
  
  return(colna)
  
}

#------------------------------------------------------------------------------------#

                # use "Phenylalanine.IS" as internal standard
                # the "wt1" was used in the paper (first sets of wide type)

#------------------------------------------------------------------------------------#
funWT(IS = 3, Met = 85, s = 72, wt = "wt1", title = "UDP-GlcNAc", acy = 0.1) # UDP-GlcNAc
funWT(IS = 3, Met = 86, s = 72, wt = "wt1", title = "UDP-Glc", acy = 0.01) # UDP-Glc
funWT(IS = 3, Met = 87, s = 72, wt = "wt1", title = "UDP-GlcA", acy = 0.01) # UDP-GlcA
funWT(IS = 3, Met = 88, s = 72, wt = "wt1", title = "GlcNAc-P", acy = 0.01) # GlcNAc-P

# second set of wide type, which was repeated after all knockdowns
funWT(IS = 3, Met = 85, s = 72, wt = "wt2", title = "UDP-GlcNAc") # UDP-GlcNAc
funWT(IS = 3, Met = 86, s = 72, wt = "wt2", title = "UDP-Glc") # UDP-Glc
funWT(IS = 3, Met = 87, s = 72, wt = "wt2", title = "UDP-GlcA") # UDP-GlcA
funWT(IS = 3, Met = 88, s = 72, wt = "wt2", title = "GlcNAc-P") # GlcNAc-P
