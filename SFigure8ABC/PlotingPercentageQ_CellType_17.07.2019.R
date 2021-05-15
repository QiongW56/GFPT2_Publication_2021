funPer <- function(Met = 85, 
                   IS = 3, 
                   KD = "GFPT2",
                   CellType = "",
                   wt = 24,
                   ht = 16,
                   s = 72,
                   title = "UDP-GlcNAc"){
  # Load data, the scramble samples were run twice with different results in MS, replaced the first scr1 with the second scr
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentRepeat_16.07.2019/QIONG")
  library(readxl)
  glycanKD <- read_excel("AcidicNeg_Knockdown_Metabolomics_3Scr_18.07.2019.xlsx", sheet = 1)
  glycanKD <- data.frame(glycanKD)
  
  # delete WT samples
  glycanKD <- glycanKD[-c(1:6, 25:30, 49:54), ]
  
  # normalize all KD data
  for (i in 10:ncol(glycanKD)){
    glycanKD[, i] <- glycanKD[, i]/glycanKD[, 2]/glycanKD[, IS]*10000
  }
  
  # Separate GFPT2, GALE, PGM2L1 and UGDH for Scr and KDs
  GFPT2 <- glycanKD[c(1:6, 7:9, 19:24, 25:27, 37:42, 43:45), ]
  GALE <- glycanKD[c(1:6, 10:12, 19:24, 28:30, 37:42, 46:48), ]
  PGM2L1 <- glycanKD[c(1:6, 13:15, 19:24, 31:33, 37:42, 49:51), ]
  UGDH <- glycanKD[c(1:6, 16:18, 19:24, 34:36, 37:42, 52:54), ]
  
  #----------------------------------------------------------#
  # function to calculate mean without taking NA into consideration
  funMean <- function(dat){
    mn <- mean(dat, na.rm = TRUE)
    return(mn)
  }
  
  #----------------------------------------------------------------------#
  # function to calculate the Scr mean values of all metabolites for three cell lines
  funMeanScr <- function(dat = glycanKD){
    datScr <- dat[c(1:6, 19:24, 37:42), ] # scramble samples from three cell lines in triplicates
    datScr$cell.line <- factor(c(rep("D492", 6), rep("D492M", 6), rep("D492HER2", 6)), 
                               levels = c("D492", "D492M", "D492HER2"))
    m <- NULL
    for (i in 10:(ncol(datScr)-1)){ # -1 is to delete the last column "cell.line"
      temp <- tapply(datScr[, i], datScr[, ncol(datScr)], funMean) # based on the "cell.line" column, calculate mean for all replicates for all metabolites
      m <- cbind(m, temp)
    }
    m <- data.frame(m)
    colnames(m) <- colnames(datScr)[10:(ncol(datScr)-1)]
    
    return(m)
    
  }
  
  ScrMean <- funMeanScr()
  
  # function to calculate Ratio of Scr/Scr(mean) for all cell lines and for all metabolites
  funRatioScr <- function(dat = glycanKD){
    datScr <- dat[c(1:6, 19:24, 37:42), ]
    for (i in 1:nrow(datScr)){
      if (i %in% 1:6){
        datScr[i, 10:ncol(datScr)] <- datScr[i, 10:ncol(datScr)]/ScrMean[1, ] # D492
      } else if (i %in% 7:12){
        datScr[i, 10:ncol(datScr)] <- datScr[i, 10:ncol(datScr)]/ScrMean[2, ] # D492M
      } else {
        datScr[i, 10:ncol(datScr)] <- datScr[i, 10:ncol(datScr)]/ScrMean[3, ] # D492HER2
      }
    }
    return(datScr)
  }
  
  RatioScr <- funRatioScr()
  
  #----------------------------------------------------------------------------#
  # Function to calculate mean of ratios of Scr/Scr(mean) for all metabolites in all cell lines
  funMeanScr <- function(dat){
    dat$cell.line <- factor(c(rep("D492", 6), rep("D492M", 6), rep("D492HER2", 6)), 
                              levels = c("D492", "D492M", "D492HER2"))
    m <- NULL
    for (i in 10:(ncol(dat)-1)){
      temp <- tapply(dat[, i], dat[, ncol(dat)], funMean)
      m <- cbind(m, temp)
    }
    m <- data.frame(m)
    colnames(m) <- colnames(dat)[10:(ncol(dat)-1)]
    
    return(m)
    
  }
  
  RatioScrMean <- funMeanScr(RatioScr)
  
  # Function to calculate Stv or SD of ratios of Scr/Scr(mean) for all metabolites in all cell lines
  funStvScr <- function(dat){
    dat$cell.line <- factor(c(rep("D492", 6), rep("D492M", 6), rep("D492HER2", 6)), 
                               levels = c("D492", "D492M", "D492HER2"))
    stv <- NULL
    for (i in 10:(ncol(dat)-1)){
      temp <- tapply(dat[, i], dat[, ncol(dat)], sd)
      stv <- cbind(stv, temp)
    }
    stv <- data.frame(stv)
    colnames(stv) <- colnames(dat)[10:(ncol(dat)-1)]
    
    return(stv)
    
  }
  
  RatioScrStv <- funStvScr(RatioScr)
  
  #-----------------------------------------------------------------------#
  # function to calculate ratios of KD/Scr for all cell lines for each knockdown of genes
  funRatioKD <- function(dat){
    datKD <- dat[c(7:9, 16:18, 25:27), ]
    for (i in 1:nrow(datKD)){
      if (i %in% 1:3){
        datKD[i, 10:ncol(datKD)] <- datKD[i, 10:ncol(datKD)]/ScrMean[1, ] # D492
      } else if (i %in% 4:6){
        datKD[i, 10:ncol(datKD)] <- datKD[i, 10:ncol(datKD)]/ScrMean[2, ] # D492M
      } else {
        datKD[i, 10:ncol(datKD)] <- datKD[i, 10:ncol(datKD)]/ScrMean[3, ] # D492HER2
      }
    }
    return(datKD)
  }
  
  RatioGFPT2 <- funRatioKD(GFPT2)
  RatioGALE <- funRatioKD(GALE)
  RatioPGM2L1 <- funRatioKD(PGM2L1)
  RatioUGDH <- funRatioKD(UGDH)
  
  # Function to calculate mean of ratios for all metabolites in all cell lines for each KD
  funMeanKD <- function(datKD){
    datKD$cell.line <- factor(c(rep("D492", 3), rep("D492M", 3), rep("D492HER2", 3)), 
                              levels = c("D492", "D492M", "D492HER2"))
    m <- NULL
    for (i in 10:(ncol(datKD)-1)){
      temp <- tapply(datKD[, i], datKD[, ncol(datKD)], funMean)
      m <- cbind(m, temp)
    }
    m <- data.frame(m)
    colnames(m) <- colnames(datKD)[10:(ncol(datKD)-1)]
    
    return(m)
    
  }
  
  RatioGFPT2mean <- funMeanKD(RatioGFPT2)
  RatioGALEmean <- funMeanKD(RatioGALE)
  RatioPGM2L1mean <- funMeanKD(RatioPGM2L1)
  RatioUGDHmean <- funMeanKD(RatioUGDH)
  
  # Function to calculate STV or SD of ratios for all metabolites in all cell lines for each KD
  funStvKD <- function(datKD){
    datKD$cell.line <- factor(c(rep("D492", 3), rep("D492M", 3), rep("D492HER2", 3)), 
                              levels = c("D492", "D492M", "D492HER2"))
    stv <- NULL
    for (i in 10:(ncol(datKD)-1)){
      temp <- tapply(datKD[, i], datKD[, ncol(datKD)], sd)
      stv <- cbind(stv, temp)
    }
    stv <- data.frame(stv)
    colnames(stv) <- colnames(datKD)[10:(ncol(datKD)-1)]
    
    return(stv)
    
  }
  
  RatioGFPT2stv <- funStvKD(RatioGFPT2)
  RatioGALEstv <- funStvKD(RatioGALE)
  RatioPGM2L1stv <- funStvKD(RatioPGM2L1)
  RatioUGDHstv <- funStvKD(RatioUGDH)
  
  #----------------------------------------------------------------------#
  # plot for specific gene knockdowns
  if (KD == "GFPT2"){
    datPlot.scr <- cbind(RatioScrMean[, Met-9], RatioScrStv[, Met-9]) # minus the "sample ID", "protein con." and all IS columns 
    datPlot.kd <- cbind(RatioGFPT2mean[, Met-9], RatioGFPT2stv[, Met-9])
    datPlot <- rbind(datPlot.scr, datPlot.kd)
  } else if (KD == "GALE"){
    datPlot.scr <- cbind(RatioScrMean[, Met-9], RatioScrStv[, Met-9])
    datPlot.kd <- cbind(RatioGALEmean[, Met-9], RatioGALEstv[, Met-9])
    datPlot <- rbind(datPlot.scr, datPlot.kd)
  } else if (KD == "PGM2L1") {
    datPlot.scr <- cbind(RatioScrMean[, Met-9], RatioScrStv[, Met-9])
    datPlot.kd <- cbind(RatioPGM2L1mean[, Met-9], RatioPGM2L1stv[, Met-9])
    datPlot <- rbind(datPlot.scr, datPlot.kd)
  } else {
    datPlot.scr <- cbind(RatioScrMean[, Met-9], RatioScrStv[, Met-9])
    datPlot.kd <- cbind(RatioUGDHmean[, Met-9], RatioUGDHstv[, Met-9])
    datPlot <- rbind(datPlot.scr, datPlot.kd)
  }
  
  # get the data ready for plotting
  datPlot <- data.frame(datPlot)
  colnames(datPlot) <- c(colnames(glycanKD)[Met], "stv")
  datPlot$cell.line <- factor(rep(c("D492", "D492M", "D492HER2"), 2), levels = c("D492", "D492M", "D492HER2"))
  datPlot$Treatment <- factor(c(rep("Scr", 3), rep("KD", 3)), levels = c("Scr", "KD"))
  datPlot$Samples <- paste0(datPlot$Treatment, "_", datPlot$cell.line)
  datPlot <- datPlot[, c(5, 4, 3, 1:2)] # columns: Samples, treatment, cell.line, mean, sd
  datPlot <- datPlot[c(1, 4, 2, 5, 3, 6), ] # rows: D492_Scr, D492_KD,  D492M_Scr, D492M_KD,  D492HER2_Scr, D492HER2_KD
  datPlot$Samples <- factor(datPlot$Samples, levels = datPlot$Samples)
  
  if (CellType == "D492"){
    datPlot <- datPlot[-3:-6, ]
    cols <- c("steelblue", "skyblue1")
    x.ax <- 2 # x axis is treatment (Scr vs. KD)
  } else if (CellType == "D492M"){
    datPlot <- datPlot[c(-1:-2, -5:-6), ]
    cols <- c("red", "pink")
    x.ax <- 2
  } else if (CellType == "D492HER2"){
    datPlot <- datPlot[-1:-4, ]
    cols <- c("orange1", "#FAE48BFF")
    x.ax <- 2
  } else {
    datPlot <- datPlot
    cols <- c("steelblue", "skyblue1", "red", "pink", "orange1", "#FAE48BFF")
    x.ax <- 3 # x axis is cell line (D492, D492M and D492HER2)
  }
  
  #------------------------------------------------------------------------------#
  # plotting
  ggplot.dat <- function(datPlot, wt, ht, s, title){
    library(ggplot2)
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
    
    # start plotting
    p <- ggplot(datPlot, aes(x = datPlot[, x.ax], y = datPlot[, 4], fill = Samples)) + # linetype = 1, 2...
      geom_bar(stat = "identity", position = position_dodge(), color = "black", 
               width = 0.8, size = 4) +
      geom_errorbar(aes(ymin = datPlot[, 4] - stv, ymax = datPlot[, 4] + stv), 
                    width = 0.2,
                    size = 4,
                    position = position_dodge(0.8)) + 
      scale_fill_manual(values = cols) + 
      scale_y_continuous(breaks = seq(0, 1.5, 0.25)) +
      theme_classic() +
      labs(title = title, 
           x = paste0("si", KD), 
           y = "KD/Scr Ratio") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 120), 
            axis.line = element_line(colour = "black", size = 4),
            legend.position = "none",
            legend.key.size = unit(1.5, "cm"),
            legend.text = element_text(size = 45, face = "bold"),
            legend.title = element_text(size = 60, face = "bold"),
            axis.title.y.left = element_text(size = 96, face = "bold", 
                                             margin = margin(t = 0, r = 50, b = 0, l = 0)),
            axis.title.x.bottom = element_text(size = 96, face = "bold",
                                             margin = margin(t = 50, r = 0, b = 0, l = 0)),
            axis.text = element_text(size = 96, face = "bold", color = "black"),
            axis.ticks.y.left = element_line(colour = "black", size = 4)) +
      geom_hline(yintercept = 1, linetype = "dashed", size = 2)
    
    print(p)
    
    setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentRepeat_16.07.2019/QIONG/Figures")
    file <- paste(format(Sys.time(), "%F %H-%M-%S"), colnames(glycanKD)[Met], colnames(glycanKD)[IS], "bar.plot.tiff", sep = " ")
    ggsave(filename = file, units = "in", dpi = s, width = wt, height = ht)
    
    print(p)
  }
  
  ggplot.dat(datPlot, wt = wt, ht = ht, s = s, title = title)
  
}

#------------------------------------------------------------------------------------#

                     # to print out numbers for each metabolite 

#------------------------------------------------------------------------------------#
# to print out all the names of metabolites
funMet <- function(){
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentRepeat_16.07.2019/QIONG")
  library(readxl)
  glycanKD <- read_excel("AcidicNeg_Knockdown_Metabolomics_3Scr_18.07.2019.xlsx", sheet = 1)
  glycanKD <- data.frame(glycanKD)
  
  colna <- colnames(glycanKD)
  
  return(colna)
}

# use "Phenylalanine.IS" as internal standard

# NOTE: first print out figures without legends, 
# then print out legends at the bottom and put them together in PPT.
# col2rgb(), col2hex()
