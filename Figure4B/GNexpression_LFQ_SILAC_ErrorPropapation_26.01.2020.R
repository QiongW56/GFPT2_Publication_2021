funGN <- function(GN){
  #-----------------------------------------------------------------------------------------#
  
                                           # LFQ
  
  #-----------------------------------------------------------------------------------------#
  # import the LFQ proteomic dataset from Dundee
  setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018")
  library(readxl)
  lfq_dundee <- read_excel("Proteomics_LFQ_summary_16.04.2018.xlsx", sheet = 1)
  lfq_dundee <- data.frame(lfq_dundee)
  lfq_dundee <- lfq_dundee[, c(5, 8:19)]
  
  # Find the expression of a specific gene
  lfq_dundee <- lfq_dundee[which(lfq_dundee[ ,1] == GN), ]
  
  # Delete one row if there are more than one rows
  lfq_dundee <- lfq_dundee[1, ]
  
  # format the dataframe
  lfq_dundee <- data.frame(t(lfq_dundee))
  lfq_dundee <- data.frame(lfq_dundee[-1, ])
  colnames(lfq_dundee) <- GN
  
  lfq_dundee$"Cell.Line" <- factor(rep(c("D492", "D492M", "D492HER2", "D492DEE"), times = c(3,3,3,3)), 
                                   levels = c("D492", "D492M", "D492HER2", "D492DEE"))
  lfq_dundee$"Quantification" <- factor(rep("LFQ", 12), levels = "LFQ")
  lfq_dundee <- lfq_dundee[, c(2:3, 1)]
  lfq_dundee[, 3] <- as.numeric(lfq_dundee[, 3])
  lfq_dundee[, 3] <- lfq_dundee[, 3]/1000000
  
  # function to calculate mean and sem
  mean.sd <- function(dat){
    # mean, sd, sem, number
    colna <- colnames(dat)[3]
    
    # since there are "NA" in the dataset, mean() will return "NA"
    funMean <- function(dat){
      
      m <- mean(dat, na.rm = TRUE)
      
      return(m)
    }
    
    funSD <- function(dat){
      
      s <- sd(dat, na.rm = TRUE)
      
      return(s)
    }
    
    m <- tapply(dat[, 3], dat[, 1], funMean)
    
    s <- tapply(dat[, 3], dat[, 1], funSD)
    # n <- tapply(dat[, 2], dat[, 1], length)
    # sem <- s/sqrt(n)
    
    dat <- data.frame(rbind(m, s))
    dat <- data.frame(t(dat))
    dat$"Cell.Line" <- factor(rownames(dat), levels = rownames(dat))
    dat$"Quantification" <- factor(rep("LFQ", 4), levels = "LFQ")
    
    colnames(dat)[1] <- colna
    colnames(dat)[2] <- "SD"
    
    dat <- dat[, c(3:4, 1:2)]
    
    return(dat)
    
  }
  
  LFQ.plot <- mean.sd(lfq_dundee)
  
  #----------------------------------------------------------------------------------------#
  
  # SILAC
  
  #----------------------------------------------------------------------------------------#
  setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018")
  library(readxl)
  silacEM <- read_excel("Proteomics_SILAC_summary_score_28.12.2018.xlsx", sheet = 1)
  silacEM <- data.frame(silacEM)
  silacEM <- silacEM[, c(3, 5, 9, 13)]
  silacEM[, 2:4] <- 1/silacEM[, 2:4] # D492M/D492
  colnames(silacEM)[2:4] <- paste0("D492M", "_", 1:3)
  silacEM.GN <- silacEM[which(silacEM[, 1] == GN), ]
  
  silacHE <- read_excel("Proteomics_SILAC_summary_score_28.12.2018.xlsx", sheet = 3)
  silacHE <- data.frame(silacHE)
  silacHE <- silacHE[, c(3, 5, 9, 13)]
  colnames(silacHE)[2:4] <- paste0("D492HER2", "_", 1:3)
  silacHE.GN <- silacHE[which(silacHE[, 1] == GN), ]
  
  # function when there are more than one rows for one gene
  funCombine <- function(dat){
    
    if (nrow(dat) >= 2){
      datCom <- 0
      for (i in 2:ncol(dat)){
        temp <- dat[, i][!is.na(dat[, i])]
        datCom <- c(datCom, temp)
      }
      
      colna <- colnames(dat)
      
      mtrx <- data.frame(matrix(1:4, nrow = 1, ncol = 4))
      
      colnames(mtrx) <- colna
      
      mtrx[1, ] <- datCom
      
      mtrx[1, 1] <- GN
      
    } else {
      mtrx <- dat
    }
    
    return(mtrx)
    
  }
  
  silacEM.GN <- funCombine(silacEM.GN)
  silacHE.GN <- funCombine(silacHE.GN)
  
  silacGN <- merge(silacEM.GN, silacHE.GN)
  silacGN[, 8:10] <- c(1, 1, 1)
  colnames(silacGN)[8:10] <- paste0("D492", "_", 1:3)
  silacGN <- silacGN[, c(1, 8:10, 2:7)]
  silacGN <- data.frame(t(silacGN))
  silacGN <- data.frame(silacGN[-1, ])
  colnames(silacGN) <- GN
  silacGN[, 1] <- as.numeric(silacGN[, 1])
  silacGN$"Cell.Line" <- factor(c(rep("D492", 3),
                                  rep("D492M", 3),
                                  rep("D492HER2", 3)), levels = c("D492", "D492M", "D492HER2"))
  silacGN$"Quantification" <- factor(rep("SILAC", 3), levels = "SILAC")
  silacGN <- silacGN[, c(2:3, 1)]
  
  # function to check if there is only "NA" for this specific gene, then skip the dataset
  funNA <- function(dat){
    if (sum(is.na(dat[3, ])) == length(dat[3, ])) {
      print("This gene was not detected")
    } else {
      dat <- dat
    }
    return(dat)
  }
  
  # function to calculate mean and sem
  mean.sd <- function(dat){
    # mean, sd, sem, number
    colna <- colnames(dat)[3]
    
    # since there are "NA" in the dataset, mean() will return "NA"
    if ((sum(!is.na(dat[4:6, 3])) == 3 & sum(!is.na(dat[7:9, 3])) == 3) | (sum(is.na(dat[4:6, 3])) == 3 & sum(is.na(dat[7:9, 3])) == 3)){
      m <- tapply(dat[, 3], dat[, 1], mean)
    } else if (sum(!is.na(dat[4:6, 3])) != 3 | sum(!is.na(dat[7:9, 3])) != 3){
      m1 <- data.frame(sum(dat[1:3, 3][which(!is.na(dat[1:3, 3]))])/sum(!is.na(dat[1:3, 3])))
      m2 <- data.frame(sum(dat[4:6, 3][which(!is.na(dat[4:6, 3]))])/sum(!is.na(dat[4:6, 3])))
      m3 <- data.frame(sum(dat[7:9, 3][which(!is.na(dat[7:9, 3]))])/sum(!is.na(dat[7:9, 3])))
      m <- cbind(m1, m2, m3)
      colnames(m) <- c(levels(dat[, 1]))
    }
    
    # since there are "NA" in the dataset, sd() will return "NA"
    if ((sum(!is.na(dat[4:6, 3])) == 3 & sum(!is.na(dat[7:9, 3])) == 3) | (sum(is.na(dat[4:6, 3])) == 3 & sum(is.na(dat[7:9, 3])) == 3)){
      s <- tapply(dat[, 3], dat[, 1], sd)
    } else if (sum(!is.na(dat[4:6, 3])) != 3 | sum(!is.na(dat[7:9, 3])) != 3){
      s1 <- data.frame(sd(dat[1:3, 3][which(!is.na(dat[1:3, 3]))]))
      s2 <- data.frame(sd(dat[4:6, 3][which(!is.na(dat[4:6, 3]))]))
      s3 <- data.frame(sd(dat[7:9, 3][which(!is.na(dat[7:9, 3]))]))
      s <- cbind(s1, s2, s3)
      colnames(s) <- c(levels(dat[, 1]))
    }
    
    # n <- tapply(dat[, 1], dat[, 2], length)
    # sem <- s/sqrt(n)
    
    dat <- data.frame(rbind(m, s))
    dat <- data.frame(t(dat))
    dat$"Cell.Line" <- factor(rownames(dat), levels = rownames(dat))
    dat$"Quantification" <- factor(rep("SILAC", 3), levels = "SILAC")
    colnames(dat)[1] <- colna
    colnames(dat)[2] <- "SD"
    dat <- dat[, c(3:4, 1:2)]
    
    return(dat)
    
  }
  
  silacGN <- funNA(silacGN)
  silacGN.plot <- mean.sd(silacGN)
  
  #--------------------------------------------------------------------------------------#
  
  # Merge LFQ and SILAC, plotting
  
  #--------------------------------------------------------------------------------------#
  datGN <- rbind(LFQ.plot, silacGN.plot)
  datGN <- datGN[-which(rownames(datGN) == "D492DEE"), ]
  
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Figures/SupplementaryFig6/Rfiles")
  source("FunErrorPropagation.R")
  
  datGN.error <- datGN
  
  if (datGN[1, 3] != 0){
    # LFQ
    datGN[1:3, 3] <- datGN[1:3, 3]/datGN[1, 3]
    
    xyz <- as.numeric(datGN.error[1:3, 3])
    xyz.error <- as.numeric(datGN.error[1:3, 4])
    
    ErrorProp <- funErrorProp(xyz[1], xyz[2], xyz[3], xyz.error[1], xyz.error[2], xyz.error[3])
    
    datGN[1:3, 4] <- as.numeric(ErrorProp)
    
    # SILAC
    datGN[4:6, 3] <- datGN[4:6, 3]/datGN[4, 3]
    
    xyz <- as.numeric(datGN.error[4:6, 3])
    xyz.error <- as.numeric(datGN.error[4:6, 4])
    
    ErrorProp <- funErrorProp(xyz[1], xyz[2], xyz[3], xyz.error[1], xyz.error[2], xyz.error[3])
    
    datGN[4:6, 4] <- as.numeric(ErrorProp)
    
  } else if (datGN[2, 3] != 0) {
    # LFQ
    datGN[1:3, 3] <- datGN[1:3, 3]/datGN[2, 3]
    
    ErrorProp <- funErrorProp(datGN.error[2, 3], datGN.error[1, 3], datGN.error[3, 3], 
                              datGN.error[2, 4], datGN.error[1, 4], datGN.error[3, 4])
    
    datGN[1:3, 4] <- c(as.numeric(ErrorProp)[2], as.numeric(ErrorProp)[1], as.numeric(ErrorProp)[3])
    
    # SILAC
    datGN[4:6, 3] <- datGN[4:6, 3]/datGN[5, 3]
    
    ErrorProp <- funErrorProp(datGN.error[5, 3], datGN.error[4, 3], datGN.error[6, 3], 
                              datGN.error[5, 4], datGN.error[4, 4], datGN.error[6, 4])
    
    datGN[4:6, 4] <- c(as.numeric(ErrorProp)[2], as.numeric(ErrorProp)[1], as.numeric(ErrorProp)[3])
    
  } else if (datGN[3, 3] != 0) {
    # LFQ
    datGN[1:3, 3] <- datGN[1:3, 3]/datGN[3, 3]
    
    ErrorProp <- funErrorProp(datGN.error[3, 3], datGN.error[1, 3], datGN.error[2, 3], 
                              datGN.error[3, 4], datGN.error[1, 4], datGN.error[2, 4])
    
    datGN[1:3, 4] <- c(as.numeric(ErrorProp)[3], as.numeric(ErrorProp)[1], as.numeric(ErrorProp)[2])
    
    # SILAC
    datGN[4:6, 3] <- datGN[4:6, 3]/datGN[6, 3]
    
    ErrorProp <- funErrorProp(datGN.error[6, 3], datGN.error[4, 3], datGN.error[5, 3], 
                              datGN.error[6, 4], datGN.error[4, 4], datGN.error[5, 4])
    
    datGN[4:6, 4] <- c(as.numeric(ErrorProp)[3], as.numeric(ErrorProp)[1], as.numeric(ErrorProp)[2])
    
  } else {
    datGN <- datGN
  }
  
  # change SD from "NaN" to "0"
  datGN$SD[which(is.na(datGN$SD))] <- 0
  
  # Plotting
  ggplot.dat <- function(dat){
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
    
    # plotting
    p <- ggplot(dat, aes(x = Quantification, y = dat[, 3], fill = Cell.Line)) + # linetype = 1, 2...
      geom_bar(stat = "identity", position = position_dodge(), width = 0.8, 
               color = "black", size = 4) +
      geom_errorbar(aes(ymin = dat[, 3] - SD, ymax = dat[, 3] + SD), 
                    width = 0.2,
                    size = 4,
                    position = position_dodge(0.8)) + 
      scale_fill_manual(values = c("steelblue", "red", "orange1")) + 
      scale_y_continuous(breaks = seq(0, max(dat[, 3]), max(dat[, 3]/5)), 
                         labels = scales::number_format(accuracy = 0.01, 
                                                        decimal.mark = ".")) +
      theme_classic() +
      labs(title = GN, 
           x = NULL, 
           y = paste0("Relative Expression", " ", "(", colnames(dat)[3], ")")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 96), 
            axis.line = element_line(colour = 'black', size = 4),
            legend.position = "none",
            legend.key.size = unit(1.3, "cm"),
            legend.text = element_text(size = 60, face = "bold"),
            legend.title = element_text(size = 60, face = "bold"),
            axis.title.y.left = element_text(size = 64, face = "bold",
                                             margin = margin(t = 0, r = 50, b = 0, l = 0)),
            axis.text = element_text(size = 96, face = "bold", color = "black"),
            axis.ticks.y.left = element_line(colour = "black", size = 4)) +
      geom_hline(yintercept = 1, linetype = "dashed", size = 2)
    
    print(p)
    
    setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Figures/SupplementaryFig6")
    file <- paste(format(Sys.time(), "%F %H-%M-%S"), colnames(datGN)[3], "LFQ_SILAC_ErrorProp_bar.plot.tiff", sep = " ")
    ggsave(filename = file, units = "in", dpi = 300, width = 20, height = 14)
    
    print(p)
    
  }
  
  ggplot.dat(dat = datGN)
  
}

# Plotting all the targets
funGN(GN = "GALE")
funGN(GN = "PGM2L1")
funGN(GN = "UGDH")
funGN(GN = "GFPT2")

