# plot different incorporation
iso1 <- function(files = "ARNAR_Targeted_UDP_NaGlu_isotopologues_All_Cell_Lines", 
                 title = "UDP-GlcNAc",
                 TP = "6h",
                 widetype = "wt2",
                 Met = 8, # For UDP-Glc, Met = 7, IS = 5
                 IS = 5,
                 sep = 10,
                 wt = 30,
                 ht = 18,
                 s = 72){
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Metabolomics_Arnar_Labelling_April2019")
  library(readxl)
  dat <- read_xlsx(paste0(files, ".xlsx"), sheet = "No knockdowns")
  dat <- data.frame(dat)
  dat <- dat[, c(1:6, 11)]
  colnames(dat) <- c("Samples", substr(colnames(dat)[2:7], 1, nchar(colnames(dat)[2:7])-4))
  
  # keep the mean enrichment column
  MeanEnrich <- dat[, c(1, 7)]
  
  # change the name system of the samples in rows
  MeanEnrich[10:27, 1] <- sub("Glc", "Gln", MeanEnrich[10:27, 1], fixed = TRUE) # fix the mistake: should be 1-Gln and 5-Gln NOT Glc
  MeanEnrich[1:27, 1] <- sub("Her2", "D492HER2", MeanEnrich[1:27, 1], fixed = TRUE)
  MeanEnrich[28:57, 1] <- sub("LC", "lc", MeanEnrich[28:57, 1], fixed = TRUE)
  MeanEnrich[28:57, 1] <- sub("LN", "ln", MeanEnrich[28:57, 1], fixed = TRUE)
  MeanEnrich[c(28:30, 43:45), 1] <- sub("-", ",", MeanEnrich[c(28:30, 43:45), 1], fixed = TRUE)
  MeanEnrich[c(28:30, 43:45), 1] <- sub(",2", ",2-", MeanEnrich[c(28:30, 43:45), 1], fixed = TRUE)
  MeanEnrich[, 1] <- c(substr(MeanEnrich[1:27, 1], 1, nchar(MeanEnrich[1:27, 1])-2), MeanEnrich[28:57, 1])
  MeanEnrich[1:27, 1] <- gsub(" ", "_", MeanEnrich[1:27, 1], fixed = TRUE)
  MeanEnrich[ , 1] <- gsub("-", paste0("-", "13C "), MeanEnrich[ , 1], fixed = TRUE)
  MeanEnrich <- MeanEnrich[c(28:42, 43:57, 1:27), ] # change order to D492, D492M and D492HER2

  # calculate mean and sd for plotting, also added columns for indexing
  MeanEnrich[, 1] <- factor(MeanEnrich[, 1], levels = unique(MeanEnrich[, 1]))
  
  # mean
  mes <- data.frame(tapply(MeanEnrich[, 2], MeanEnrich[, 1], mean)) # Means
  colnames(mes) <- "Mean"
  mes$"Samples" <- rownames(mes)
  mes <- mes[, c(2, 1)]
  
  # sd
  stv <- data.frame(tapply(MeanEnrich[, 2], MeanEnrich[, 1], sd)) # SD
  colnames(stv) <- "SD"
  stv$"Samples" <- rownames(stv)
  stv <- stv[, c(2, 1)]
  
  # function to swap D492HER2_1,2Glc_6h and D492HER2_1Gln_6h
  funChg <- function(dat){
    rowna13 <- dat$Samples[13]
    rowna16 <- dat$Samples[16]
    
    dat$Samples[13] <- "a"
    dat$Samples[16] <- "b"
    
    dat$Samples[13] <- rowna16
    dat$Samples[16] <- rowna13
    
    dat <- dat[c(1:12, 16, 14:15, 13, 17:19), ]
    
    return(dat)
  }
  
  mes <- funChg(mes)
  stv <- funChg(stv)
  
  # delete unlabelled samples
  mes <- mes[-c(4:5, 9:10, 12, 15, 18), ]
  stv <- stv[-c(4:5, 9:10, 12, 15, 18), ]
  
  # Merge (mean) and (stv) together
  datSum <- cbind(mes, stv[, 2])
  colnames(datSum)[3] <- "SD"
  datSum$"Cell.Line" <- c(rep("D492", 3), rep("D492M", 3), rep("D492HER2", 6))
  datSum <- datSum[, c(1, 4, 2:3)]
  
  # prepare data for plotting
  datSum_6h <- datSum[c(1:6, 8, 10, 12), ]
  datSum_24h <- datSum[c(1:6, 7, 9, 11), ]
  
  # import the metabolite expression data
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentBasic_03.12.2019/Rfiles")
  source("KDexperiment_Basic_DataOrg.R")
  
  glycanKD <- glycanKD[grep(widetype, glycanKD$Sample.ID, ignore.case = TRUE), ]
  
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
  
  datWT <- GlyFun(glycanKD, IS, Met)
  
  datWT <- rbind(datWT, datWT, datWT)
  
  datWT <- datWT[c(1, 4, 7, 2, 5, 8, 3, 6, 9), ]
  
  funErrorProp <- function(x, y, x.error, y.error){
    
    R <- x*y
    
    # absolute values of the ratios
    R.abs <- abs(R)
    
    r.x <- (x.error/x)^2
    r.y <- (y.error/y)^2
    
    new.error <- R.abs*(r.x + r.y)^(1/2)
    
    return(new.error)
    
  }
  
  datSD <- NULL
  for (i in 1:nrow(datSum_6h)){
    temp <- funErrorProp(datSum_6h[i, 3], datWT[i, 2], datSum_6h[i, 4], datWT[i, 3])
    datSD <- c(datSD, temp)
  }
  
  datSD[which(is.na(datSD))] <- 0
  
  datSD <- data.frame(datSD)
  
  datSum_6h[, 3] <- datSum_6h[, 3]*datWT[, 2]
  datSum_6h[, 4] <- datSD
  
  TP <- tolower(TP)
  if (TP == "6h"){
    mPlot <- datSum_6h
    mPlot[, 1] <- substr(mPlot[, 1], 1, nchar(mPlot[, 1])-3)
    mPlot$"Treatment" <- rep(c("1,2-13C Glc", "1-13C Gln", "5-13C Gln"), 3) 
    mPlot <- mPlot[, c(1:2, 5, 3:4)]
    
  } else if (TP == "24h") {
    mPlot <- datSum_24h
    mPlot[, 1] <- gsub("\\_6.*","", mPlot[, 1])
    mPlot[, 1] <- gsub("\\_24.*","", mPlot[, 1])
    mPlot$"Treatment" <- rep(c("1,2-13C Glc", "1-13C Gln", "5-13C Gln"), 3) 
    mPlot <- mPlot[, c(1:2, 5, 3:4)]
    
  } else {
    cat("please use \"6h\" or \"24h\".")
  }
  
  # function to change the columns into factor for plotting
  funFac <- function(dat){
    dat[, 2] <- factor(dat[, 2], levels = c("D492", "D492M", "D492HER2"))
    dat[, 3] <- factor(dat[, 3], levels = c("1,2-13C Glc", "1-13C Gln", "5-13C Gln"))
    
    return(dat)
  }
  
  mPlot <- funFac(mPlot)

  # plotting
  funPlot <- function(datPlot, sep){
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
    
    # start to plot
    p <- ggplot(datPlot, aes(x = Treatment, y = datPlot[, 4], fill = Cell.Line)) + # linetype = 1, 2...
      geom_bar(stat = "identity", position = position_dodge(), width = 0.9, 
               color = "black", size = 4) +
      geom_errorbar(aes(ymin = datPlot[, 4] - SD, ymax = datPlot[, 4] + SD), 
                    width = 0.2,
                    size = 4,
                    position = position_dodge(0.9)) + 
      scale_fill_manual(values = c("steelblue", "red", "orange1")) + 
      scale_y_continuous(breaks = seq(0, (max(datPlot[, 4])+max(datPlot[, 5])), sep)) +
      theme_classic() +
      labs(title = title, 
           x = "Treatments", 
           y = expression('Relative '^"13"*"C Incorporation")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 120), 
            axis.line = element_line(colour = 'black', size = 4),
            legend.position = "right",
            legend.key.size = unit(2, "cm"),
            legend.text = element_text(size = 64, face = "bold"),
            legend.title = element_text(size = 64, face = "bold"),
            axis.title.y.left = element_text(size = 96, face = "bold",
                                             margin = margin(t = 0, r = 50, b = 0, l = 0)),
            axis.title.x.bottom = element_text(size = 108, face = "bold",
                                               margin = margin(t = 50, r = 0, b = 0, l = 0)),
            axis.text = element_text(size = 82, face = "bold", color = "black"),
            axis.text.x.bottom = element_text(angle = 0, hjust = 0.5),
            axis.ticks.y.left = element_line(colour = "black", size = 4))
    
    setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Metabolomics_Arnar_Labelling_April2019/Figures")
    file <- paste(format(Sys.time(), "%F %H-%M-%S"), "TotalIncorporation", title, "bar.plot.tiff", sep = " ")
    ggsave(filename = file, units = "in", dpi = s, width = wt, height = ht)
    
    print(p)
  }
  
funPlot(mPlot, sep)

}

#------------------------------------------------------------------------#

                 # plot for UDP-GlcNAc and UDP-Glc

#------------------------------------------------------------------------#
iso1()
iso1(files = "ARNAR_Targeted_UDP_Glu_isotopologues_All_Cell_Lines",
     title = "UDP-Glc",
     Met = 7,
     IS = 5,
     sep = 2)

